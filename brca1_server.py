#!/usr/bin/env python3
"""Shared store for the BRCA1 variant explorer.

Serves the page and its data release, and keeps the laboratory's records in one SQLite file for everyone:
notes, laboratory classifications (proposed by one reviewer, signed off by a second), snapshots, pedigrees,
the reclassification watch list and in-house data. Every change is written to an audit trail.

Standard library only (Python 3.8+). Run it on a machine on the hospital network; colleagues open
http(s)://<host>:<port>/ in their browser. Without this server the page still works and keeps records in each
browser only.

    python3 brca1_server.py adduser alice --display "Alice Tremblay" --admin   # asks for a password
    python3 brca1_server.py serve --port 8443 --cert cert.pem --key key.pem
    python3 brca1_server.py users | passwd NAME | disable NAME | enable NAME
    python3 brca1_server.py export > backup.json                              # everything, for backups

The database holds patient-related information (pedigrees, sample identifiers). Keep it on hospital storage,
back it up with the laboratory's other records, and never commit it to git (.gitignore excludes *.sqlite).
"""
import argparse, base64, getpass, gzip, hashlib, hmac, http.server, json, os, re, secrets, socketserver, sqlite3, ssl, sys, threading, time, urllib.parse
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone

VERSION = '1.0'
HERE = os.path.dirname(os.path.abspath(__file__))
COLLS = {'note', 'classification', 'snapshot', 'pedigree', 'watch', 'inhouse'}
PER_VARIANT = ['note', 'classification', 'snapshot', 'pedigree']
CLASSES = ['Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'Likely benign', 'Benign']
SESSION_HOURS = 12
MAX_BODY = 25 * 1024 * 1024
STATIC = {'.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8', '.css': 'text/css; charset=utf-8',
          '.json': 'application/json', '.svg': 'image/svg+xml', '.png': 'image/png', '.ico': 'image/x-icon'}
SCHEMA = '''
CREATE TABLE IF NOT EXISTS users (name TEXT PRIMARY KEY, display TEXT NOT NULL, pw TEXT, role TEXT NOT NULL DEFAULT 'reviewer',
  active INTEGER NOT NULL DEFAULT 1, created TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS sessions (token TEXT PRIMARY KEY, user TEXT NOT NULL, created TEXT NOT NULL, expires TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS docs (id INTEGER PRIMARY KEY AUTOINCREMENT, coll TEXT NOT NULL, variant TEXT, data TEXT NOT NULL,
  version INTEGER NOT NULL DEFAULT 1, created_by TEXT NOT NULL, created_at TEXT NOT NULL, updated_by TEXT, updated_at TEXT,
  deleted INTEGER NOT NULL DEFAULT 0);
CREATE INDEX IF NOT EXISTS docs_cv ON docs(coll, variant, deleted);
CREATE TABLE IF NOT EXISTS audit (id INTEGER PRIMARY KEY AUTOINCREMENT, at TEXT NOT NULL, user TEXT NOT NULL, action TEXT NOT NULL,
  coll TEXT, doc INTEGER, variant TEXT, before TEXT, after TEXT);
CREATE INDEX IF NOT EXISTS audit_v ON audit(variant);
'''
DB_PATH = None
CFG = {'trust_header': None, 'secure_cookie': False, 'root': HERE, 'page': 'brca1-variant-explorer.html'}
FAILS, FAILS_LOCK = {}, threading.Lock()
GZ, GZ_LOCK = {}, threading.Lock()


def now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace('+00:00', 'Z')


def connect():
    con = sqlite3.connect(DB_PATH, timeout=10)
    con.row_factory = sqlite3.Row
    con.execute('PRAGMA journal_mode=WAL')
    con.execute('PRAGMA busy_timeout=10000')
    con.execute('PRAGMA foreign_keys=ON')
    return con


@contextmanager
def db():
    """One connection per request: committed when the block succeeds, rolled back otherwise, always closed."""
    con = connect()
    try:
        yield con
        con.commit()
    except Exception:
        con.rollback()
        raise
    finally:
        con.close()


def init_db():
    with db() as con:
        con.executescript(SCHEMA)


# ---------------------------------------------------------------- passwords and sessions
def hash_pw(pw, salt=None, n=310000):
    salt = salt or secrets.token_bytes(16)
    dk = hashlib.pbkdf2_hmac('sha256', pw.encode(), salt, n)
    return f'pbkdf2${n}${base64.b64encode(salt).decode()}${base64.b64encode(dk).decode()}'


def check_pw(pw, stored):
    try:
        _, n, salt, dk = stored.split('$')
        return hmac.compare_digest(hash_pw(pw, base64.b64decode(salt), int(n)).split('$')[3], dk)
    except Exception:
        return False


def tok_hash(tok):
    return hashlib.sha256(tok.encode()).hexdigest()


def audit(con, user, action, coll=None, doc=None, variant=None, before=None, after=None):
    con.execute('INSERT INTO audit (at, user, action, coll, doc, variant, before, after) VALUES (?,?,?,?,?,?,?,?)',
                (now(), user, action, coll, doc, variant, None if before is None else json.dumps(before), None if after is None else json.dumps(after)))


def row_doc(r, names):
    d = {k: r[k] for k in ('id', 'coll', 'variant', 'version', 'created_by', 'created_at', 'updated_by', 'updated_at')}
    d['data'] = json.loads(r['data'])
    d['created_by_display'] = names.get(r['created_by'], r['created_by'])
    if r['updated_by']:
        d['updated_by_display'] = names.get(r['updated_by'], r['updated_by'])
    return d


class ApiError(Exception):
    def __init__(self, code, msg, extra=None):
        super().__init__(msg)
        self.code, self.msg, self.extra = code, msg, extra or {}


# ---------------------------------------------------------------- HTTP
class Handler(http.server.BaseHTTPRequestHandler):
    server_version = f'brca1-store/{VERSION}'
    protocol_version = 'HTTP/1.1'

    def log_message(self, fmt, *args):  # one line per request, without query strings (they can hold variant names only, but keep logs lean)
        sys.stderr.write('%s %s %s\n' % (time.strftime('%H:%M:%S'), self.address_string(), fmt % args))

    def headers_common(self):
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Referrer-Policy', 'no-referrer')
        self.send_header('X-Frame-Options', 'DENY')

    def send_json(self, code, obj, cookie=None):
        body = json.dumps(obj, separators=(',', ':')).encode()
        self.send_response(code)
        self.headers_common()
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Cache-Control', 'no-store')
        if cookie is not None:
            self.send_header('Set-Cookie', cookie)
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    # ------------------------------------------------ static files: the page and its data files, from --root only
    def do_GET(self):
        path = urllib.parse.urlsplit(self.path).path
        if path.startswith('/api/'):
            return self.api('GET')
        name = CFG['page'] if path in ('/', '/index.html') else urllib.parse.unquote(path.lstrip('/'))
        ext = os.path.splitext(name)[1].lower()
        if '/' in name or '\\' in name or name.startswith('.') or ext not in STATIC:
            return self.send_error(404)
        full = os.path.join(CFG['root'], name)
        if not os.path.isfile(full):
            return self.send_error(404)
        st = os.stat(full)
        etag = f'"{int(st.st_mtime)}-{st.st_size}"'
        if self.headers.get('If-None-Match') == etag:
            self.send_response(304); self.headers_common(); self.send_header('ETag', etag); self.send_header('Content-Length', '0'); self.end_headers()
            return
        gz = 'gzip' in (self.headers.get('Accept-Encoding') or '') and st.st_size > 256 * 1024
        if gz:
            with GZ_LOCK:  # compressed once per file version, kept in memory (the 70 MB data release becomes ~10 MB)
                c = GZ.get(full)
                if not c or c[0] != etag:
                    with open(full, 'rb') as f:
                        GZ[full] = c = (etag, gzip.compress(f.read(), 6))
            body = c[1]
        else:
            with open(full, 'rb') as f:
                body = f.read()
        self.send_response(200)
        self.headers_common()
        self.send_header('Content-Type', STATIC[ext])
        self.send_header('Cache-Control', 'no-cache')
        self.send_header('ETag', etag)
        if gz:
            self.send_header('Content-Encoding', 'gzip')
        self.send_header('Vary', 'Accept-Encoding')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        self.api('POST')

    def do_PUT(self):
        self.api('PUT')

    def do_DELETE(self):
        self.api('DELETE')

    # ------------------------------------------------ API
    def body(self):
        n = int(self.headers.get('Content-Length') or 0)
        if n > MAX_BODY:
            raise ApiError(413, 'request too large')
        raw = self.rfile.read(n) if n else b''
        if not raw:
            return {}
        try:
            return json.loads(raw)
        except ValueError:
            raise ApiError(400, 'invalid JSON')

    def cookie_token(self):
        for part in (self.headers.get('Cookie') or '').split(';'):
            k, _, v = part.strip().partition('=')
            if k == 'brca1sid':
                return v
        return None

    def current_user(self, con):
        if CFG['trust_header']:
            name = self.headers.get(CFG['trust_header'])
            if name:
                u = con.execute('SELECT * FROM users WHERE name=? AND active=1', (name.strip().lower(),)).fetchone()
                return u
        tok = self.cookie_token()
        if not tok:
            return None
        s = con.execute('SELECT * FROM sessions WHERE token=?', (tok_hash(tok),)).fetchone()
        if not s or s['expires'] < now():
            return None
        return con.execute('SELECT * FROM users WHERE name=? AND active=1', (s['user'],)).fetchone()

    def api(self, method):
        try:
            u = urllib.parse.urlsplit(self.path)
            q = dict(urllib.parse.parse_qsl(u.query))
            parts = [p for p in u.path.split('/') if p][1:]  # after 'api'
            if method != 'GET':
                # cross-site request protection: same-origin JavaScript sets this header, a form on another site cannot
                if self.headers.get('X-Requested-With') != 'brca1':
                    raise ApiError(403, 'missing request header')
                origin = self.headers.get('Origin')
                if origin and urllib.parse.urlsplit(origin).netloc != self.headers.get('Host'):
                    raise ApiError(403, 'cross-origin request refused')
            with db() as con:
                if parts == ['login'] and method == 'POST':
                    res = self.login(con, self.body())
                else:
                    user = self.current_user(con)
                    if not user:
                        raise ApiError(401, 'not signed in')
                    if parts == ['me']:
                        res = (200, {'name': user['name'], 'display': user['display'], 'role': user['role'], 'server': VERSION})
                    else:
                        res = self.route(con, user, method, parts, q)
            self.send_json(*res)  # after the commit
        except ApiError as e:
            self.send_json(e.code, {'error': e.msg, **e.extra})
        except Exception as e:  # never leak a traceback to the browser
            sys.stderr.write(f'error: {type(e).__name__}: {e}\n')
            self.send_json(500, {'error': 'server error'})

    def login(self, con, b):
        name, pw = str(b.get('name') or '').strip().lower(), str(b.get('password') or '')
        key = (name, self.client_address[0])
        with FAILS_LOCK:
            recent = [t for t in FAILS.get(key, []) if t > time.time() - 900]
            FAILS[key] = recent
            if len(recent) >= 5:
                raise ApiError(429, 'too many failed attempts: wait 15 minutes')
        u = con.execute('SELECT * FROM users WHERE name=? AND active=1', (name,)).fetchone()
        if not u or not u['pw'] or not check_pw(pw, u['pw']):
            with FAILS_LOCK:
                FAILS.setdefault(key, []).append(time.time())
            audit(con, name or '?', 'login_failed')
            con.commit()  # keep the failed attempt in the audit trail
            raise ApiError(401, 'wrong name or password')
        tok = secrets.token_urlsafe(32)
        exp = (datetime.now(timezone.utc) + timedelta(hours=SESSION_HOURS)).replace(microsecond=0).isoformat().replace('+00:00', 'Z')
        con.execute('DELETE FROM sessions WHERE expires < ?', (now(),))
        con.execute('INSERT INTO sessions (token, user, created, expires) VALUES (?,?,?,?)', (tok_hash(tok), u['name'], now(), exp))
        audit(con, u['name'], 'login')
        cookie = f'brca1sid={tok}; HttpOnly; SameSite=Strict; Path=/; Max-Age={SESSION_HOURS * 3600}' + ('; Secure' if CFG['secure_cookie'] else '')
        return 200, {'name': u['name'], 'display': u['display'], 'role': u['role'], 'server': VERSION}, cookie

    def route(self, con, user, method, parts, q):
        me = user['name']
        names = {r['name']: r['display'] for r in con.execute('SELECT name, display FROM users')}
        if parts == ['logout'] and method == 'POST':
            tok = self.cookie_token()
            if tok:
                con.execute('DELETE FROM sessions WHERE token=?', (tok_hash(tok),))
            audit(con, me, 'logout')
            return (200, {'ok': True}, 'brca1sid=; HttpOnly; SameSite=Strict; Path=/; Max-Age=0')
        if parts == ['variant'] and method == 'GET':
            key = q.get('key', '')
            out = {c: [] for c in PER_VARIANT}
            for r in con.execute(f"SELECT * FROM docs WHERE variant=? AND deleted=0 AND coll IN ({','.join('?' * len(PER_VARIANT))}) ORDER BY id", (key, *PER_VARIANT)):
                out[r['coll']].append(row_doc(r, names))
            return (200, out)
        if parts == ['docs'] and method == 'GET':
            coll = q.get('coll')
            if coll not in COLLS:
                raise ApiError(400, 'unknown collection')
            sql, args = 'SELECT * FROM docs WHERE coll=? AND deleted=0', [coll]
            if 'variant' in q:
                sql += ' AND variant=?'; args.append(q['variant'])
            return (200, {'docs': [row_doc(r, names) for r in con.execute(sql + ' ORDER BY id', args)]})
        if parts == ['docs'] and method == 'POST':
            b = self.body()
            d = self.create(con, me, b.get('coll'), b.get('variant'), b.get('data'))
            return (201, row_doc(d, names))
        if parts == ['docs', 'bulk'] and method == 'POST':
            b = self.body()
            items = b.get('items') or []
            if not isinstance(items, list) or len(items) > 20000:
                raise ApiError(400, 'items must be a list of at most 20,000')
            ids = [self.create(con, me, b.get('coll'), it.get('variant'), it.get('data'), log=False)['id'] for it in items]
            audit(con, me, 'bulk_create', b.get('coll'), None, None, None, {'n': len(ids)})
            return (201, {'n': len(ids)})
        if parts == ['docs'] and method == 'DELETE':
            coll = q.get('coll')
            if coll != 'inhouse':
                raise ApiError(400, 'only in-house rows can be cleared together')
            if user['role'] != 'admin':
                raise ApiError(403, 'clearing all in-house rows needs an administrator')
            n = con.execute('UPDATE docs SET deleted=1, updated_by=?, updated_at=? WHERE coll=? AND deleted=0', (me, now(), coll)).rowcount
            audit(con, me, 'bulk_delete', coll, None, None, None, {'n': n})
            return (200, {'n': n})
        if len(parts) >= 2 and parts[0] == 'docs' and parts[1].isdigit():
            r = con.execute('SELECT * FROM docs WHERE id=? AND deleted=0', (int(parts[1]),)).fetchone()
            if not r:
                raise ApiError(404, 'not found')
            if len(parts) == 2 and method == 'PUT':
                return (200, row_doc(self.update(con, user, r, self.body()), names))
            if len(parts) == 2 and method == 'DELETE':
                self.delete(con, user, r)
                return (200, {'ok': True})
            if len(parts) == 3 and method == 'POST' and parts[2] in ('check', 'return', 'withdraw'):
                return (200, row_doc(self.signoff(con, user, r, parts[2], self.body()), names))
        if parts == ['audit'] and method == 'GET':
            if 'variant' in q:
                rows = con.execute('SELECT * FROM audit WHERE variant=? ORDER BY id DESC LIMIT 500', (q['variant'],))
            else:
                rows = con.execute('SELECT * FROM audit ORDER BY id DESC LIMIT ?', (min(int(q.get('limit') or 200), 2000),))
            return (200, {'audit': [{'at': a['at'], 'user': a['user'], 'display': names.get(a['user'], a['user']), 'action': a['action'], 'coll': a['coll'],
                                                   'doc': a['doc'], 'variant': a['variant'], 'after': json.loads(a['after']) if a['after'] else None} for a in rows]})
        if parts == ['export'] and method == 'GET':
            if user['role'] != 'admin':
                raise ApiError(403, 'export needs an administrator')
            audit(con, me, 'export')
            return (200, export(con))
        raise ApiError(404, 'unknown endpoint')

    # ------------------------------------------------ document rules
    def create(self, con, me, coll, variant, data, log=True):
        if coll not in COLLS:
            raise ApiError(400, 'unknown collection')
        if not isinstance(data, dict):
            raise ApiError(400, 'data must be an object')
        variant = None if variant is None else str(variant)[:200]
        if coll in PER_VARIANT and not variant:
            raise ApiError(400, 'variant required')
        if coll == 'classification':
            if data.get('cls') not in CLASSES:
                raise ApiError(400, 'unknown classification')
            if not str(data.get('rationale') or '').strip():
                raise ApiError(400, 'a rationale is required')
            data = {k: data.get(k) for k in ('cls', 'modifier', 'rationale', 'snapshot')}
            data['status'] = 'proposed'
        if coll in ('pedigree', 'watch'):  # one per variant
            ex = con.execute('SELECT * FROM docs WHERE coll=? AND variant=? AND deleted=0', (coll, variant)).fetchone()
            if ex:
                raise ApiError(409, 'already exists', {'doc': row_doc(ex, {})})
        text = json.dumps(data)
        if len(text) > 4 * 1024 * 1024:
            raise ApiError(413, 'document too large')
        t = now()
        cur = con.execute('INSERT INTO docs (coll, variant, data, created_by, created_at) VALUES (?,?,?,?,?)', (coll, variant, text, me, t))
        if log:
            audit(con, me, 'create', coll, cur.lastrowid, variant, None, data)
        return con.execute('SELECT * FROM docs WHERE id=?', (cur.lastrowid,)).fetchone()

    def update(self, con, user, r, b):
        me, coll = user['name'], r['coll']
        if coll in ('snapshot', 'classification', 'inhouse'):
            raise ApiError(403, f'{coll} records cannot be edited')
        if coll == 'note' and r['created_by'] != me and user['role'] != 'admin':
            raise ApiError(403, 'only the author can edit a note')
        if b.get('version') != r['version']:
            raise ApiError(409, 'changed by someone else', {'doc': row_doc(r, {})})
        data = b.get('data')
        if not isinstance(data, dict):
            raise ApiError(400, 'data must be an object')
        con.execute('UPDATE docs SET data=?, version=version+1, updated_by=?, updated_at=? WHERE id=?', (json.dumps(data), me, now(), r['id']))
        audit(con, me, 'update', coll, r['id'], r['variant'], json.loads(r['data']), data)
        return con.execute('SELECT * FROM docs WHERE id=?', (r['id'],)).fetchone()

    def delete(self, con, user, r):
        me, coll = user['name'], r['coll']
        if coll in ('snapshot', 'classification'):
            raise ApiError(403, f'{coll} records are kept (audit trail)')
        if coll == 'note' and r['created_by'] != me and user['role'] != 'admin':
            raise ApiError(403, 'only the author can delete a note')
        con.execute('UPDATE docs SET deleted=1, updated_by=?, updated_at=? WHERE id=?', (me, now(), r['id']))
        audit(con, me, 'delete', coll, r['id'], r['variant'], json.loads(r['data']), None)

    def signoff(self, con, user, r, action, b):
        me = user['name']
        if r['coll'] != 'classification':
            raise ApiError(400, 'only classifications are signed off')
        data = json.loads(r['data'])
        if data.get('status') != 'proposed':
            raise ApiError(409, 'this classification is no longer awaiting sign-off')
        comment = str(b.get('comment') or '').strip()
        if action == 'withdraw':
            if r['created_by'] != me:
                raise ApiError(403, 'only the proposer can withdraw')
            data.update(status='withdrawn', closed_by=me, closed_by_display=user['display'], closed_at=now(), comment=comment)
        else:
            if r['created_by'] == me:
                raise ApiError(403, 'a second reviewer must sign off')
            if action == 'return' and not comment:
                raise ApiError(400, 'say why it is returned')
            data.update(status='checked' if action == 'check' else 'returned', closed_by=me, closed_by_display=user['display'], closed_at=now(), comment=comment)
        con.execute('UPDATE docs SET data=?, version=version+1, updated_by=?, updated_at=? WHERE id=?', (json.dumps(data), me, now(), r['id']))
        audit(con, me, action, 'classification', r['id'], r['variant'], json.loads(r['data']), data)
        return con.execute('SELECT * FROM docs WHERE id=?', (r['id'],)).fetchone()


def export(con):
    names = {r['name']: r['display'] for r in con.execute('SELECT name, display FROM users')}
    return {'exported': now(), 'server': VERSION,
            'users': [{k: r[k] for k in ('name', 'display', 'role', 'active', 'created')} for r in con.execute('SELECT * FROM users')],
            'docs': [{**row_doc(r, names), 'deleted': r['deleted']} for r in con.execute('SELECT * FROM docs ORDER BY id')],
            'audit': [dict(a) for a in con.execute('SELECT * FROM audit ORDER BY id')]}


class Server(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True
    allow_reuse_address = True


# ---------------------------------------------------------------- command line
def ask_pw(args):
    if getattr(args, 'password_stdin', False):
        return sys.stdin.readline().rstrip('\n')
    a, b = getpass.getpass('Password: '), getpass.getpass('Again: ')
    if a != b:
        sys.exit('The passwords differ.')
    return a


def main():
    global DB_PATH
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--db', default=os.path.join(HERE, 'brca1-store.sqlite'), help='SQLite file (default: brca1-store.sqlite next to this script)')
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('serve', help='serve the page, the data release and the shared store')
    s.add_argument('--host', default='127.0.0.1', help='address to listen on (0.0.0.0 for the whole network; default 127.0.0.1)')
    s.add_argument('--port', type=int, default=8000)
    s.add_argument('--root', default=HERE, help='folder with the page and brca1-data.js (default: this folder)')
    s.add_argument('--page', default='brca1-variant-explorer.html', help='page served at /')
    s.add_argument('--cert', help='TLS certificate (PEM); with --key, serves HTTPS')
    s.add_argument('--key', help='TLS private key (PEM)')
    s.add_argument('--trust-header', help='sign users in from this request header (e.g. X-Remote-User) set by a single-sign-on proxy; '
                                         'only safe when the proxy is the only way to reach this server')
    s.add_argument('--secure-cookie', action='store_true', help='mark the session cookie Secure (when an HTTPS proxy sits in front)')
    a = sub.add_parser('adduser', help='create an account')
    a.add_argument('name'); a.add_argument('--display'); a.add_argument('--admin', action='store_true'); a.add_argument('--password-stdin', action='store_true')
    p = sub.add_parser('passwd', help='set a password'); p.add_argument('name'); p.add_argument('--password-stdin', action='store_true')
    for c in ('disable', 'enable'):
        x = sub.add_parser(c, help=f'{c} an account'); x.add_argument('name')
    sub.add_parser('users', help='list accounts')
    sub.add_parser('export', help='print every record and the audit trail as JSON')
    args = ap.parse_args()
    DB_PATH = os.path.abspath(args.db)
    init_db()
    if args.cmd == 'adduser':
        name = args.name.strip().lower()
        if not re.match(r'^[a-z0-9._-]{2,40}$', name):
            sys.exit('Names: 2-40 characters, letters, digits, dot, dash or underscore.')
        with db() as con:
            if con.execute('SELECT 1 FROM users WHERE name=?', (name,)).fetchone():
                sys.exit(f'{name} already exists.')
            con.execute('INSERT INTO users (name, display, pw, role, created) VALUES (?,?,?,?,?)',
                        (name, args.display or name, hash_pw(ask_pw(args)), 'admin' if args.admin else 'reviewer', now()))
            audit(con, 'cli', 'adduser', after={'name': name, 'role': 'admin' if args.admin else 'reviewer'})
        print(f'Created {name}.')
    elif args.cmd == 'passwd':
        with db() as con:
            n = con.execute('UPDATE users SET pw=? WHERE name=?', (hash_pw(ask_pw(args)), args.name.lower())).rowcount
            con.execute('DELETE FROM sessions WHERE user=?', (args.name.lower(),))
            audit(con, 'cli', 'passwd', after={'name': args.name.lower()})
        print('Password set.' if n else 'No such account.')
    elif args.cmd in ('disable', 'enable'):
        with db() as con:
            n = con.execute('UPDATE users SET active=? WHERE name=?', (1 if args.cmd == 'enable' else 0, args.name.lower())).rowcount
            if args.cmd == 'disable':
                con.execute('DELETE FROM sessions WHERE user=?', (args.name.lower(),))
            audit(con, 'cli', args.cmd, after={'name': args.name.lower()})
        print('Done.' if n else 'No such account.')
    elif args.cmd == 'users':
        with db() as con:
            for r in con.execute('SELECT * FROM users ORDER BY name'):
                print(f"{r['name']:20s} {r['role']:9s} {'active' if r['active'] else 'disabled':9s} {r['display']}")
    elif args.cmd == 'export':
        with db() as con:
            json.dump(export(con), sys.stdout, indent=1)
    elif args.cmd == 'serve':
        CFG.update(root=os.path.abspath(args.root), page=args.page, trust_header=args.trust_header, secure_cookie=bool(args.secure_cookie or args.cert))
        if not os.path.isfile(os.path.join(CFG['root'], CFG['page'])):
            sys.exit(f"{CFG['page']} not found in {CFG['root']}")
        with db() as con:
            if not con.execute('SELECT 1 FROM users WHERE active=1').fetchone() and not args.trust_header:
                print('No account yet: create one with  python3 brca1_server.py adduser NAME --admin', file=sys.stderr)
        httpd = Server((args.host, args.port), Handler)
        scheme = 'http'
        if args.cert and args.key:
            ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
            ctx.minimum_version = ssl.TLSVersion.TLSv1_2
            ctx.load_cert_chain(args.cert, args.key)
            httpd.socket = ctx.wrap_socket(httpd.socket, server_side=True)
            scheme = 'https'
        elif args.host not in ('127.0.0.1', 'localhost'):
            print('WARNING: no --cert/--key: passwords and records travel unencrypted on the network. Use HTTPS (or a hospital HTTPS proxy).', file=sys.stderr)
        print(f"Serving {CFG['page']} and the shared store on {scheme}://{args.host}:{args.port}/  (database {DB_PATH})", file=sys.stderr)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == '__main__':
    main()
