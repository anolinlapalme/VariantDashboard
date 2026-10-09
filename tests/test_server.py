"""Shared-store rules: python3 -m unittest discover tests  (standard library only; starts brca1_server.py on a free port)."""
import http.cookiejar, json, os, socket, subprocess, sys, tempfile, time, unittest, urllib.error, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SERVER = os.path.join(ROOT, 'brca1_server.py')


def free_port():
    with socket.socket() as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


class StoreTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.db = os.path.join(cls.tmp, 'test.sqlite')
        for name, pw, extra in (('alice', 'pw-alice', ['--admin']), ('bob', 'pw-bob', [])):
            subprocess.run([sys.executable, SERVER, '--db', cls.db, 'adduser', name, '--password-stdin', *extra], input=pw + '\n', text=True, check=True, capture_output=True)
        cls.port = free_port()
        cls.proc = subprocess.Popen([sys.executable, SERVER, '--db', cls.db, 'serve', '--port', str(cls.port), '--root', ROOT], stderr=subprocess.DEVNULL)
        for _ in range(50):
            try:
                urllib.request.urlopen(f'http://127.0.0.1:{cls.port}/api/me')
            except urllib.error.HTTPError:
                break
            except OSError:
                time.sleep(0.1)

    @classmethod
    def tearDownClass(cls):
        cls.proc.terminate()
        cls.proc.wait()

    def client(self, name=None, pw=None):
        op = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()), urllib.request.ProxyHandler({}))
        if name:
            self.assertEqual(self.call(op, 'POST', '/api/login', {'name': name, 'password': pw})[0], 200)
        return op

    def call(self, op, method, path, body=None, header=True):
        req = urllib.request.Request(f'http://127.0.0.1:{self.port}{path}', method=method, data=None if body is None else json.dumps(body).encode(),
                                     headers={'Content-Type': 'application/json', **({'X-Requested-With': 'brca1'} if header else {})})
        try:
            with op.open(req) as r:
                return r.status, json.loads(r.read() or b'{}')
        except urllib.error.HTTPError as e:
            return e.code, json.loads(e.read() or b'{}')

    def test_login_required_and_csrf_header(self):
        anon = self.client()
        self.assertEqual(self.call(anon, 'GET', '/api/variant?key=c.1A>G')[0], 401)
        self.assertEqual(self.call(anon, 'POST', '/api/login', {'name': 'alice', 'password': 'pw-alice'}, header=False)[0], 403)
        self.assertEqual(self.call(anon, 'POST', '/api/login', {'name': 'alice', 'password': 'wrong'})[0], 401)

    def test_two_person_signoff(self):
        a, b = self.client('alice', 'pw-alice'), self.client('bob', 'pw-bob')
        st, d = self.call(a, 'POST', '/api/docs', {'coll': 'classification', 'variant': 'c.2A>G', 'data': {'cls': 'Likely pathogenic', 'rationale': 'test', 'status': 'checked'}})
        self.assertEqual(st, 201)
        self.assertEqual(d['data']['status'], 'proposed')  # the client cannot set the status
        self.assertEqual(self.call(a, 'POST', f"/api/docs/{d['id']}/check", {})[0], 403)  # not by the proposer
        self.assertEqual(self.call(b, 'POST', f"/api/docs/{d['id']}/return", {})[0], 400)  # a return needs a reason
        st, d2 = self.call(b, 'POST', f"/api/docs/{d['id']}/check", {'comment': 'agree'})
        self.assertEqual((st, d2['data']['status'], d2['data']['closed_by']), (200, 'checked', 'bob'))
        self.assertEqual(self.call(b, 'POST', f"/api/docs/{d['id']}/check", {})[0], 409)  # once only
        self.assertEqual(self.call(a, 'PUT', f"/api/docs/{d['id']}", {'data': {}, 'version': d2['version']})[0], 403)  # kept as is
        self.assertEqual(self.call(a, 'DELETE', f"/api/docs/{d['id']}")[0], 403)
        self.assertEqual(self.call(a, 'POST', '/api/docs', {'coll': 'classification', 'variant': 'c.2A>G', 'data': {'cls': 'Benign', 'rationale': ' '}})[0], 400)

    def test_versions_notes_and_admin(self):
        a, b = self.client('alice', 'pw-alice'), self.client('bob', 'pw-bob')
        st, p = self.call(a, 'POST', '/api/docs', {'coll': 'pedigree', 'variant': 'c.3A>G', 'data': {'families': []}})
        self.assertEqual(self.call(b, 'POST', '/api/docs', {'coll': 'pedigree', 'variant': 'c.3A>G', 'data': {'families': []}})[0], 409)  # one per variant
        self.assertEqual(self.call(b, 'PUT', f"/api/docs/{p['id']}", {'data': {'families': [1]}, 'version': 1})[0], 200)
        self.assertEqual(self.call(a, 'PUT', f"/api/docs/{p['id']}", {'data': {'families': [2]}, 'version': 1})[0], 409)  # stale copy
        st, n = self.call(a, 'POST', '/api/docs', {'coll': 'note', 'variant': 'c.3A>G', 'data': {'text': 'x'}})
        self.assertEqual(self.call(b, 'DELETE', f"/api/docs/{n['id']}")[0], 403)
        self.call(b, 'POST', '/api/docs/bulk', {'coll': 'inhouse', 'items': [{'variant': 'c.3A>G', 'data': {'hgvsc': 'c.3A>G'}}]})
        self.assertEqual(self.call(b, 'DELETE', '/api/docs?coll=inhouse')[0], 403)  # administrators only
        self.assertEqual(self.call(a, 'DELETE', '/api/docs?coll=inhouse')[0], 200)
        self.assertEqual(self.call(b, 'GET', '/api/export')[0], 403)
        st, au = self.call(a, 'GET', '/api/audit?variant=c.3A>G')
        self.assertIn('update', [x['action'] for x in au['audit']])

    def test_static_files_only(self):
        op = self.client()
        for bad in ('/test.sqlite', '/brca1_server.py', '/../README.md', '/%2e%2e/x.html', '/.git/config'):
            with self.assertRaises(urllib.error.HTTPError) as e:
                op.open(f'http://127.0.0.1:{self.port}{bad}')
            self.assertEqual(e.exception.code, 404)
        with op.open(f'http://127.0.0.1:{self.port}/') as r:
            self.assertEqual(r.status, 200)


if __name__ == '__main__':
    unittest.main()
