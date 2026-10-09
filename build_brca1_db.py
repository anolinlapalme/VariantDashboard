#!/usr/bin/env python3
"""Build brca1-data.js, a dated data release for the BRCA1 variant explorer.

The explorer page reads this file instead of querying databases each time a variant is opened, so every
reviewer sees the same values for the same release and nothing is fetched or interpreted at viewing time.

    python3 build_brca1_db.py                      # writes brca1-data.js next to this script
    python3 build_brca1_db.py --out ~/Downloads/brca1-data.js

Standard library only (Python 3.8+). Put brca1-data.js in the same folder as brca1-variant-explorer.html.

What it collects (all public sources, with the version and date of each recorded in the file):
  Ensembl     transcript, CDS and genomic sequence; VEP annotation (SpliceAI, REVEL, BayesDel, AlphaMissense,
              CADD, MaxEntScan, mutfunc, GERP, BLOSUM62, dbSNP/PubMed links) for every variant in the release;
              Compara orthologs
  gnomAD v4.1 gene variants, per-variant filters/FAF/ancestry groups, read depth around every exon
  ClinVar     BRCA1 variants (via gnomAD), aggregate records (E-utilities) and every submission (SCV) from
              the ClinVar submission_summary file
  BRCA Exchange  ENIGMA expert-panel classifications
  MaveDB      the six BRCA1 multiplexed assay score sets with their calibrations
  ClinGen CSpec  ENIGMA VCEP v1.2 Table 4, Table 9 (+ Dace & Findlay sheet) and supplementary tables
  UniProt, AlphaFold DB  features and per-residue confidence
  SpliceVault mis-splicing events seen in 300K RNA-seq samples at the splice site a variant is predicted to remove (Ensembl FTP)
  NCBI ALFA   allele counts per population from aggregated dbGaP studies

The variant set is: every gnomAD and ClinVar BRCA1 variant, every possible SNV in the coding sequence and in
the 20 intronic bases next to each exon, every saturation-genome-editing SNV, and every variant named in the
ENIGMA tables, BRCA Exchange (ENIGMA) and the French-Canadian founder list.

Time: the first run takes roughly 45-90 minutes, mostly Ensembl VEP. VEP results are cached per Ensembl release
in --cache, so later runs (to refresh ClinVar, gnomAD, BRCA Exchange) take about 10-20 minutes. Interrupted runs
resume from the cache. Downloads ClinVar's submission_summary file (~100 MB) unless --skip-scv.

Literature (unless --skip-literature): the PubMed IDs ClinVar records for each variant (var_citations.txt), the papers
NCBI LitVar2 finds for each rsID (text mining of PubMed and PMC), and PubMed title, year, journal, first author and
publication type for every paper, including those Ensembl links to the rsID. Refreshed every --lit-max-age days.
The first literature run adds roughly 5-15 minutes (LitVar2 is paced at 3 requests a second, only for rsIDs that
LitVar2's BRCA1 index lists with papers); later runs refetch only rsIDs whose paper count changed.

Changes between releases: each build compares the new release with the previous one (a small digest of every release is kept
in --releases, by default brca1-releases/ next to the output) and lists, per variant, what changed: ClinVar classification,
review status, date evaluated and submissions, the ENIGMA classification and Table 9 code, gnomAD counts and filtering allele
frequency, SGE scores and newly indexed papers. The list goes into the release (shown on the page) and into a CSV file.
With --reported FILE (your reported variants: an hgvsc column, optional classification, date and sample/report columns; the
page's in-house CSV works as is), the reported variants that changed are written as a reclassification queue (CSV).

    python3 build_brca1_db.py --reported ~/lab/brca1_reported.csv
    python3 build_brca1_db.py --changes-only --previous old/brca1-data.js   # recompute the list without downloading anything
"""
import argparse, ast, csv, gzip, hashlib, html, io, json, math, os, re, ssl, struct, sys, threading, time, urllib.error, urllib.parse, urllib.request, zipfile, zlib
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from xml.etree import ElementTree as ET

VERSION = '1.2'
GENE, CHROM = 'BRCA1', '17'
ENST, ENSP, REFSEQ, UNIPROT = 'ENST00000357654', 'ENSP00000350283', 'NM_007294.4', 'P38398'
MAVE = ['urn:mavedb:00000097-0-2', 'urn:mavedb:00000081-a-2', 'urn:mavedb:00001222-b-2', 'urn:mavedb:00001222-a-2', 'urn:mavedb:00000003-a-2', 'urn:mavedb:00000003-b-2']
SGE_URN = MAVE[0]
ENIGMA_FILES = {'table4': 'ca5cf57b-94df-4ad6-a001-c62ceccb3845', 'table9': '0a35d6a8-5050-44b6-8a9d-babe8cdc06b2', 'supp': 'cb4a09fe-30f4-4aa8-9d76-d7ea407c9754'}
ENIGMA_SHEETS = {'table4': [r'Annotated Exons'], 'table9': [r'Table9', r'Dace'], 'supp': [r'^ST2 ', r'^ST4 ', r'^ST7 ', r'^ST10 ', r'^ST16 ']}
SPECIES = ['pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii', 'macaca_mulatta', 'callithrix_jacchus', 'mus_musculus', 'rattus_norvegicus', 'oryctolagus_cuniculus',
           'canis_lupus_familiaris', 'felis_catus', 'equus_caballus', 'bos_taurus', 'sus_scrofa', 'loxodonta_africana', 'dasypus_novemcinctus', 'monodelphis_domestica',
           'ornithorhynchus_anatinus', 'gallus_gallus', 'taeniopygia_guttata', 'chrysemys_picta_bellii', 'anolis_carolinensis', 'xenopus_tropicalis', 'latimeria_chalumnae',
           'lepisosteus_oculatus', 'danio_rerio', 'oryzias_latipes', 'takifugu_rubripes', 'callorhinchus_milii', 'eptatretus_burgeri']
# BRCA1 rows of the French-Canadian founder list (Behl et al. 2020), HGVS as published
FC_BRCA1 = ['c.962G>A', 'c.1016dupA', 'c.1961dupA', 'c.2125_2126insA', 'c.2834_2836delGTAinsC', 'c.3649_3650insA', 'c.3756_3759delGTCT', 'c.4327C>T', 'c.5102_5103delTG', 'c.1054G>T', 'c.5536C>T']
VEP_PARAMS = {'transcript_id': ENST, 'numbers': 1, 'hgvs': 1, 'canonical': 1, 'mane': 1, 'protein': 1, 'AlphaMissense': 1, 'CADD': 1, 'REVEL': 1, 'SpliceAI': 2, 'MaxEntScan': 1,
              'Conservation': 1, 'mutfunc': 1, 'Blosum62': 1, 'dbNSFP': 'BayesDel_noAF_score', 'pubmed': 1, 'variant_class': 1}
TC_KEEP = ['transcript_id', 'hgvsc', 'hgvsp', 'consequence_terms', 'exon', 'intron', 'protein_start', 'protein_end', 'amino_acids', 'spliceai', 'revel', 'cadd_phred',
           'alphamissense', 'bayesdel_noaf_score', 'maxentscan_ref', 'maxentscan_alt', 'maxentscan_diff', 'mutfunc', 'conservation', 'blosum62']
BX_KEEP = ['id', 'Gene_Symbol', 'HGVS_cDNA', 'Clinical_significance_ENIGMA', 'Date_last_evaluated_ENIGMA', 'Comment_on_clinical_significance_ENIGMA',
           'Assertion_method_citation_ENIGMA', 'Pathogenicity_all', 'Pathogenicity_expert']
INTRON_FLANK = 20
LITVAR = 'https://www.ncbi.nlm.nih.gov/research/litvar2-api'
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
# PubMed publication types worth showing next to a title
PUBTYPES = {'Review': 'review', 'Systematic Review': 'review', 'Meta-Analysis': 'meta-analysis', 'Case Reports': 'case report', 'Clinical Trial': 'trial',
            'Practice Guideline': 'guideline', 'Guideline': 'guideline', 'Letter': 'letter', 'Comment': 'comment', 'Editorial': 'editorial', 'Preprint': 'preprint',
            'Retracted Publication': 'retracted', 'Published Erratum': 'erratum'}
UA = f'brca1-explorer-build/{VERSION} (python-urllib)'
try:  # certifi, when installed, avoids certificate errors with python.org builds on macOS
    import certifi
    SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    SSL_CTX = ssl.create_default_context()
CERT_HELP = ('Python could not verify HTTPS certificates. On macOS with the python.org installer, run '
             '"Install Certificates.command" in /Applications/Python 3.x/, or run: pip3 install certifi. Then rerun this script.')

def now_iso():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace('+00:00', 'Z')

def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)

# ---------------------------------------------------------------- HTTP with retries and per-host pacing
_last, _pace_lock = {}, threading.Lock()
PACE = {'gnomad.broadinstitute.org': 6.5, 'eutils.ncbi.nlm.nih.gov': 0.4, 'rest.ensembl.org': 0.1, 'brcaexchange.org': 0.5, 'www.ncbi.nlm.nih.gov': 0.34}

def http(url, data=None, headers=None, timeout=180, tries=6, binary=False, meta=None):
    host = urllib.parse.urlparse(url).hostname or ''
    for attempt in range(tries):
        with _pace_lock:  # reserve the next slot for this host
            slot = max(time.time(), _last.get(host, 0) + PACE.get(host, 0))
            _last[host] = slot
        if slot > time.time():
            time.sleep(slot - time.time())
        try:
            req = urllib.request.Request(url, data=data, headers={'User-Agent': UA, **(headers or {})})
            with urllib.request.urlopen(req, timeout=timeout, context=SSL_CTX) as r:
                body = r.read()
                if meta is not None:
                    meta.update({k.lower(): v for k, v in r.headers.items()})
                return body if binary else body.decode('utf-8')
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < tries - 1:
                wait = float(e.headers.get('Retry-After') or 0) or min(120, 5 * 2 ** attempt)
                log(f'  HTTP {e.code} from {host}, retrying in {wait:.0f} s')
                time.sleep(wait)
                continue
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            if isinstance(getattr(e, 'reason', e), ssl.SSLCertVerificationError):
                sys.exit(CERT_HELP)
            if attempt < tries - 1:
                log(f'  {type(e).__name__} from {host} ({e}), retrying')
                time.sleep(min(120, 5 * 2 ** attempt))
                continue
            raise

def get_json(url, **kw):
    return json.loads(http(url, headers={'Accept': 'application/json'}, **kw))

def post_json(url, obj, **kw):
    return json.loads(http(url, data=json.dumps(obj).encode(), headers={'Content-Type': 'application/json', 'Accept': 'application/json'}, **kw))

def gnomad(query):
    j = post_json('https://gnomad.broadinstitute.org/api', {'query': query}, timeout=240)
    if j.get('errors') and not j.get('data'):
        raise RuntimeError('gnomAD: ' + json.dumps(j['errors'])[:300])
    return j

# ---------------------------------------------------------------- cache
# Volatile sources (ClinVar, BRCA Exchange, gnomAD's ClinVar list, MaveDB, ENIGMA files, UniProt) are refetched when
# their cached copy is older than --max-age days; release-bound data (Ensembl/VEP per release, gnomAD v4.1 variant
# details and coverage) are kept. --refresh NAME... forces groups again ("all" for everything).
class Cache:
    def __init__(self, root, refresh, max_age_days):
        self.root, self.refresh, self.max_age = root, set(refresh or []), max_age_days * 86400
        os.makedirs(root, exist_ok=True)

    def path(self, name, ext='.json'):
        return os.path.join(self.root, re.sub(r'[^A-Za-z0-9._-]+', '_', name) + ext)

    def fresh(self, p, group, volatile):
        if not os.path.exists(p) or self.refresh & {'all', group}:
            return False
        return not volatile or time.time() - os.path.getmtime(p) < self.max_age

    def get(self, name, fn, group=None, volatile=False):
        p = self.path(name)
        if self.fresh(p, group or name, volatile):
            with open(p) as f:
                return json.load(f)
        val = fn()
        tmp = p + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(val, f)
        os.replace(tmp, p)
        return val

    def store(self, name, group):
        """Append-only per-key store (JSON lines), so interrupted runs resume and new inputs only cost new calls."""
        p = self.path(name, '.jsonl')
        if self.refresh & {'all', group} and os.path.exists(p):
            os.replace(p, p + '.old')
        return Store(p)


class Store(dict):
    def __init__(self, path):
        super().__init__()
        self.path = path
        if os.path.exists(path):
            with open(path, encoding='utf-8') as f:
                for line in f:
                    try:
                        self.update(json.loads(line))
                    except ValueError:
                        pass  # a line cut short by an interruption

    _lock = threading.Lock()

    def add(self, items):
        with self._lock, open(self.path, 'a', encoding='utf-8') as f:
            self.update(items)
            for k, v in items.items():
                f.write(json.dumps({k: v}, separators=(',', ':')) + '\n')

# ---------------------------------------------------------------- transcript model (mirrors the page's buildTx)
def build_tx(lk):
    minus = lk['strand'] == -1
    ex = sorted(({'start': e['start'], 'end': e['end'], 'len': e['end'] - e['start'] + 1} for e in lk['Exon']), key=lambda e: -e['start'] if minus else e['start'])
    cd = 0
    for i, e in enumerate(ex):
        e['n'] = i + 1; e['cdnaStart'] = cd + 1; cd += e['len']; e['cdnaEnd'] = cd

    def g_to_cdna(g):
        for e in ex:
            if e['start'] <= g <= e['end']:
                return e['cdnaStart'] + (e['end'] - g if minus else g - e['start'])
        return None
    c1g = lk['Translation']['end'] if minus else lk['Translation']['start']
    lastg = lk['Translation']['start'] if minus else lk['Translation']['end']
    c1 = g_to_cdna(c1g)
    cds_len = g_to_cdna(lastg) - c1 + 1

    def cdna_to_g(cdna):
        for e in ex:
            if e['cdnaStart'] <= cdna <= e['cdnaEnd']:
                return e['end'] - (cdna - e['cdnaStart']) if minus else e['start'] + (cdna - e['cdnaStart'])
        return None

    def c_to_g(c):
        m = re.match(r'c\.([-*]?)(\d+)(?:([+-])(\d+))?', c)
        if not m:
            return None
        n = int(m.group(2))
        cdna = c1 - n if m.group(1) == '-' else c1 + cds_len - 1 + n if m.group(1) == '*' else c1 + n - 1
        g = cdna_to_g(cdna)
        if g is None:
            return None
        if m.group(3):
            off = int(m.group(4)) * (1 if m.group(3) == '+' else -1)
            g += -off if minus else off
        return g
    return {'minus': minus, 'ex': ex, 'c1': c1, 'cds_len': cds_len, 'cdna_to_g': cdna_to_g, 'c_to_g': c_to_g}

COMP = str.maketrans('ACGT', 'TGCA')

# ---------------------------------------------------------------- tabix over HTTP: only the gene region of a genome-wide file is downloaded
SV_URL = 'https://ftp.ensembl.org/pub/current_variation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz'
ALFA_URL = 'https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/freq.vcf.gz'
# ALFA population BioSample IDs (NCBI ALFA documentation)
ALFA_POPS = {'SAMN10492705': 'Total', 'SAMN10492695': 'European', 'SAMN10492703': 'African', 'SAMN10492698': 'African American', 'SAMN10492696': 'African Others',
             'SAMN10492704': 'Asian', 'SAMN10492697': 'East Asian', 'SAMN10492701': 'Other Asian', 'SAMN10492702': 'South Asian', 'SAMN10492699': 'Latin American 1',
             'SAMN10492700': 'Latin American 2', 'SAMN11605645': 'Other'}


def range_get(url, a, b):
    """Bytes a..b (inclusive) of a remote file. Refuses a server that ignores Range, so a genome-wide file is never downloaded whole."""
    host = urllib.parse.urlparse(url).hostname or ''
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA, 'Range': f'bytes={a}-{b}'})
            with urllib.request.urlopen(req, timeout=300, context=SSL_CTX) as r:
                if r.status != 206:
                    raise RuntimeError(f'{host} ignored the byte-range request (HTTP {r.status})')
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            if isinstance(getattr(e, 'reason', e), ssl.SSLCertVerificationError):
                sys.exit(CERT_HELP)
            if attempt == 4:
                raise
            log(f'  {type(e).__name__} from {host}, retrying')
            time.sleep(5 * 2 ** attempt)


def bgzf_blocks(raw, base=0, stop=None):
    """Decompresses consecutive BGZF blocks; returns the text and {compressed offset: offset in the text}."""
    data, starts, q = bytearray(), {}, 0
    while q + 18 <= len(raw) and (stop is None or base + q <= stop):
        if raw[q:q + 4] != b'\x1f\x8b\x08\x04':
            break
        bsize = struct.unpack_from('<H', raw, q + 16)[0] + 1
        if q + bsize > len(raw):
            break
        starts[base + q] = len(data)
        data += zlib.decompress(raw[q + 18:q + bsize - 8], -15)
        q += bsize
    return data, starts


def tabix_fetch(url, seqnames, beg, end, fetch=None, tbi=None):
    """Tab-split rows of a bgzipped, tabix-indexed file overlapping beg..end (1-based) on the first of seqnames present in the index."""
    fetch = fetch or range_get
    t = gzip.decompress(tbi if tbi is not None else http(url + '.tbi', binary=True, timeout=600))
    magic, n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from('<4s8i', t, 0)
    if magic != b'TBI\x01':
        raise RuntimeError('not a tabix index')
    names = [n.decode() for n in t[36:36 + l_nm].split(b'\0')[:-1]]
    p, refs = 36 + l_nm, []
    for _ in range(n_ref):
        n_bin, = struct.unpack_from('<i', t, p); p += 4
        bins = {}
        for _ in range(n_bin):
            b, n_chunk = struct.unpack_from('<Ii', t, p); p += 8
            bins[b] = [struct.unpack_from('<QQ', t, p + 16 * k) for k in range(n_chunk)]; p += 16 * n_chunk
        n_intv, = struct.unpack_from('<i', t, p); p += 4
        refs.append((bins, struct.unpack_from(f'<{n_intv}Q', t, p))); p += 8 * n_intv
    rid = next((names.index(s) for s in seqnames if s in names), None)
    if rid is None:
        raise RuntimeError(f'none of {seqnames} in the index (it has {", ".join(names[:4])}…)')
    bins, ioff = refs[rid]
    b0, e0 = beg - 1, end  # 0-based, half-open
    want = [0]
    for shift, off in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        want.extend(range(off + (b0 >> shift), off + ((e0 - 1) >> shift) + 1))
    min_off = ioff[b0 >> 14] if (b0 >> 14) < len(ioff) else 0
    merged = []
    for c in sorted(c for w in want for c in bins.get(w, []) if c[1] > min_off):
        if merged and c[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], c[1]))
        else:
            merged.append(c)
    rows, seen = [], set()
    for cb, ce in merged:
        cs, cend = cb >> 16, ce >> 16
        data, starts = bgzf_blocks(fetch(url, cs, cend + 65535), cs, cend)  # a BGZF block is at most 64 KiB
        a = starts.get(cs, 0) + (cb & 0xFFFF)
        z = starts[cend] + (ce & 0xFFFF) if cend in starts else len(data)
        for ln in bytes(data[a:z]).decode('utf-8', 'replace').split('\n'):
            if not ln or ln.startswith('#') or ln in seen:
                continue
            f = ln.split('\t')
            try:
                pos = int(f[col_beg - 1])
            except (IndexError, ValueError):
                continue
            if f[col_seq - 1] == names[rid] and beg <= pos <= end:
                seen.add(ln); rows.append(f)
    return rows


def vcf_header(url, fetch=None):
    fetch = fetch or range_get
    for size in (1 << 18, 1 << 21, 1 << 23):
        text = bytes(bgzf_blocks(fetch(url, 0, size - 1))[0]).decode('utf-8', 'replace')
        lines = text.split('\n')
        for i, ln in enumerate(lines):
            if ln.startswith('#CHROM'):
                return lines[:i + 1]
    raise RuntimeError('VCF header not found')


def _num(x):
    try:
        v = float(x)
        return int(v) if v.is_integer() and abs(v) > 1 else round(v, 4)
    except (TypeError, ValueError):
        return None if x in (None, '', 'NA', '.') else x


def splicevault_region(g0, g1, fetch=None, tbi=None):
    """SpliceVault (Dawes et al. 2023) rows for the MANE transcript: the mis-splicing events seen at the splice site a variant is predicted to remove."""
    rows = tabix_fetch(SV_URL, [CHROM, 'chr' + CHROM], g0, g1, fetch, tbi)
    out = {}
    for f in rows:
        if len(f) < 12 or f[4].split('.')[0] != ENST:
            continue
        # chromosome, position, ref, alt, transcript, site type, site position, SpliceAI delta, out-of-frame fraction, top events, GTEx samples, max depth
        out[f'{CHROM}-{f[1]}-{f[2]}-{f[3]}'] = [f[5], _num(f[6]), _num(f[7]), f[8], f[9], _num(f[10]), _num(f[11])]
    return {'v': out}


def alfa_region(g0, g1, fetch=None, tbi=None):
    """NCBI ALFA allele counts per population (aggregated dbGaP studies) for the gene region."""
    head = vcf_header(ALFA_URL, fetch)
    samples = head[-1].lstrip('#').split('\t')[9:]
    rel = next((ln.split('=', 1)[1] for ln in head if re.match(r'##(fileDate|release|source)=', ln, re.I)), '')
    order = [s for s in ALFA_POPS if s in samples] + [s for s in samples if s not in ALFA_POPS]
    idx = {s: samples.index(s) for s in order}
    out = {}
    for f in tabix_fetch(ALFA_URL, ['NC_0000%s.11' % CHROM, CHROM, 'chr' + CHROM], g0, g1, fetch, tbi):
        if len(f) < 10:
            continue
        fmt = f[8].split(':')
        if 'AN' not in fmt or 'AC' not in fmt:
            continue
        ian, iac = fmt.index('AN'), fmt.index('AC')
        cells = [f[9 + idx[s]].split(':') for s in order]
        for k, alt in enumerate(f[4].split(',')):
            counts = []
            for c in cells:
                an = int(c[ian]) if len(c) > ian and c[ian].isdigit() else 0
                acs = c[iac].split(',') if len(c) > iac else []
                counts.append([int(acs[k]) if k < len(acs) and acs[k].isdigit() else 0, an])
            if any(ac for ac, _ in counts):
                out[f'{CHROM}-{f[1]}-{f[3]}-{alt}'] = [f[2], counts]
    return {'release': rel, 'pops': [ALFA_POPS.get(s, s) for s in order], 'v': out}


# ---------------------------------------------------------------- xlsx reader (values as SheetJS sheet_to_json(header: 1, blankrows: false, defval: null))
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main', 'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}

def col_index(ref):
    n = 0
    for ch in re.match(r'[A-Z]+', ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1

def read_xlsx(data, sheet_patterns):
    z = zipfile.ZipFile(io.BytesIO(data))
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    target = {r.get('Id'): r.get('Target') for r in rels}
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS):
            shared.append(''.join(t.text or '' for t in si.iter('{%s}t' % NS['m'])))
    out = {}
    for sh in wb.find('m:sheets', NS):
        name = sh.get('name')
        if not any(re.search(p, name) for p in sheet_patterns):
            continue
        t = target[sh.get('{%s}id' % NS['r'])].lstrip('/')
        path = t if t.startswith('xl/') else 'xl/' + t
        root = ET.fromstring(z.read(path))
        rows = []
        for row in root.iter('{%s}row' % NS['m']):
            vals = {}
            for c in row.findall('m:c', NS):
                typ, v = c.get('t'), c.find('m:v', NS)
                if typ == 'inlineStr':
                    val = ''.join(x.text or '' for x in c.iter('{%s}t' % NS['m']))
                elif v is None or v.text is None:
                    continue
                elif typ == 's':
                    val = shared[int(v.text)]
                elif typ in ('str', 'e', 'd'):
                    val = v.text
                elif typ == 'b':
                    val = v.text == '1'
                else:
                    f = float(v.text)
                    val = int(f) if f.is_integer() and abs(f) < 2 ** 53 else f
                if val == '' or val is None:
                    continue
                vals[col_index(c.get('r'))] = val
            if vals:
                width = max(vals) + 1
                rows.append([vals.get(i) for i in range(width)])
        out[name] = rows
    return out

# ---------------------------------------------------------------- helpers
def norm_c(c):
    c = re.sub(r'^.*:(?=c\.)', '', str(c).strip())
    return re.sub(r'(del|dup)[ACGT]+$', r'\1', c, flags=re.I)

def iso_date(s):
    s = (s or '').strip()
    if not s or s == '-':
        return None
    for fmt in ('%b %d, %Y', '%Y-%m-%d', '%Y/%m/%d %H:%M', '%Y/%m/%d'):
        try:
            return datetime.strptime(s, fmt).strftime('%Y-%m-%d')
        except ValueError:
            pass
    return s

def slim_vep(r):
    tc = next((t for t in r.get('transcript_consequences') or [] if t.get('transcript_id') == ENST), None)
    if not tc:
        return None
    t = {k: tc[k] for k in TC_KEEP if k in tc}
    if isinstance(t.get('spliceai'), dict):
        t['spliceai'] = {k: v for k, v in t['spliceai'].items() if k[:3] in ('DS_', 'DP_')}
    if isinstance(t.get('mutfunc'), dict):
        t['mutfunc'] = {k: {kk: vv for kk, vv in v.items() if kk in ('ddG', 'evidence')} for k, v in t['mutfunc'].items() if isinstance(v, dict) and k in ('exp', 'mod', 'int')}
    if 'hgvsp' in t:
        t['hgvsp'] = urllib.parse.unquote(t['hgvsp'])
    coloc = [{'id': c['id'], **({'pubmed': c['pubmed']} if c.get('pubmed') else {})} for c in r.get('colocated_variants') or [] if str(c.get('id', '')).startswith('rs')]
    out = {'input': r.get('input'), 'allele_string': r.get('allele_string'), 'start': r.get('start'), 'end': r.get('end'), 'strand': r.get('strand'), 'transcript_consequences': [t]}
    if coloc:
        out['colocated_variants'] = coloc
    return out

def gid_of(v):
    ra, aa = (v.get('allele_string') or '/').split('/')[:2]
    if len(ra) == 1 and len(aa) == 1 and '-' not in (ra, aa):
        if v.get('strand') == -1:
            ra, aa = ra.translate(COMP), aa.translate(COMP)
        return f"{CHROM}-{v['start']}-{ra}-{aa}"
    return None

# ---------------------------------------------------------------- steps
def vep_annotate(items, kind, store):
    """items: VCF-style lines (kind='region') or HGVS (kind='hgvs'). Results are kept per input in `store` (slimmed)."""
    url = f'https://rest.ensembl.org/vep/human/{kind}'
    todo = [x for x in items if x not in store]
    size, done = 200, 0

    def run(batch):
        key = 'variants' if kind == 'region' else 'hgvs_notations'
        try:
            res = post_json(url, {key: batch, **VEP_PARAMS}, timeout=600)
        except urllib.error.HTTPError as e:
            if e.code == 400 and len(batch) > 1:  # one unparsable notation fails the whole batch: split until it is isolated
                half = len(batch) // 2
                run(batch[:half]); run(batch[half:])
                return
            if e.code == 400:
                log(f'  VEP could not annotate {batch[0]}')
                store.add({batch[0]: None})
                return
            raise
        got = {r.get('input'): slim_vep(r) for r in res}
        store.add({x: got.get(x) for x in batch})
    if todo:
        log(f'  VEP {kind}: {len(todo)} to annotate ({len(items) - len(todo)} cached), 3 requests at a time')
    batches = [todo[i:i + size] for i in range(0, len(todo), size)]
    with ThreadPoolExecutor(max_workers=3) as ex:
        for k, _ in enumerate(ex.map(run, batches)):
            if (k + 1) % 5 == 0 or k + 1 == len(batches):
                log(f'  VEP {kind}: {min((k + 1) * size, len(todo))}/{len(todo)}')
    return {x: store.get(x) for x in items}

def parse_lines(text):
    """LitVar2 returns one record per line (JSON, or Python-style dicts); a JSON array is accepted too."""
    t = text.strip()
    if t.startswith('['):
        return json.loads(t)
    out = []
    for line in t.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except ValueError:
            try:
                out.append(ast.literal_eval(line))
            except (ValueError, SyntaxError):
                pass
    return out


def literature(a, V, vid_of, key):
    """Adds rec['lit'] = {'cv': ClinVar-cited PMIDs, 'lv': LitVar2 PMIDs, 'rsShared': other alleles of the rsID} and returns PubMed summaries."""
    cache = Cache(a.cache, a.refresh, a.lit_max_age)
    # ClinVar: PubMed IDs per VariationID (allele-specific)
    log('Literature: ClinVar citations (var_citations.txt, streamed and filtered)')

    def fetch_cit():
        want, out, meta = set(vid_of.values()), {}, {}
        req = urllib.request.Request('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt', headers={'User-Agent': UA})
        with urllib.request.urlopen(req, timeout=600, context=SSL_CTX) as r:
            meta['modified'] = r.headers.get('Last-Modified')
            iv = isrc = iid = None
            for raw in io.TextIOWrapper(r, encoding='utf-8', errors='replace'):
                if raw.startswith('#'):
                    hdr = raw.lstrip('#').rstrip('\n').split('\t')
                    iv, isrc, iid = (hdr.index(k) if k in hdr else None for k in ('VariationID', 'citation_source', 'citation_id'))
                    continue
                if iv is None:
                    continue
                f = raw.rstrip('\n').split('\t')
                if len(f) > max(iv, isrc, iid) and f[iv] in want and f[isrc] == 'PubMed' and f[iid].isdigit():
                    out.setdefault(f[iv], set()).add(int(f[iid]))
        if iv is None:
            raise RuntimeError('var_citations.txt: header line not found')
        return {'cit': {k: sorted(v) for k, v in out.items()}, 'modified': meta.get('modified')}
    try:
        res = cache.get('clinvar-citations', fetch_cit, group='literature', volatile=True)
        cit, cit_mod = res['cit'], res.get('modified')
        log(f'  {sum(len(x) for x in cit.values())} citations for {len(cit)} variants')
    except Exception as e:
        log(f'  WARNING: ClinVar citations unavailable ({e})')
        cit, cit_mod = {}, None

    # LitVar2: the gene index lists every BRCA1 rsID with its paper count; papers are fetched only for rsIDs with papers
    rs_of = {}
    for c, rec in V.items():
        for x in rec['vep'].get('colocated_variants') or []:
            rs_of.setdefault(x['id'].lower(), []).append(c)
    log('Literature: LitVar2 gene index')
    try:
        idx = cache.get(f'litvar-gene-{GENE}', lambda: parse_lines(http(f'{LITVAR}/variant/search/gene/{GENE}', timeout=300)), group='literature', volatile=True)
        counts = {str(x.get('rsid') or '').lower(): int(x.get('pmids_count') or 0) for x in idx if isinstance(x, dict) and x.get('rsid')}
        todo = sorted(r for r in rs_of if counts.get(r, 0) > 0)
        log(f'  {len(counts)} BRCA1 rsIDs in LitVar2; {len(todo)} of this release have papers')
    except Exception as e:  # fall back to rsIDs that ClinVar or Ensembl already link to papers
        log(f'  LitVar2 gene index unavailable ({e}); querying rsIDs already linked to papers')
        counts = {}
        todo = sorted(r for r, cs in rs_of.items() if any(cit.get(vid_of.get(c)) or any(x.get('pubmed') for x in V[c]['vep'].get('colocated_variants') or [] if x['id'].lower() == r) for c in cs))
    lvstore = cache.store('litvar-publications', 'literature')
    kf = lambda r: f'{r}#{counts.get(r, "?")}'  # a changed paper count refetches that rsID
    need = [r for r in todo if kf(r) not in lvstore]
    if need:
        log(f'  LitVar2: {len(need)} rsIDs to fetch ({len(todo) - len(need)} cached; about {len(need) * PACE.get("www.ncbi.nlm.nih.gov", 0.34) / 60:.0f} min)')
    lv_errors = []

    def fetch_pub(r):
        try:
            j = get_json(f'{LITVAR}/variant/get/{urllib.parse.quote("litvar@" + r + "##", safe="")}/publications', timeout=120)
            return r, sorted({int(p) for p in j.get('pmids') or [] if str(p).isdigit()})
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return r, []
            lv_errors.append(e); return r, None
        except Exception as e:
            lv_errors.append(e); return r, None
    with ThreadPoolExecutor(max_workers=3) as ex:
        for n, (r, pm) in enumerate(ex.map(fetch_pub, need), 1):
            if pm is not None:
                lvstore.add({kf(r): pm})
            if n % 250 == 0:
                log(f'  {n}/{len(need)}')
    if lv_errors:
        log(f'  WARNING: {len(lv_errors)} LitVar2 requests failed ({lv_errors[0]}); rerun to complete them')

    # attach to records
    allp = set()
    for c, rec in V.items():
        rs = [x['id'].lower() for x in rec['vep'].get('colocated_variants') or []]
        lv = sorted({p for r in rs for p in lvstore.get(kf(r)) or []})
        cv_ = cit.get(vid_of.get(c)) or []
        lit = {k: v for k, v in (('cv', cv_), ('lv', lv)) if v}
        shared = sorted({o for r in rs for o in rs_of.get(r, []) if o != c})
        if lit and shared:
            lit['rsShared'] = shared[:8]
        if lit:
            rec['lit'] = lit
        allp.update(cv_, lv)
        for x in rec['vep'].get('colocated_variants') or []:
            allp.update(int(p) for p in x.get('pubmed') or [] if str(p).isdigit())

    # PubMed summaries (kept; titles do not change)
    log(f'Literature: PubMed summaries for {len(allp)} papers')
    pstore = cache.store('pubmed-summaries', 'pubmed')
    need = sorted(p for p in allp if str(p) not in pstore)
    for i in range(0, len(need), 200):
        chunk = need[i:i + 200]
        body = urllib.parse.urlencode({'db': 'pubmed', 'id': ','.join(map(str, chunk)), 'retmode': 'json'}).encode()
        try:
            j = json.loads(http(f'{EUTILS}/esummary.fcgi?tool=brca1-explorer' + key, data=body, headers={'Content-Type': 'application/x-www-form-urlencoded'}))
        except Exception as e:
            log(f'  WARNING: PubMed summaries incomplete ({e}); rerun to complete them'); break
        got = {}
        for p in chunk:
            r = (j.get('result') or {}).get(str(p))
            if not r or r.get('error'):
                got[str(p)] = None
                continue
            title = html.unescape(re.sub(r'<[^>]+>', '', r.get('title') or '')).strip()
            year = (r.get('sortpubdate') or r.get('pubdate') or '')[:4]
            au = ((r.get('authors') or [{}])[0] or {}).get('name') or ''
            types = sorted({PUBTYPES[t] for t in r.get('pubtype') or [] if t in PUBTYPES})
            got[str(p)] = [title[:400], int(year) if year.isdigit() else None, r.get('source') or '', au, ';'.join(types)]
        pstore.add(got)
        if (i // 200 + 1) % 25 == 0:
            log(f'  {min(i + 200, len(need))}/{len(need)}')
    pubs = {str(p): pstore[str(p)] for p in sorted(allp) if pstore.get(str(p))}
    n_lv = sum(1 for r in todo if lvstore.get(kf(r)))
    note = (f'ClinVar var_citations.txt{f" ({cit_mod})" if cit_mod else ""}; LitVar2 ({n_lv} BRCA1 rsIDs with papers); '
            f'Ensembl rsID citations; PubMed summaries for {len(pubs)} papers')
    return pubs, note


# ---------------------------------------------------------------- changes between releases (and the reclassification queue)
# Every build keeps a small digest of the values that matter for reclassification in --releases and compares the new release
# with the previous one. The result goes into the release (db['changes'], listed on the page) and into CSV files; with
# --reported, the variants your laboratory has reported that changed are written as a reclassification queue.
FAF_LEVELS = (2e-5, 1e-4, 1e-3)  # frequency cut-offs in the ENIGMA / CanVIG-UK specifications: crossing one is flagged
SGE_BINS = (-1.328, -0.748)      # Findlay et al. 2018: non-functional <= -1.328, functional >= -0.748
SRC_NAMES = {'vep': 'Ensembl VEP', 'compara': 'Ensembl Compara', 'uniprot': 'UniProt', 'alphafold': 'AlphaFold DB', 'gnomad': 'gnomAD', 'clinvar': 'ClinVar',
             'bx': 'BRCA Exchange (ENIGMA)', 'mavedb': 'MaveDB', 'enigma': 'ENIGMA VCEP tables', 'literature': 'Literature indexes',
             'splicevault': 'SpliceVault', 'alfa': 'NCBI ALFA'}
SUPS = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')


def sci(x):
    if x is None:
        return '—'
    if x == 0:
        return '0'
    e = math.floor(math.log10(abs(x)))
    return f'{x:.2g}' if e >= -2 else f'{x / 10 ** e:.1f}×10{str(e).translate(SUPS)}'


def natkey(c):
    return [int(x) if x.isdigit() else x for x in re.split(r'(\d+)', c)]


def tier(cls):
    s = (cls or '').lower()
    if 'conflicting' in s:
        return 'conflicting'
    if 'pathogenic' in s and 'benign' not in s:
        return 'P/LP'
    if 'uncertain' in s or 'vus' in s:
        return 'VUS'
    if 'benign' in s:
        return 'B/LB'
    return None


def faf_of(gn):
    """gnomAD v4.1 joint (exomes + genomes) grpmax FAF95; only without a joint value, the higher of exomes and genomes. The page uses the same rule."""
    if not gn:
        return None, None
    j = (gn.get('joint') or {}).get('faf95') or {}
    if j.get('popmax') is not None:
        return j['popmax'], j.get('popmax_population')
    best = (None, None)
    for k in ('exome', 'genome'):
        f = (gn.get(k) or {}).get('faf95') or {}
        if f.get('popmax') is not None and (best[0] is None or f['popmax'] > best[0]):
            best = (f['popmax'], f.get('popmax_population'))
    return best


def enigma_t9(enigma):
    """ENIGMA Table 9 code and weight per variant, read as the page reads them."""
    rows = next((r for n, r in (enigma.get('table9') or {}).items() if re.search('Table9', n)), [])
    hi = next((i for i, r in enumerate(rows) if r and r[0] == 'Gene' and len(r) > 1 and 'HGVS Nucleotide' in str(r[1] or '')), -1)
    out = {}
    for r in rows[hi + 1:]:
        if r and r[0] == GENE and len(r) > 1 and r[1]:
            out[re.sub(r'^[^:\s]*:(?=c\.)', '', str(r[1]).strip())] = [r[3] if len(r) > 3 else None, r[4] if len(r) > 4 else None]
    return out


def release_digest(db):
    """The values compared between releases, per variant (a few MB, kept gzipped in --releases)."""
    meta = db['meta']
    gene = ((db.get('gene') or {}).get('gnomadGene') or {}).get('data', {}).get('gene') or {}
    stars = {str(v['clinvar_variation_id']): v.get('gold_stars') for v in gene.get('clinvar_variants') or [] if v.get('clinvar_variation_id')}
    t9 = enigma_t9(db.get('enigma') or {})
    V = {}
    for c, r in db['v'].items():
        d = {}
        if r.get('cv'):
            gc = r['cv'].get('gc') or {}
            acc = (r.get('vcv') or {}).get('acc') or ''
            d['cv'] = [gc.get('description'), gc.get('review_status'), iso_date(gc.get('last_evaluated')), stars.get(str(int(acc[3:])) if acc[3:].isdigit() else None)]
        if r.get('scv'):
            d['scv'] = [[s.get('scv'), s.get('submitter'), s.get('cls'), s.get('evaluated')] for s in r['scv']]
        bx = r.get('bx') or {}
        if bx.get('Clinical_significance_ENIGMA'):
            d['en'] = [bx['Clinical_significance_ENIGMA'], iso_date(str(bx.get('Date_last_evaluated_ENIGMA') or ''))]
        if c in t9:
            d['t9'] = t9[c]
        gn = r.get('gn')
        if gn:
            f, pop = faf_of(gn)
            d['gn'] = [sum((gn.get(k) or {}).get('ac') or 0 for k in ('exome', 'genome')), sum((gn.get(k) or {}).get('an') or 0 for k in ('exome', 'genome')), f, pop]
        lit = r.get('lit') or {}
        if lit.get('cv') or lit.get('lv'):
            d['lit'] = [lit.get('cv') or [], lit.get('lv') or []]
        if d:
            V[c] = d
    mave = {}
    for urn, m in (db.get('mave') or {}).items():
        sc = {}
        for row in csv.DictReader(io.StringIO(m.get('csv') or '')):
            try:
                score = float(row.get('score'))
            except (TypeError, ValueError):
                continue
            if score != score:  # NaN
                continue
            nt = (row.get('hgvs_nt') or '').split(':')[-1]
            key = nt if nt.startswith('c.') else (row.get('hgvs_pro') or '')
            if key and key != 'NA':
                sc[key] = round(score, 4)
        mave[urn] = {'mod': (m.get('meta') or {}).get('modificationDate'), 's': sc}
    src = meta.get('sources') or {}
    return {'built': meta['built'], 'builder': meta.get('builder'), 'sources': {k: (x or {}).get('version') for k, x in src.items()},
            'has': {'scv': bool((meta.get('counts') or {}).get('clinvar_submissions')), 'lit': 'literature' in src}, 'v': V, 'mave': mave}


def read_release_file(path):
    with open(path, encoding='utf-8') as f:
        t = f.read()
    i = t.index('window.BRCA1_DB = ') + len('window.BRCA1_DB = ')
    return json.loads(t[i:].rstrip().rstrip(';'))


def write_release(path, db):
    text = '/* BRCA1 variant explorer data release, built ' + db['meta']['built'] + ' by build_brca1_db.py. Public data; see meta.sources. */\nwindow.BRCA1_DB = ' + json.dumps(db, separators=(',', ':')) + ';\n'
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        f.write(text)
    os.replace(tmp, path)
    return len(text)


def digest_path(rdir, built):
    return os.path.join(rdir, f"brca1-release-{built.replace(':', '')}.json.gz")


def save_digest(rdir, dg):
    os.makedirs(rdir, exist_ok=True)
    p = digest_path(rdir, dg['built'])
    with gzip.open(p + '.tmp', 'wt', encoding='utf-8') as f:
        json.dump(dg, f, separators=(',', ':'))
    os.replace(p + '.tmp', p)
    return p


def load_digest(path):
    if path.endswith('.json.gz'):
        with gzip.open(path, 'rt', encoding='utf-8') as f:
            return json.load(f)
    if path.endswith('.json'):
        with open(path, encoding='utf-8') as f:
            return json.load(f)
    return release_digest(read_release_file(path))  # a brca1-data.js


def previous_digest(a, built, rdir):
    """--previous if given; otherwise the newest digest in --releases older than this build; otherwise the release file about to be replaced."""
    if a.previous:
        return load_digest(os.path.expanduser(a.previous))
    me = os.path.basename(digest_path(rdir, built))
    olds = sorted(n for n in (os.listdir(rdir) if os.path.isdir(rdir) else []) if re.match(r'brca1-release-.+\.json\.gz$', n) and n < me)
    if olds:
        return load_digest(os.path.join(rdir, olds[-1]))
    if not a.changes_only and os.path.exists(a.out):
        log(f'  reading the previous release {a.out}')
        dg = release_digest(read_release_file(a.out))
        if dg['built'] < built:
            save_digest(rdir, dg)
            return dg
    return None


def release_changes(prev, cur, keys):
    """Per variant: [kind, 1 if classification-level else 0, text]. keys = the variants in the new release."""
    P, C = prev.get('v') or {}, cur.get('v') or {}
    both = lambda k: (prev.get('has') or {}).get(k) and (cur.get('has') or {}).get(k)
    sge_a = ((prev.get('mave') or {}).get(SGE_URN) or {}).get('s') or {}
    sge_b = ((cur.get('mave') or {}).get(SGE_URN) or {}).get('s') or {}
    sge_cls = lambda x: 'non-functional' if x <= SGE_BINS[0] else 'functional' if x >= SGE_BINS[1] else 'intermediate'
    band = lambda f: 0 if f is None else sum(f > L for L in FAF_LEVELS)
    stars = lambda n: '★' * (n or 0) + '☆' * (4 - (n or 0))
    base = lambda acc: (acc or '').split('.')[0]
    sdesc = lambda s: f"{s[1]} ({s[2] or 'no classification'}{', ' + s[3] if s[3] else ''})"
    out = {}
    for c in sorted((set(P) | set(C) | set(sge_a) | set(sge_b)) & keys, key=natkey):
        a, b, items = P.get(c, {}), C.get(c, {}), []
        add = lambda kind, major, text: items.append([kind, 1 if major else 0, text])
        # ClinVar aggregate record
        ca, cb = a.get('cv'), b.get('cv')
        if cb and not ca:
            add('clinvar', True, f"New in ClinVar: {cb[0] or 'no classification'}{f' ({cb[1]})' if cb[1] else ''}")
        elif ca and not cb:
            add('clinvar', True, f'No ClinVar record in this release (was {ca[0]})')
        elif ca and cb:
            if (ca[0] or '') != (cb[0] or ''):
                add('clinvar', True, f'ClinVar: {ca[0]} → {cb[0]}')
            if (ca[1] or '') != (cb[1] or '') or ca[3] != cb[3]:
                add('review', ca[3] != cb[3], f"ClinVar review status: {ca[1] or '—'} → {cb[1] or '—'}" + (f' ({stars(ca[3])} → {stars(cb[3])})' if ca[3] != cb[3] else ''))
            if (ca[2] or '') != (cb[2] or '') and (ca[0] or '') == (cb[0] or ''):
                add('evaluated', False, f"ClinVar date last evaluated: {ca[2] or '—'} → {cb[2] or '—'}")
        # ClinVar submissions
        if both('scv'):
            A, B = {base(s[0]): s for s in a.get('scv') or []}, {base(s[0]): s for s in b.get('scv') or []}
            new, gone = [B[k] for k in B if k not in A], [A[k] for k in A if k not in B]
            if new:
                kinds = {tier(s[2]) for s in B.values()} - {tier(s[2]) for s in A.values()} - {None}  # a kind of classification no submitter gave before
                add('scv_new', bool(kinds), f"{len(new)} new ClinVar submission{'s' if len(new) > 1 else ''}: " + '; '.join(map(sdesc, new[:4])) + ('…' if len(new) > 4 else ''))
            for k in B:
                if k in A and (A[k][2] or '') != (B[k][2] or ''):
                    add('scv_changed', tier(A[k][2]) != tier(B[k][2]), f'ClinVar submission {k} ({B[k][1]}): {A[k][2]} → {B[k][2]}')
            if gone:
                add('scv_gone', False, f"{len(gone)} ClinVar submission{'s' if len(gone) > 1 else ''} no longer listed: " + '; '.join(map(sdesc, gone[:4])) + ('…' if len(gone) > 4 else ''))
        # ENIGMA (BRCA Exchange) and ENIGMA Table 9
        ea, eb = a.get('en'), b.get('en')
        if eb and not ea:
            add('enigma', True, f"New ENIGMA classification: {eb[0]}{f' ({eb[1]})' if eb[1] else ''}")
        elif ea and not eb:
            add('enigma', True, f'ENIGMA classification no longer in BRCA Exchange (was {ea[0]})')
        elif ea and eb:
            if ea[0] != eb[0]:
                add('enigma', True, f"ENIGMA: {ea[0]} → {eb[0]}{f' ({eb[1]})' if eb[1] else ''}")
            elif (ea[1] or '') != (eb[1] or ''):
                add('enigma_date', False, f"ENIGMA date evaluated: {ea[1] or '—'} → {eb[1] or '—'}")
        if a.get('t9') != b.get('t9'):
            f9 = lambda t: ' '.join(str(x) for x in t if x) if t else 'not listed'
            add('t9', True, f"ENIGMA Table 9: {f9(a.get('t9'))} → {f9(b.get('t9'))}")
        # gnomAD
        ga, gb = a.get('gn'), b.get('gn')
        if gb and not ga:
            add('gnomad', band(gb[2]) > 0, f'Now in gnomAD: AC {gb[0]:,}/{gb[1]:,}' + (f", FAF95 {sci(gb[2])} {(gb[3] or '').upper()}" if gb[2] is not None else ''))
        elif ga and not gb:
            add('gnomad', False, f'No longer in gnomAD (was AC {ga[0]:,})')
        elif ga and gb:
            if ga[0] != gb[0]:
                add('gnomad_ac', False, f'gnomAD AC {ga[0]:,} → {gb[0]:,} (AN {gb[1]:,})')
            ba, bb = band(ga[2]), band(gb[2])
            if ba != bb:
                lvl = FAF_LEVELS[bb - 1] if bb > ba else FAF_LEVELS[bb]
                add('faf', True, f"gnomAD FAF95 {sci(ga[2])} → {sci(gb[2])} {(gb[3] or '').upper()}: now {'above' if bb > ba else 'below'} {sci(lvl)}")
        # SGE (Findlay 2018)
        sa, sb = sge_a.get(c), sge_b.get(c)
        if sb is not None and sa is None:
            add('sge', True, f'New SGE score (Findlay 2018): {sb:.3f}, {sge_cls(sb)}')
        elif sa is not None and sb is None:
            add('sge', True, f'SGE score no longer in MaveDB (was {sa:.3f})')
        elif sa is not None and (abs(sa - sb) >= 0.05 or sge_cls(sa) != sge_cls(sb)):
            add('sge', sge_cls(sa) != sge_cls(sb), f'SGE score (Findlay 2018) {sa:.3f} → {sb:.3f}' + (f': {sge_cls(sa)} → {sge_cls(sb)}' if sge_cls(sa) != sge_cls(sb) else ''))
        # papers
        if both('lit'):
            la, lb = a.get('lit') or [[], []], b.get('lit') or [[], []]
            new = sorted((set(lb[0]) | set(lb[1])) - set(la[0]) - set(la[1]))
            new_cv = set(lb[0]) - set(la[0])
            if new:
                add('lit', False, f"{len(new)} new paper{'s' if len(new) > 1 else ''} indexed" + (f' ({len(new_cv)} cited in ClinVar for this allele)' if new_cv else '') + f": PMID {', '.join(map(str, new[:5]))}{'…' if len(new) > 5 else ''}")
        if items:
            out[c] = items
    # source versions and MaveDB score sets (set level)
    sp, sc_ = prev.get('sources') or {}, cur.get('sources') or {}
    sources = [[SRC_NAMES.get(k, k), sp.get(k), sc_.get(k)] for k in sorted(set(sp) | set(sc_)) if sp.get(k) != sc_.get(k)]
    assays, mp, mc = [], prev.get('mave') or {}, cur.get('mave') or {}
    for urn in sorted(set(mp) | set(mc)):
        ma, mb, short = mp.get(urn), mc.get(urn), urn.split(':')[-1]
        if not ma:
            assays.append(f"MaveDB {short}: new score set ({len(mb['s']):,} scores)")
        elif not mb:
            assays.append(f'MaveDB {short}: not in this release')
        else:
            n_ch = sum(1 for k, v in mb['s'].items() if k in ma['s'] and abs(v - ma['s'][k]) >= 0.05)
            n_add, n_rem = len(set(mb['s']) - set(ma['s'])), len(set(ma['s']) - set(mb['s']))
            if n_ch or n_add or n_rem or ma.get('mod') != mb.get('mod'):
                assays.append(f"MaveDB {short} (modified {ma.get('mod') or '—'} → {mb.get('mod') or '—'}): scores changed by ≥ 0.05: {n_ch}, added: {n_add}, removed: {n_rem}")
    return {'from': prev['built'], 'to': cur['built'], 'sources': sources, 'assays': assays, 'v': out,
            'n': {'variants': len(out), 'major': sum(1 for x in out.values() if any(i[1] for i in x))}}


def read_reported(path, db):
    """Your laboratory's reported variants (CSV/TSV with an hgvsc column; optional classification, date, sample/report), matched to the release."""
    with open(path, encoding='utf-8-sig', newline='') as f:
        text = f.read()
    first = text.splitlines()[0] if text.strip() else ''
    delim = '\t' if '\t' in first else ';' if first.count(';') > first.count(',') else ','
    idx_c = (db.get('idx') or {}).get('c') or {}
    out = []
    for r in csv.DictReader(io.StringIO(text), delimiter=delim):
        r = {(k or '').strip().lower(): (v or '').strip() for k, v in r.items() if k}
        if re.search(r'artefact|artifact|freq', r.get('type', ''), re.I):
            continue  # the page's in-house CSV can be used as is
        given = r.get('hgvsc') or r.get('variant') or r.get('hgvs') or r.get('c.') or ''
        if not given:
            continue
        n = norm_c(given)
        key = n if n in db['v'] else idx_c.get(n) or idx_c.get(given.split(':')[-1])
        out.append({'given': given, 'key': key, 'cls': r.get('classification') or r.get('class') or '', 'date': r.get('date') or r.get('reported') or r.get('report_date') or '',
                    'ref': r.get('sample') or r.get('report') or r.get('report_id') or r.get('id') or ''})
    return out


def write_change_reports(rdir, db, ch, cur, reported):
    os.makedirs(rdir, exist_ok=True)
    prot = lambda c: (((db['v'].get(c) or {}).get('vep') or {}).get('transcript_consequences') or [{}])[0].get('hgvsp', '').split(':')[-1]
    rep_keys = {r['key'] for r in reported or [] if r['key']}
    order = sorted(ch['v'].items(), key=lambda kv: (-max(i[1] for i in kv[1]), natkey(kv[0])))
    p1 = os.path.join(rdir, f"brca1-changes_{ch['from'][:10]}_to_{ch['to'][:10]}.csv")
    with open(p1, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(['variant', 'protein', 'classification_level', 'kind', 'change'] + (['reported_by_you'] if reported is not None else []))
        for c, items in order:
            for i in items:
                w.writerow([c, prot(c), 'yes' if i[1] else '', i[0], i[2]] + (['yes' if c in rep_keys else ''] if reported is not None else []))
    out = {'changes': p1}
    if reported is None:
        return out
    p2 = os.path.join(rdir, f"brca1-reclassification-queue_{ch['to'][:10]}.csv")
    rows = []
    for r in reported:
        items = ch['v'].get(r['key']) if r['key'] else None
        if r['key'] and not items:
            continue
        d = (cur.get('v') or {}).get(r['key']) or {}
        cv, en, t9 = d.get('cv'), d.get('en'), d.get('t9')
        expert = ('ENIGMA', en[0]) if en else ('ClinVar', cv[0]) if cv and (cv[3] or 0) >= 2 else None
        vs = '' if not (r['cls'] and expert and tier(r['cls']) and tier(expert[1])) else f"{'same as' if tier(r['cls']) == tier(expert[1]) else 'differs from'} {expert[0]} ({expert[1]})"
        rows.append([r['given'], r['key'] or '', prot(r['key']) if r['key'] else '', r['cls'], r['date'], r['ref'],
                     'yes' if items and any(i[1] for i in items) else '',
                     ' | '.join(i[2] for i in items or [] if i[1]), ' | '.join(i[2] for i in items or [] if not i[1]) if items else 'not found in the release: check the HGVS (exon-level CNVs are not in the release)',
                     f"{cv[0]} ({'★' * (cv[3] or 0)})" if cv else '', en[0] if en else '', ' '.join(str(x) for x in t9 if x) if t9 else '', vs])
    rows.sort(key=lambda x: (x[1] == '', x[6] != 'yes', natkey(x[1] or x[0])))
    with open(p2, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(['variant_as_reported', 'release_key', 'protein', 'your_classification', 'your_date', 'sample_or_report', 'classification_level_change',
                    'classification_level_changes', 'other_changes', 'clinvar_now', 'enigma_now', 'enigma_table9_now', 'your_class_vs_expert_now'])
        w.writerows(rows)
    out.update(queue=p2, queued=sum(1 for x in rows if x[1]), unmatched=sum(1 for x in rows if not x[1]))
    return out


def apply_changes(a, db):
    """Compares this release with the previous one, stores the result in db['changes'], writes the CSV files and keeps this release's digest."""
    rdir = os.path.expanduser(a.releases or os.path.join(os.path.dirname(os.path.abspath(a.out)), 'brca1-releases'))
    log('Changes since the previous release')
    cur = release_digest(db)
    try:
        prev = previous_digest(a, cur['built'], rdir)
    except Exception as e:
        log(f'  WARNING: previous release could not be read ({e}); no change list')
        prev = None
    reported = None
    if a.reported:
        try:
            reported = read_reported(os.path.expanduser(a.reported), db)
            log(f"  {len(reported)} reported variants read from {a.reported} ({sum(1 for r in reported if not r['key'])} not matched to the release)")
        except Exception as e:
            log(f'  WARNING: could not read {a.reported} ({e})')
    if prev and prev['built'] < cur['built']:
        ch = release_changes(prev, cur, set(db['v']))
        db['changes'] = ch
        log(f"  release {cur['built'][:10]} vs {prev['built'][:10]}: {ch['n']['variants']} variants changed ({ch['n']['major']} at classification level); "
            f"{len(ch['sources'])} source versions changed")
        files = write_change_reports(rdir, db, ch, cur, reported)
        log(f"  wrote {files['changes']}")
        if 'queue' in files:
            log(f"  reclassification queue: {files['queued']} of your reported variants changed, {files['unmatched']} not found -> {files['queue']}")
    else:
        db['changes'] = {'from': None}
        log('  no earlier release to compare with: the next build will list what changed since this one')
    log(f'  kept this release\'s digest: {save_digest(rdir, cur)}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'brca1-data.js'))
    ap.add_argument('--cache', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '.brca1-build-cache'))
    ap.add_argument('--max-age', type=float, default=3, help='refetch ClinVar, BRCA Exchange, MaveDB, ENIGMA and UniProt data cached more than this many days ago (default 3)')
    ap.add_argument('--refresh', nargs='*', default=[], help='force groups to be fetched again: ensembl vep gnomad gnomad-details clinvar scv bx mavedb enigma uniprot alphafold splicevault alfa literature pubmed, or all')
    ap.add_argument('--skip-literature', action='store_true', help='skip ClinVar citations, LitVar2 and PubMed titles (no paper list in the page)')
    ap.add_argument('--lit-max-age', type=float, default=30, help='refetch ClinVar citations and LitVar2 results cached more than this many days ago (default 30)')
    ap.add_argument('--skip-scv', action='store_true', help='skip the ClinVar submission_summary download (no per-submission table)')
    ap.add_argument('--ncbi-api-key', default=os.environ.get('NCBI_API_KEY'), help='optional NCBI API key (faster E-utilities)')
    ap.add_argument('--releases', default=None, help='folder for release digests and the change and reclassification-queue CSV files (default: brca1-releases next to --out)')
    ap.add_argument('--previous', default=None, help='release to compare with: a brca1-data.js or a digest from --releases (default: the newest earlier release)')
    ap.add_argument('--reported', default=None, help="CSV of your laboratory's reported variants (hgvsc; optional classification, date, sample): writes the reclassification queue")
    ap.add_argument('--changes-only', action='store_true', help='download nothing: recompute the change list of the existing --out release against --previous (or the newest earlier digest)')
    ap.add_argument('--max-snv-positions', type=int, default=0, help=argparse.SUPPRESS)  # testing only
    a = ap.parse_args()
    if a.changes_only:
        log(f'Reading {a.out}')
        db = read_release_file(a.out)
        db.pop('changes', None)
        apply_changes(a, db)
        n = write_release(a.out, db)
        log(f'Wrote {a.out} ({n / 1e6:.1f} MB) with the change list')
        return
    if a.ncbi_api_key:
        PACE['eutils.ncbi.nlm.nih.gov'] = 0.11
    cache = Cache(a.cache, a.refresh, a.max_age)
    t0 = time.time()
    src = {}

    def stamp(key, version):
        src[key] = {'version': version, 'retrieved': now_iso()}

    # 1. Ensembl transcript, sequences, release
    log('Ensembl: release, transcript, sequences')
    info = cache.get('ensembl-info', lambda: get_json('https://rest.ensembl.org/info/software?content-type=application/json'), group='ensembl', volatile=True)
    release = info.get('release')
    lk = cache.get(f'lookup-{release}', lambda: get_json(f'https://rest.ensembl.org/lookup/id/{ENST}?expand=1;content-type=application/json'), group='ensembl')
    cds = cache.get(f'cds-{release}', lambda: http(f'https://rest.ensembl.org/sequence/id/{ENST}?type=cds;content-type=text/plain').strip(), group='ensembl')
    g0, g1 = lk['start'] - 3000, lk['end'] + 3000
    seq = cache.get(f'region-{release}-{g0}-{g1}', lambda: http(f'https://rest.ensembl.org/sequence/region/human/{CHROM}:{g0}..{g1}:1?content-type=text/plain').strip().upper(), group='ensembl')
    tx = build_tx(lk)
    stamp('vep', f'Ensembl release {release} (VEP REST, GRCh38)')
    stamp('compara', f'Ensembl release {release}')

    # 2. UniProt, AlphaFold
    log('UniProt and AlphaFold')
    fasta = cache.get('uniprot-fasta', lambda: http(f'https://rest.uniprot.org/uniprotkb/{UNIPROT}.fasta'), group='uniprot', volatile=True)
    up = cache.get('uniprot-features', lambda: get_json(f'https://rest.uniprot.org/uniprotkb/{UNIPROT}?fields=ft_domain,ft_region,ft_zn_fing,ft_coiled,ft_motif,ft_mutagen,ft_mod_res&format=json'), group='uniprot', volatile=True)
    up_ver = (up.get('entryAudit') or {}).get('entryVersion')
    up = {'features': up.get('features', [])}
    try:
        af = cache.get('alphafold', lambda: get_json(f'https://alphafold.ebi.ac.uk/files/AF-{UNIPROT}-F1-confidence_v6.json'), group='alphafold', volatile=True)
        af = {'confidenceScore': af.get('confidenceScore')} if isinstance(af, dict) else None
    except Exception as e:
        log('  AlphaFold unavailable:', e); af = None
    stamp('uniprot', f'{UNIPROT}' + (f', entry version {up_ver}' if up_ver else ''))
    stamp('alphafold', 'v6')

    # 3. gnomAD gene variants and ClinVar variants (the same query the page makes)
    log('gnomAD: gene variants and ClinVar variants')
    gq = f'''{{ gene(gene_symbol: "{GENE}", reference_genome: GRCh38) {{
      variants(dataset: gnomad_r4) {{ variant_id pos ref alt consequence hgvsc hgvsp transcript_id exome {{ ac an homozygote_count }} genome {{ ac an homozygote_count }} }}
      clinvar_variants {{ variant_id pos ref alt clinical_significance gold_stars hgvsc hgvsp major_consequence review_status clinvar_variation_id }} }} }}'''
    gene = cache.get('gnomad-gene', lambda: gnomad(gq), group='gnomad', volatile=True)
    gvars = [v for v in gene['data']['gene']['variants'] if v.get('transcript_id') in (ENST, None)]
    cvars = gene['data']['gene']['clinvar_variants']
    gene['data']['gene']['variants'] = gvars
    log(f'  {len(gvars)} gnomAD variants, {len(cvars)} ClinVar variants')

    # 4. gnomAD per-variant details (filters, flags, FAF, ancestry groups), in batches
    log('gnomAD: per-variant details (paced; ~6 s per batch)')
    VQ = 'flags exome { ac an homozygote_count filters faf95 { popmax popmax_population } populations { id ac an } } genome { ac an homozygote_count filters faf95 { popmax popmax_population } populations { id ac an } } joint { ac an homozygote_count faf95 { popmax popmax_population } }'
    gstore = cache.store('gnomad-v4.1-variants', 'gnomad-details')
    gnomad_errors = []

    def fetch_details():
        try:
            todo = [v['variant_id'] for v in gvars if v['variant_id'] not in gstore]
            log(f'  {len(todo)} to fetch ({len(gvars) - len(todo)} cached)')
            B = 20
            for i in range(0, len(todo), B):
                chunk = todo[i:i + B]
                q = '{ ' + ' '.join(f'v{k}: variant(variantId: "{vid}", dataset: gnomad_r4) {{ {VQ} }}' for k, vid in enumerate(chunk)) + ' }'
                d = post_json('https://gnomad.broadinstitute.org/api', {'query': q}, timeout=240).get('data') or {}
                got = {}
                for k, vid in enumerate(chunk):
                    x = d.get(f'v{k}')
                    if x:
                        for s_ in ('exome', 'genome'):
                            if x.get(s_) and x[s_].get('populations'):
                                x[s_]['populations'] = [p for p in x[s_]['populations'] if '_' not in p['id'] and p['id'] not in ('XX', 'XY')]
                    got[vid] = x
                gstore.add(got)
                if (i // B + 1) % 25 == 0:
                    log(f'  {min(i + B, len(todo))}/{len(todo)}')
        except Exception as e:  # reported at the end; the release is still written
            gnomad_errors.append(e)
    gthread = threading.Thread(target=fetch_details, daemon=True)
    gthread.start()
    log('  (continuing with other sources while gnomAD details download)')

    def details_now():
        gthread.join()
        if gnomad_errors:
            log(f'  WARNING: gnomAD details incomplete ({gnomad_errors[0]}); rerun to complete them')
        return {v['variant_id']: gstore.get(v['variant_id']) for v in gvars}
    # read depth around every exon
    log('gnomAD: coverage around exons')
    cov = {'exome': {}, 'genome': {}}
    for e in tx['ex']:
        a_, b_ = e['start'] - 60, e['end'] + 60
        c = cache.get(f'gnomad-cov-{a_}-{b_}', lambda a_=a_, b_=b_: gnomad(f'{{ region(chrom: "{CHROM}", start: {a_}, stop: {b_}, reference_genome: GRCh38) {{ coverage(dataset: gnomad_r4) {{ exome {{ pos mean median over_20 }} genome {{ pos mean median over_20 }} }} }} }}'), group='gnomad')
        cc = ((c.get('data') or {}).get('region') or {}).get('coverage') or {}
        for k in ('exome', 'genome'):
            for x in cc.get(k) or []:
                cov[k][x['pos']] = [x['pos'], round(x['mean'], 1) if x.get('mean') is not None else None, x.get('median'), round(x['over_20'], 3) if x.get('over_20') is not None else None]
    cov = {k: [v[p] for p in sorted(v)] for k, v in cov.items()}
    stamp('gnomad', 'v4.1 (exomes + genomes)')

    # 4b. SpliceVault (Ensembl FTP) and NCBI ALFA: tabix range requests, so only the BRCA1 region of these genome-wide files is downloaded
    log('SpliceVault and ALFA (BRCA1 region only)')
    sv = alfa = None
    try:
        sv = cache.get(f'splicevault-{g0}-{g1}', lambda: splicevault_region(g0, g1), group='splicevault')  # a fixed dataset: kept
        stamp('splicevault', f"SpliceVault 300K-RNA (Ensembl VEP plugin file), {len(sv['v'])} BRCA1 variants")
        log(f"  SpliceVault: {len(sv['v'])} variants")
    except Exception as e:
        log(f'  WARNING: SpliceVault unavailable ({e})')
    try:
        alfa = cache.get(f'alfa-{g0}-{g1}', lambda: alfa_region(g0, g1), group='alfa', volatile=True)
        stamp('alfa', f"ALFA {alfa.get('release') or ''}".strip() + f", {len(alfa['v'])} BRCA1 variants")
        log(f"  ALFA: {len(alfa['v'])} variants")
    except Exception as e:
        log(f'  WARNING: ALFA unavailable ({e})')

    # 5. ClinVar aggregate records (E-utilities esummary) and submissions (submission_summary)
    log('ClinVar: aggregate records')
    vids = sorted({str(v['clinvar_variation_id']) for v in cvars if v.get('clinvar_variation_id')}, key=int)
    key = f'&api_key={a.ncbi_api_key}' if a.ncbi_api_key else ''
    cv = {}
    for i in range(0, len(vids), 300):
        chunk = vids[i:i + 300]

        def fetch(chunk=chunk):
            body = urllib.parse.urlencode({'db': 'clinvar', 'id': ','.join(chunk), 'retmode': 'json'}).encode()
            j = json.loads(http('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?tool=brca1-explorer' + key, data=body, headers={'Content-Type': 'application/x-www-form-urlencoded'}))
            out = {}
            for vid in chunk:
                r = (j.get('result') or {}).get(vid)
                if not r:
                    continue
                gc = r.get('germline_classification') or {}
                le = gc.get('last_evaluated')
                out[vid] = {'cv': {'gc': {'description': gc.get('description') or None, 'last_evaluated': None if not le or le.startswith('1/01/01') else le, 'review_status': gc.get('review_status'),
                                          'trait_set': [{'trait_name': t.get('trait_name')} for t in gc.get('trait_set') or []][:6]},
                                   'nscv': len((r.get('supporting_submissions') or {}).get('scv') or [])},
                            'vcv': {'acc': r.get('accession'), 'version': r.get('accession_version', '').split('.')[-1] if r.get('accession_version') else None}}
            return out
        cv.update(cache.get(f'clinvar-esummary-{hashlib.sha1(",".join(chunk).encode()).hexdigest()[:16]}', fetch, group='clinvar', volatile=True))
    scv, clinvar_release = {}, None
    if not a.skip_scv:
        log('ClinVar: submission_summary.txt.gz (about 100 MB, streamed and filtered)')

        def fetch_scv():
            want, out, hdr, meta = set(vids), {}, None, {}
            req = urllib.request.Request('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz', headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=600, context=SSL_CTX) as r:
                meta['last-modified'] = r.headers.get('Last-Modified')
                for raw in io.TextIOWrapper(gzip.GzipFile(fileobj=r), encoding='utf-8', errors='replace'):
                    if raw.startswith('#VariationID'):
                        hdr = raw.lstrip('#').rstrip('\n').split('\t'); continue
                    if raw.startswith('#') or not hdr:
                        continue
                    f = raw.rstrip('\n').split('\t')
                    if f[0] not in want:
                        continue
                    d = dict(zip(hdr, f))
                    cm = d.get('Description', '')
                    out.setdefault(f[0], []).append({'scv': d.get('SCV'), 'submitter': d.get('Submitter'), 'cls': d.get('ClinicalSignificance'),
                        'evaluated': iso_date(d.get('DateLastEvaluated')), 'review': d.get('ReviewStatus'), 'comment': '' if cm in ('-', 'na') else cm[:1500],
                        'origin': d.get('OriginCounts', ''), 'method': d.get('CollectionMethod', '')})
            if hdr is None:
                raise RuntimeError('submission_summary: header line not found')
            for k in out:
                out[k].sort(key=lambda s: s['evaluated'] or '', reverse=True)
            return {'scv': out, 'modified': meta.get('last-modified')}
        try:
            res = cache.get('clinvar-scv', fetch_scv, group='scv', volatile=True)
            scv, clinvar_release = res['scv'], res.get('modified')
            log(f'  {sum(len(v) for v in scv.values())} submissions for {len(scv)} variants')
        except Exception as e:
            log(f'  WARNING: ClinVar submissions unavailable ({e}); the release will have no per-submission table')
    stamp('clinvar', f'ClinVar via gnomAD + E-utilities' + (f'; submission_summary {clinvar_release}' if clinvar_release else ''))

    # 6. BRCA Exchange: ENIGMA classifications
    log('BRCA Exchange: ENIGMA classifications')

    def fetch_bx():
        rows, page = [], 0
        while True:
            q = urllib.parse.urlencode([('format', 'json'), ('filter', 'Gene_Symbol'), ('filterValue', GENE), ('include', 'Variant_in_ENIGMA'), ('page_size', 1000), ('page_num', page)])
            j = get_json('https://brcaexchange.org/backend/data/?' + q, timeout=240)
            data = j.get('data') or []
            rows += [{k: d.get(k) for k in BX_KEEP if d.get(k) not in (None, '', '-')} for d in data if d.get('Gene_Symbol') == GENE]
            page += 1
            if not data or len(rows) >= int(j.get('count') or 0) or page > 100:
                return rows
    bx_rows = cache.get('bx-enigma', fetch_bx, group='bx', volatile=True)
    log(f'  {len(bx_rows)} BRCA1 rows with an ENIGMA record')
    stamp('bx', f'BRCA Exchange API, {len(bx_rows)} BRCA1 ENIGMA records')

    # 7. MaveDB score sets (meta + scores, trimmed to the columns the page reads)
    log('MaveDB score sets')
    mave = {}
    for urn in MAVE:
        def fetch(urn=urn):
            meta = get_json(f'https://api.mavedb.org/api/v1/score-sets/{urn}')
            text = http(f'https://api.mavedb.org/api/v1/score-sets/{urn}/scores', timeout=300)
            rd = csv.reader(io.StringIO(text))
            head = next(rd)
            keep = [c for c in ('accession', 'hgvs_nt', 'hgvs_splice', 'hgvs_pro', 'score') if c in head]
            ix = [head.index(c) for c in keep]
            buf = io.StringIO()
            w = csv.writer(buf, lineterminator='\n')
            w.writerow(keep)
            for row in rd:
                if row:
                    w.writerow([row[i] if i < len(row) else '' for i in ix])
            tg = (meta.get('targetGenes') or [{}])[0]
            slim = {'urn': meta.get('urn'), 'title': meta.get('title'), 'publishedDate': meta.get('publishedDate'), 'modificationDate': meta.get('modificationDate'),
                    'license': {'shortName': (meta.get('license') or {}).get('shortName')}, 'targetGenes': [{'targetSequence': tg.get('targetSequence')}],
                    'scoreCalibrations': meta.get('scoreCalibrations')}
            return {'meta': slim, 'csv': buf.getvalue()}
        try:
            mave[urn] = cache.get(f'mave-{urn}', fetch, group='mavedb', volatile=True)
        except Exception as e:
            log(f'  {urn} unavailable: {e}')
    stamp('mavedb', '; '.join(f"{u.split(':')[-1]} ({(m['meta'].get('modificationDate') or m['meta'].get('publishedDate') or '')[:10]})" for u, m in mave.items()))

    # 8. ENIGMA VCEP workbooks from the ClinGen CSpec registry
    log('ENIGMA VCEP v1.2 workbooks (ClinGen CSpec)')
    enigma, hashes = {}, {}
    for k, fid in ENIGMA_FILES.items():
        def fetch(fid=fid, k=k):
            data = http(f'https://cspec.genome.network/cspec/File/id/{fid}/data', binary=True, timeout=300)
            return {'sheets': read_xlsx(data, ENIGMA_SHEETS[k]), 'sha256': hashlib.sha256(data).hexdigest()}
        try:
            r = cache.get(f'enigma-{k}', fetch, group='enigma', volatile=True)
            enigma[k], hashes[k] = r['sheets'], r['sha256'][:12]
        except Exception as e:
            if k == 'supp':
                log('  supplementary tables unavailable:', e)
            else:
                raise
    stamp('enigma', 'v1.2, 2024-11-18 (Tables 4, 9; Suppl. Tables 2, 4, 7, 10, 16); files ' + ', '.join(f'{k} sha256 {h}' for k, h in hashes.items()))

    # 9. Compara orthologs, projected onto human positions (as the page does)
    log('Ensembl Compara orthologs')
    protein = ''.join(fasta.split('\n')[1:]).strip()

    def fetch_hom():
        u = (f'https://rest.ensembl.org/homology/symbol/human/{GENE}?type=orthologues;aligned=1;sequence=protein;cigar_line=0;content-type=application/json;'
             + ';'.join('target_species=' + s for s in SPECIES))
        j = get_json(u, timeout=300)
        out = {}
        for hm in j['data'][0]['homologies']:
            s, t = hm['source']['align_seq'], hm['target']['align_seq']
            proj, i = ['-'] * len(protein), 0
            for k in range(len(s)):
                if s[k] == '-':
                    continue
                if i < len(protein):
                    proj[i] = '-' if t[k] == '-' else t[k]
                i += 1
            sp = hm['target']['species']
            if sp not in out or hm['type'] == 'ortholog_one2one':
                out[sp] = {'seq': ''.join(proj), 'pid': hm['source']['perc_id'], 'id': hm['target']['protein_id']}
        return out
    homology = cache.get(f'compara-{release}', fetch_hom, group='ensembl')

    # 10. variant set
    log('Assembling the variant set')
    region_inputs, hgvs_inputs = {}, set()

    def add_vcf(pos, ref, alt):
        region_inputs[f'{CHROM}-{pos}-{ref}-{alt}'] = f'{CHROM} {pos} . {ref} {alt} . . .'
    for v in gvars + cvars:
        add_vcf(v['pos'], v['ref'], v['alt'])
    # every coding SNV and every SNV within INTRON_FLANK bases of an exon
    positions = set()
    for c in range(1, tx['cds_len'] + 1):
        g = tx['c_to_g'](f'c.{c}')
        if g:
            positions.add(g)
    for i, e in enumerate(tx['ex']):
        for k in range(1, INTRON_FLANK + 1):
            if i > 0:
                positions.add(e['end'] + k if tx['minus'] else e['start'] - k)
            if i < len(tx['ex']) - 1:
                positions.add(e['start'] - k if tx['minus'] else e['end'] + k)
    if a.max_snv_positions:
        positions = set(sorted(positions)[:a.max_snv_positions])
    for g in positions:
        ref = seq[g - g0] if 0 <= g - g0 < len(seq) else None
        if ref in ('A', 'C', 'G', 'T'):
            for alt in 'ACGT':
                if alt != ref:
                    add_vcf(g, ref, alt)
    # named variants: SGE, ENIGMA tables, BRCA Exchange, founders
    named = set()
    if SGE_URN in mave:
        for row in csv.DictReader(io.StringIO(mave[SGE_URN]['csv'])):
            c = (row.get('hgvs_nt') or '').split(':')[-1]
            if c.startswith('c.'):
                named.add(c)
    t9 = next((rows for n, rows in (enigma.get('table9') or {}).items() if re.search('Table9', n)), [])
    for r in t9:
        if r and r[0] == GENE and len(r) > 1 and isinstance(r[1], str) and r[1].strip().startswith('c.'):
            named.add(r[1].strip())
    for n, rows in (enigma.get('table4') or {}).items():
        for r in rows:
            if not r or r[0] != GENE:
                continue
            for ip, ia in ((4, 5), (15, 16)):
                pre, alts = (r[ip] if len(r) > ip else None), (r[ia] if len(r) > ia else None)
                if isinstance(pre, str) and isinstance(alts, str) and re.match(r'^c\.[-\d]+[+-][12][ACGT]>$', pre.strip()):
                    for alt in re.split(r'[ ,]+', alts.strip()):
                        if alt:
                            named.add(pre.strip() + alt)
    full = re.compile(r'^c\.[-*]?\d+(?:[+-]\d+)?(?:_[-*]?\d+(?:[+-]\d+)?)?(?:[ACGT]>[ACGT]|delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|inv)$')
    for n, rows in (enigma.get('supp') or {}).items():
        for r in rows:
            if not r:
                continue
            if n.startswith('ST16') and r[0] == GENE and len(r) > 2 and isinstance(r[1], str) and isinstance(r[2], str):
                for alt in re.split(r'[ ,]+', r[2].strip()):
                    if full.match(r[1].strip() + alt):
                        named.add(r[1].strip() + alt)
                continue
            for x in r[:4]:
                if isinstance(x, str) and (GENE in r[:3] or n.startswith('ST10')):
                    c = re.sub(r'^.*:(?=c\.)', '', x.strip())
                    if full.match(c):
                        named.add(c)
    for d in bx_rows:
        c = str(d.get('HGVS_cDNA', '')).split(':')[-1]
        if c.startswith('c.'):
            named.add(c)
    named |= set(FC_BRCA1)
    # named SNVs at positions already in the region set need no HGVS call
    for c in sorted(named):
        m = re.match(r'^c\.([-*]?\d+(?:[+-]\d+)?)([ACGT])>([ACGT])$', c)
        if m:
            g = tx['c_to_g']('c.' + m.group(1))
            if g and 0 <= g - g0 < len(seq):
                ref, alt = (m.group(2), m.group(3)) if not tx['minus'] else (m.group(2).translate(COMP), m.group(3).translate(COMP))
                if seq[g - g0] == ref:
                    add_vcf(g, ref, alt)
                    continue
        hgvs_inputs.add(f'{REFSEQ}:{c}')
    log(f'  {len(region_inputs)} genomic inputs, {len(hgvs_inputs)} HGVS inputs')

    # 11. VEP (cached per Ensembl release)
    log(f'Ensembl VEP: annotating (cached per release {release}; first run is the slow part)')
    vstore = cache.store(f'vep-release-{release}', 'vep')
    reg = vep_annotate(sorted(region_inputs.values()), 'region', vstore)
    hg = vep_annotate(sorted(hgvs_inputs), 'hgvs', vstore)

    # 12. records
    log('Waiting for gnomAD details' if gthread.is_alive() else 'gnomAD details done')
    details = details_now()
    log('Assembling records')
    g_by_c = {v['hgvsc']: v['variant_id'] for v in gvars if v.get('hgvsc')}
    g_by_id = {v['variant_id'] for v in gvars}
    c_by_c = {v['hgvsc']: str(v['clinvar_variation_id']) for v in cvars if v.get('hgvsc') and v.get('clinvar_variation_id')}
    c_by_id = {v['variant_id']: str(v['clinvar_variation_id']) for v in cvars if v.get('clinvar_variation_id')}
    bx = {}
    for d in bx_rows:
        c = str(d.get('HGVS_cDNA', '')).split(':')[-1]
        bx.setdefault(c, d); bx.setdefault(norm_c(c), d)
    V, idx = {}, {'g': {}, 'rs': {}, 'p': {}, 'c': {}}
    fc_keys, vid_of = {}, {}
    region_set = set(region_inputs.values())
    for inp, s in list(reg.items()) + list(hg.items()):
        if not s:
            continue
        t = s['transcript_consequences'][0]
        c = (t.get('hgvsc') or '').split(':')[-1]
        if not c.startswith('c.'):
            continue
        if c not in V:
            gid = gid_of(s)
            rec = {'vep': s}
            vcf = None
            if inp in region_set:
                p = inp.split()
                vcf = f'{p[0]}-{p[1]}-{p[3]}-{p[4]}'
            gv = g_by_c.get(c) or (gid if gid in g_by_id else None) or (vcf if vcf in g_by_id else None)
            if gv and details.get(gv):
                rec['gn'] = details[gv]
            vid = c_by_c.get(c) or c_by_id.get(gid or '') or c_by_id.get(vcf or '')
            if vid:
                vid_of[c] = vid
            if vid and vid in cv:
                rec['cv'], rec['vcv'] = cv[vid]['cv'], cv[vid]['vcv']
            if vid and vid in scv:
                rec['scv'] = scv[vid]
            if c in bx or norm_c(c) in bx:
                rec['bx'] = bx.get(c) or bx.get(norm_c(c))
            V[c] = rec
            for key in filter(None, [gid, vcf]):
                idx['g'].setdefault(key, c)
            for cv_ in s.get('colocated_variants') or []:
                idx['rs'].setdefault(cv_['id'].lower(), c)
            hp = t.get('hgvsp', '')
            if ':p.' in hp:
                idx['p'].setdefault(hp.split(':p.')[1], [])
                if c not in idx['p'][hp.split(':p.')[1]]:
                    idx['p'][hp.split(':p.')[1]].append(c)
        if inp and inp.startswith(REFSEQ + ':'):
            given = inp.split(':', 1)[1]
            if given != c:
                idx['c'][given] = c
            if given in FC_BRCA1:
                fc_keys[c] = given
    for c in FC_BRCA1:
        if c in V:
            fc_keys.setdefault(c, c)
    log(f'  {len(V)} variant records')

    # 13. literature: ClinVar citations, LitVar2 papers per rsID, PubMed summaries (titles are indexed, not summarised)
    pubs, lit_note = {}, None
    if not a.skip_literature:
        try:
            pubs, lit_note = literature(a, V, vid_of, key)
            stamp('literature', lit_note)
        except Exception as e:
            log(f'  WARNING: literature unavailable ({e}); the release will have no paper list')

    db = {
        'meta': {'built': now_iso(), 'builder': f'build_brca1_db.py {VERSION}', 'transcript': f'{REFSEQ} / {ENST}', 'sources': src,
                 'counts': {'variants': len(V), 'gnomad': len(gvars), 'clinvar': len(cvars), 'clinvar_submissions': sum(len(x) for x in scv.values()), 'bx_enigma': len(bx_rows), 'papers': len(pubs)},
                 'seconds': round(time.time() - t0)},
        'gene': {'lookup': lk, 'cds': cds, 'fasta': fasta, 'uniprot': up, 'alphafold': af, 'gnomadGene': gene, 'homology': homology},
        'mave': mave, 'enigma': enigma, 'seq': {'from': g0, 'seq': seq}, 'cov': cov, 'sv': sv, 'alfa': alfa, 'fcKeys': fc_keys, 'idx': idx, 'v': V, 'pubs': pubs,
    }
    # 14. what changed since the previous release (read before the old file is replaced), and the reclassification queue
    try:
        apply_changes(a, db)
    except Exception as e:  # the release itself is still written
        log(f'  WARNING: change list failed ({e})')
        db['changes'] = {'from': None}
    n = write_release(a.out, db)
    log(f'Wrote {a.out} ({n / 1e6:.1f} MB, {len(V)} variants) in {(time.time() - t0) / 60:.1f} min')
    log('Put it next to brca1-variant-explorer.html and reload the page; the header then shows the release date.')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit('Interrupted: rerun to resume from the cache.')
