"""Change list between releases: python3 -m unittest discover tests  (standard library only)."""
import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import build_brca1_db as B


def digest(built, v, mave=None, sources=None):
    return {'built': built, 'sources': sources or {'clinvar': 'x'}, 'has': {'scv': True, 'lit': True}, 'v': v, 'mave': mave or {}}


class ChangesTest(unittest.TestCase):
    def test_classification_level_changes(self):
        prev = digest('2026-09-01T00:00:00Z', {
            'c.1A>G': {'cv': ['Uncertain significance', 'criteria provided, single submitter', '2024-01-01', 1], 'scv': [['SCV1.1', 'Lab A', 'Uncertain significance', '2024-01-01']]},
            'c.2A>G': {'gn': [3, 1000000, 5e-5, 'nfe']},
            'c.3A>G': {'lit': [[1], [2]]},
        }, {B.SGE_URN: {'mod': 'a', 's': {'c.4A>G': -0.5}}}, {'clinvar': 'old'})
        cur = digest('2026-10-01T00:00:00Z', {
            'c.1A>G': {'cv': ['Likely pathogenic', 'criteria provided, multiple submitters, no conflicts', '2026-08-01', 2],
                       'scv': [['SCV1.2', 'Lab A', 'Likely pathogenic', '2026-08-01'], ['SCV2.1', 'Lab B', 'Likely pathogenic', '2026-07-01']], 'en': ['Pathogenic', '2026-09-01']},
            'c.2A>G': {'gn': [30, 1000000, 2e-4, 'nfe']},
            'c.3A>G': {'lit': [[1, 5], [2]]},
        }, {B.SGE_URN: {'mod': 'b', 's': {'c.4A>G': -1.5}}}, {'clinvar': 'new'})
        ch = B.release_changes(prev, cur, {'c.1A>G', 'c.2A>G', 'c.3A>G', 'c.4A>G'})
        kinds = {c: {(i[0], i[1]) for i in items} for c, items in ch['v'].items()}
        self.assertIn(('clinvar', 1), kinds['c.1A>G'])
        self.assertIn(('review', 1), kinds['c.1A>G'])
        self.assertIn(('scv_new', 1), kinds['c.1A>G'])      # a kind of classification no submitter gave before
        self.assertIn(('scv_changed', 1), kinds['c.1A>G'])  # same submission, new class
        self.assertIn(('enigma', 1), kinds['c.1A>G'])
        self.assertIn(('faf', 1), kinds['c.2A>G'])          # crossed 1e-4
        self.assertIn(('gnomad_ac', 0), kinds['c.2A>G'])
        self.assertEqual(kinds['c.3A>G'], {('lit', 0)})
        self.assertIn(('sge', 1), kinds['c.4A>G'])          # functional -> non-functional
        self.assertEqual(ch['sources'], [['ClinVar', 'old', 'new']])
        self.assertEqual(ch['n'], {'variants': 4, 'major': 3})

    def test_faf_rule_matches_page(self):
        self.assertEqual(B.faf_of({'joint': {'faf95': {'popmax': 1e-4, 'popmax_population': 'afr'}}, 'exome': {'faf95': {'popmax': 5e-4}}}), (1e-4, 'afr'))
        self.assertEqual(B.faf_of({'exome': {'faf95': {'popmax': 2e-5, 'popmax_population': 'nfe'}}, 'genome': {'faf95': {'popmax': 3e-5, 'popmax_population': 'eas'}}}), (3e-5, 'eas'))
        self.assertEqual(B.faf_of(None), (None, None))


if __name__ == '__main__':
    unittest.main()
