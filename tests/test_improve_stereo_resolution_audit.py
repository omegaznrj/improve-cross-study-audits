"""Synthetic, structure-only stereo triage tests."""
import unittest
from research.oncology.improve_stereo_resolution_audit import stereo_triage
from research.oncology.improve_structure_equivalence_audit import _rdkit


class StereoTriageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.Chem, _ = _rdkit()
        except RuntimeError as exc:
            raise unittest.SkipTest(str(exc))

    def test_same_connectivity_and_stereo(self):
        self.assertEqual(stereo_triage("CCO", "OCC", self.Chem),
                         "same_isomeric_structure")

    def test_unspecified_stereo(self):
        self.assertEqual(stereo_triage("F[C@H](Cl)Br", "FC(Cl)Br", self.Chem),
                         "stereo_annotation_incomplete_or_different")

    def test_specified_disagreement(self):
        self.assertEqual(stereo_triage("F[C@H](Cl)Br", "F[C@@H](Cl)Br", self.Chem),
                         "specified_stereo_disagreement_needs_atom_mapping")

    def test_connectivity_difference(self):
        self.assertEqual(stereo_triage("CCO", "CCN", self.Chem),
                         "other_structural_discrepancy")


if __name__ == "__main__":
    unittest.main()
