"""Synthetic RDKit equivalence classifications; no source data."""
import unittest
from research.oncology.improve_structure_equivalence_audit import classify, _rdkit


class StructureEquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.Chem, _ = _rdkit()
        except RuntimeError as exc:
            raise unittest.SkipTest(str(exc))

    def test_ordering_equivalence(self):
        self.assertEqual(classify("CCO", "OCC", self.Chem),
                         "same_canonical_isomeric_structure")

    def test_stereochemical_difference(self):
        self.assertEqual(
            classify("F[C@H](Cl)Br", "F[C@@H](Cl)Br", self.Chem),
            "same_connectivity_stereo_difference_or_unspecified"
        )

    def test_fragment_discrepancy(self):
        self.assertEqual(
            classify("CCO", "CCO.[Na+]", self.Chem),
            "different_fragment_counts_review_salt_or_mixture"
        )

    def test_different_structure(self):
        self.assertEqual(
            classify("CCO", "CCN", self.Chem),
            "different_connectivity_or_other_structure_difference"
        )

    def test_invalid(self):
        self.assertEqual(classify("CCO", "not_a_smiles", self.Chem),
                         "invalid_or_blank_structure")


if __name__ == "__main__":
    unittest.main()
