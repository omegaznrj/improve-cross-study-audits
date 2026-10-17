"""Synthetic RDKit molecular correspondence tests; no biological records."""
import unittest
from research.oncology.improve_atom_mapped_stereo_audit import mapped_stereo_relation
from research.oncology.improve_structure_equivalence_audit import _rdkit


class AtomMappedStereoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.Chem, _ = _rdkit()
        except RuntimeError as exc:
            raise unittest.SkipTest(str(exc))

    def test_equivalent_atom_order(self):
        self.assertEqual(mapped_stereo_relation("CCO", "OCC", self.Chem),
                         "same_canonical_isomeric_structure")

    def test_unassigned_stereo(self):
        self.assertEqual(mapped_stereo_relation("F[C@H](Cl)Br", "FC(Cl)Br", self.Chem),
                         "incomplete_stereo_annotation")

    def test_conflicting_chiral_assignment(self):
        self.assertEqual(mapped_stereo_relation("F[C@H](Cl)Br", "F[C@@H](Cl)Br", self.Chem),
                         "assigned_stereo_conflict")

    def test_connectivity_distinct(self):
        self.assertEqual(mapped_stereo_relation("CCN", "CCO", self.Chem),
                         "different_connectivity")

    def test_invalid(self):
        self.assertEqual(mapped_stereo_relation("C(", "CCO", self.Chem),
                         "invalid_structure")


if __name__ == "__main__":
    unittest.main()
