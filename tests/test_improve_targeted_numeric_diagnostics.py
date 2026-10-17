"""Synthetic Phase 1.5C targeted checks."""
import io
import unittest
import zipfile

from research.oncology.improve_split_audit import BASE
from research.oncology.improve_targeted_numeric_diagnostics import mask_comparison, mordred_extremes


def fixture(tables):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, text in tables.items():
            z.writestr(BASE + "x_data/" + name, text)
    return zipfile.ZipFile(io.BytesIO(buf.getvalue()))


class TargetedQualityTests(unittest.TestCase):
    def test_copy_mask_exact(self):
        with fixture({
            "cancer_copy_number.tsv": "\tA\tB\nS1\tNA\t2\nS2\t1\t2\n",
            "cancer_discretized_copy_number.tsv": "\tA\tB\nS1\tNaN\t1\nS2\t0\t-1\n",
        }) as z:
            report = mask_comparison(z, {"S1", "S2"})
        self.assertEqual(report["entities_with_equal_missing_masks"], 2)
        self.assertEqual(report["missing_cells_in_either_but_not_both"], 0)

    def test_copy_mask_difference(self):
        with fixture({
            "cancer_copy_number.tsv": "\tA\nS1\tNA\n",
            "cancer_discretized_copy_number.tsv": "\tA\nS1\t2\n",
        }) as z:
            report = mask_comparison(z, {"S1"})
        self.assertEqual(report["missing_cells_in_either_but_not_both"], 1)

    def test_mordred_extremes(self):
        with fixture({
            "drug_mordred.tsv": "improve_chem_id\tA\tB\nD1\t2000000000\t-1\nD2\t0\t0\n",
        }) as z:
            report = mordred_extremes(z, {"D1", "D2"})
        self.assertEqual(report["features_with_abs_value_ge_1e9"], 1)
        self.assertEqual(report["features_with_negative_values"], 1)


if __name__ == "__main__":
    unittest.main()
