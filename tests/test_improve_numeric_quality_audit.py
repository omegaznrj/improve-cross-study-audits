"""Entirely synthetic numerical feature QC tests; no original study data."""
import io
import unittest
import zipfile

from research.oncology.improve_split_audit import BASE
from research.oncology.improve_numeric_quality_audit import _scan_table, _value_kind


class NumericQualityTests(unittest.TestCase):
    def test_cell_categories(self):
        self.assertEqual(_value_kind("NaN")[0], "missing")
        self.assertEqual(_value_kind("inf")[0], "nonfinite")
        self.assertEqual(_value_kind("word")[0], "nonnumeric")
        self.assertEqual(_value_kind("1")[0], "valid")

    def test_binary_domain(self):
        self.assertEqual(_value_kind("2", binary=True)[0], "domain_invalid")
        self.assertEqual(_value_kind("1.0", binary=True)[0], "valid")

    def test_count_domain(self):
        self.assertEqual(_value_kind("-1", count=True)[0], "domain_invalid")
        self.assertEqual(_value_kind("1.5", count=True)[0], "domain_invalid")
        self.assertEqual(_value_kind("0", count=True)[0], "valid")

    def test_mock_molecular_table(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr(BASE + "x_data/cancer_gene_expression.tsv",
                       "\tgeneA\tgeneB\nCell_A\t1\t2\nCell_B\t1\tNA\nOther\t8\t9\n")
        with zipfile.ZipFile(io.BytesIO(buf.getvalue())) as z:
            stats, matched, complete = _scan_table(
                z, "cancer_gene_expression.tsv", {"Cell_A", "Cell_B"})
        self.assertEqual(stats["valid_cells"], 3)
        self.assertEqual(stats["missing_cells"], 1)
        self.assertEqual(stats["response_linked_ids_matched"], 2)
        self.assertEqual(stats["columns_constant_among_valid_observations"], 1)
        self.assertEqual(stats["columns_with_only_one_valid_observation"], 1)
        self.assertEqual(complete, {"Cell_A"})
        self.assertEqual(matched, {"Cell_A", "Cell_B"})

    def test_missing_table(self):
        with zipfile.ZipFile(io.BytesIO(), "w") as z:
            status, matched, complete = _scan_table(z, "drug_mordred.tsv", {"Drug_A"})
        self.assertEqual(status["status"], "missing_file")
        self.assertFalse(matched)
        self.assertFalse(complete)


if __name__ == "__main__":
    unittest.main()
