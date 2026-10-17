"""Synthetic missingness geometry cases, independent of clinical data."""
import io
import unittest
import zipfile

from research.oncology.improve_split_audit import BASE
from research.oncology.improve_missingness_geometry_audit import scan_geometry, _bin


class GeometryTests(unittest.TestCase):
    def test_fraction_bins(self):
        self.assertEqual(_bin(0), "0")
        self.assertEqual(_bin(0.04), "0.05")
        self.assertEqual(_bin(1), "1.0")

    def test_missing_values_not_zeros(self):
        b = io.BytesIO()
        with zipfile.ZipFile(b, "w") as z:
            z.writestr(BASE + "x_data/cancer_DNA_methylation.tsv",
                       "\tA\tB\nS1\t0\tNA\nS2\t1\t2\nOTHER\t3\t4\n")
        with zipfile.ZipFile(io.BytesIO(b.getvalue())) as z:
            a = scan_geometry(z, "cancer_DNA_methylation.tsv",
                              {"S1","S2"}, {"CCLE":{"S1"}, "gCSI":{"S2"}})
        self.assertEqual(a["response_linked_entity_count"], 2)
        self.assertEqual(a["row_missing_fraction_histogram_upper_bounds"]["0.5"], 1)
        self.assertEqual(a["row_missing_fraction_histogram_upper_bounds"]["0"], 1)
        self.assertEqual(a["column_missing_fraction_histogram_upper_bounds"]["0.5"], 1)
        self.assertEqual(a["per_study_unique_entity_descriptives"]["CCLE"]["missing_cells"], 1)
        self.assertEqual(a["per_study_unique_entity_descriptives"]["gCSI"]["missing_cells"], 0)
        self.assertEqual(a["global_valid_min"], 0)
        self.assertEqual(a["global_valid_max"], 2)

    def test_duplicate_linked_key_is_rejected(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr(BASE+"x_data/cancer_gene_expression.tsv",
                       "\tF\nS1\t2\nS1\t3\n")
        with zipfile.ZipFile(io.BytesIO(b.getvalue())) as z:
            with self.assertRaisesRegex(ValueError,"duplicate"):
                scan_geometry(z,"cancer_gene_expression.tsv",{"S1"},{"CCLE":{"S1"}})


if __name__ == "__main__":
    unittest.main()
