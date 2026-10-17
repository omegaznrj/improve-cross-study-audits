"""Phase 1.5: streaming, read-only numerical QC of IMPROVE molecular matrices.

No response target values, fitting, row-level exports or changes to original ZIP.
Cell-feature keys must match response sample IDs exactly; drug feature keys must
match response chemical IDs exactly. Reports table-level and per-study/partition
row *coverage* (not replicated pseudo-measurements).
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import io
import json
import math
from pathlib import Path
import zipfile

from research.oncology.improve_split_audit import BASE, STUDIES, _indices, _read_response
from research.oncology.improve_feature_key_audit import SAMPLE_FEATURES

DRUG_NUMERIC = ("drug_ecfp4_nbits512.tsv", "drug_mordred.tsv")
TABLES = SAMPLE_FEATURES + DRUG_NUMERIC
PARTITIONS = ("train", "val", "test")
MAX_FEATURES = 100000
MAX_TABLE_ROWS = 2000000
MAX_LINE_BYTES = 32 * 1024 * 1024
MISSING = {"", "na", "nan", "null", "none", "n/a"}


def _value_kind(value: str, *, binary: bool = False, count: bool = False):
    s = value.strip()
    if s.lower() in MISSING:
        return "missing", None
    try:
        num = float(s)
    except ValueError:
        return "nonnumeric", None
    if not math.isfinite(num):
        return "nonfinite", None
    if binary and num not in (0.0, 1.0):
        return "domain_invalid", num
    if count and (num < 0 or not num.is_integer()):
        return "domain_invalid", num
    return "valid", num


def _scan_table(z: zipfile.ZipFile, table: str, expected: set[str]) -> tuple[dict, set[str], set[str]]:
    path = BASE + "x_data/" + table
    if path not in z.namelist():
        return {"status": "missing_file"}, set(), set()
    is_binary = table == "drug_ecfp4_nbits512.tsv"
    is_count = table == "cancer_mutation_count.tsv"
    seen, complete = set(), set()
    counts = Counter()
    # Column-level constancy only among valid measurements; report counts,
    # not features or values. Avoid full matrix materialization.
    first = []
    changed = []
    valid_per_col = []
    missing_per_col = []
    with z.open(path) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as fh:
        header_raw = fh.readline(MAX_LINE_BYTES + 1)
        if not header_raw or len(header_raw.encode("utf-8")) > MAX_LINE_BYTES:
            raise ValueError(table + ": absent/oversized header")
        header = next(csv.reader([header_raw], delimiter="\t", strict=True))
        if len(header) < 2 or len(header) > MAX_FEATURES + 1 or len(header) != len(set(header)):
            raise ValueError(table + ": malformed/duplicate feature headings")
        if table in DRUG_NUMERIC and header[0] != "improve_chem_id":
            raise ValueError(table + ": wrong drug identifier heading")
        if is_binary and header[1:] != [f"ecfp4.{i}" for i in range(512)]:
            raise ValueError("ECFP4 columns must be 512 ordered bits")
        n = len(header) - 1
        first = [None] * n
        changed = [False] * n
        valid_per_col = [0] * n
        missing_per_col = [0] * n
        reader = csv.reader(fh, delimiter="\t", strict=True)
        for row_index, cells in enumerate(reader):
            if row_index >= MAX_TABLE_ROWS:
                raise ValueError(table + ": safety row cap")
            if len(cells) != n + 1:
                raise ValueError(table + ": inconsistent field count at row " + str(row_index))
            counts["data_rows"] += 1
            key = cells[0].strip()
            if not key:
                counts["blank_id_rows"] += 1
            elif key in seen:
                counts["duplicate_id_rows"] += 1
            seen.add(key)
            if key not in expected:
                continue
            counts["matching_rows"] += 1
            row_good = True
            for i, raw_value in enumerate(cells[1:]):
                kind, num = _value_kind(raw_value, binary=is_binary, count=is_count)
                counts[kind + "_cells"] += 1
                if kind != "valid":
                    row_good = False
                    missing_per_col[i] += 1
                else:
                    valid_per_col[i] += 1
                    if first[i] is None:
                        first[i] = num
                    elif num != first[i]:
                        changed[i] = True
            if row_good:
                complete.add(key)
    matched = seen & expected
    return {
        "status": "scanned",
        "feature_columns": n,
        "rows_in_table": counts["data_rows"],
        "response_linked_ids_matched": len(matched),
        "blank_key_rows": counts["blank_id_rows"],
        "duplicate_key_rows": counts["duplicate_id_rows"],
        "matched_data_rows": counts["matching_rows"],
        "valid_cells": counts["valid_cells"],
        "missing_cells": counts["missing_cells"],
        "nonfinite_cells": counts["nonfinite_cells"],
        "nonnumeric_cells": counts["nonnumeric_cells"],
        "domain_invalid_cells": counts["domain_invalid_cells"],
        "columns_with_no_valid_observations": sum(v == 0 for v in valid_per_col),
        "columns_with_only_one_valid_observation": sum(v == 1 for v in valid_per_col),
        "columns_constant_among_valid_observations": sum(v >= 2 and not ch for v, ch in zip(valid_per_col, changed)),
        "columns_with_any_missing_invalid": sum(m > 0 for m in missing_per_col),
        "response_linked_ids_with_complete_numeric_rows": len(complete & expected),
        "domain_check": "binary_bits" if is_binary else ("nonnegative_integer" if is_count else "finite_numeric_only"),
    }, matched, complete & expected


def audit_numeric_quality(path: Path, fold: int = 0):
    if type(fold) is not int or fold not in range(10):
        raise ValueError("fold must be 0..9")
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate ZIP member names")
        response = _read_response(z)
        splits = {study: {
            part: _indices(z, f"{study}_split_{fold}_{part}.txt", response["rows"])
            for part in PARTITIONS
        } for study in STUDIES}
        output = {}
        for table in TABLES:
            entity = "drug" if table in DRUG_NUMERIC else "sample"
            id_map = response["drug_code_by_id"] if entity == "drug" else response["sample_code_by_id"]
            encoded = response["drugs"] if entity == "drug" else response["samples"]
            metric, present, complete = _scan_table(z, table, set(id_map))
            present_codes = {id_map[k] for k in present}
            complete_codes = {id_map[k] for k in complete}
            coverage = {}
            for study in STUDIES:
                coverage[study] = {}
                for part in PARTITIONS:
                    ix = splits[study][part]
                    coverage[study][part] = {
                        "response_rows": len(ix),
                        "rows_with_matching_entity": sum(encoded[i] in present_codes for i in ix),
                        "rows_with_complete_numeric_entity": sum(encoded[i] in complete_codes for i in ix),
                    }
            metric["per_study_partition"] = coverage
            output[table] = metric
    return {
        "schema": "tumorlab-improve-numeric-qc-phase15-v1",
        "fold": fold,
        "response_rows": response["rows"],
        "feature_tables": output,
        "outcomes_read": False,
        "model_training": False,
        "limitations": [
            "Matched entity-row counts reuse features across many response records; not independent measurements.",
            "Constant-feature results are pooled across response-linked entities, not train-only.",
            "Numerically valid does not imply experimentally comparable, normalized or biologically valid.",
            "No missingness pattern, quantiles, assay batch provenance, source preparation, or endpoint analysis.",
            "One-to-many feature keys must never be naively joined to response rows.",
        ],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=Path(__file__).resolve().parents[2] / "data/raw/improve_cross_study/csa_data.zip")
    p.add_argument("--fold", type=int, default=0)
    args = p.parse_args(argv)
    try:
        print(json.dumps(audit_numeric_quality(args.archive, args.fold), indent=2))
    except (OSError, ValueError, UnicodeDecodeError, csv.Error, zipfile.BadZipFile) as exc:
        p.exit(2, "IMPROVE numeric QC: " + str(exc) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
