"""Phase 1.5B: aggregate missingness geometry and numerical range diagnostics.

Streams source TSVs in original ZIP. Reports feature/sample histograms only.
Study membership and train partition are descriptive, not independent samples
or leakage-safe fitted preprocessing statistics.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import io
import json
import math
from pathlib import Path
import zipfile

from research.oncology.improve_split_audit import BASE, STUDIES, _read_response, _indices
from research.oncology.improve_numeric_quality_audit import (
    TABLES, DRUG_NUMERIC, _value_kind, MAX_FEATURES, MAX_TABLE_ROWS,
)

BINS = (0, 0.001, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0)


def _bin(fraction):
    for upper in BINS:
        if fraction <= upper + 1e-12:
            return str(upper)
    raise ValueError("fraction outside [0,1]")


def scan_geometry(z, table, keys, study_entities):
    member = BASE + "x_data/" + table
    if member not in z.namelist():
        return {"status": "missing_file"}
    binary = table == "drug_ecfp4_nbits512.tsv"
    count = table == "cancer_mutation_count.tsv"
    col_bad, col_valid, col_zero = [], [], []
    row_missing = Counter()
    per_study = {study: {"entity_ids_in_matrix": 0, "missing_cells": 0,
                          "valid_cells": 0, "entities_with_any_missing": 0}
                 for study in STUDIES}
    seen = set()
    minval = math.inf
    maxval = -math.inf
    with z.open(member) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as fh:
        header = next(csv.reader([fh.readline()], delimiter="\t"))
        n = len(header) - 1
        if n <= 0 or n > MAX_FEATURES or len(header) != len(set(header)):
            raise ValueError(table + ": invalid headings")
        col_bad, col_valid, col_zero = [0]*n, [0]*n, [0]*n
        for index, row in enumerate(csv.reader(fh, delimiter="\t", strict=True)):
            if index >= MAX_TABLE_ROWS or len(row) != n + 1:
                raise ValueError(table + ": invalid/oversized rows")
            key = row[0].strip()
            if key not in keys:
                continue
            if key in seen:
                raise ValueError(table + ": duplicate response-linked feature key")
            seen.add(key)
            bad = 0
            zeros = 0
            for j, rawval in enumerate(row[1:]):
                kind, value = _value_kind(rawval, binary=binary, count=count)
                if kind == "valid":
                    col_valid[j] += 1
                    if value == 0:
                        col_zero[j] += 1
                        zeros += 1
                    minval = min(minval, value)
                    maxval = max(maxval, value)
                else:
                    col_bad[j] += 1
                    bad += 1
            row_missing[_bin(bad/n)] += 1
            for study, entity_ids in study_entities.items():
                if key in entity_ids:
                    d = per_study[study]
                    d["entity_ids_in_matrix"] += 1
                    d["missing_cells"] += bad
                    d["valid_cells"] += n-bad
                    d["entities_with_any_missing"] += bad > 0
    nrows = len(seen)
    column_missing_bins = Counter(_bin(b/nrows) for b in col_bad) if nrows else Counter()
    zero_bins = Counter(_bin(z/v) for z,v in zip(col_zero,col_valid) if v)
    return {
        "status": "scanned", "response_linked_entity_count": nrows,
        "feature_columns": n,
        "row_missing_fraction_histogram_upper_bounds": dict(sorted(row_missing.items(),key=lambda kv:float(kv[0]))),
        "column_missing_fraction_histogram_upper_bounds": dict(sorted(column_missing_bins.items(),key=lambda kv:float(kv[0]))),
        "column_zero_fraction_among_valid_histogram_upper_bounds": dict(sorted(zero_bins.items(),key=lambda kv:float(kv[0]))),
        "columns_no_valid": sum(v==0 for v in col_valid),
        "columns_missing_more_than_half": sum(b/nrows>0.5 for b in col_bad) if nrows else 0,
        "columns_zero_more_than_95pct_of_valid": sum(v>0 and z/v>0.95 for z,v in zip(col_zero,col_valid)),
        "global_valid_min": minval if nrows and minval != math.inf else None,
        "global_valid_max": maxval if nrows and maxval != -math.inf else None,
        "per_study_unique_entity_descriptives": per_study,
    }


def audit_geometry(path, fold=0):
    if type(fold) is not int or fold not in range(10):
        raise ValueError("fold must be 0..9")
    with zipfile.ZipFile(path) as z:
        if len(z.namelist()) != len(set(z.namelist())):
            raise ValueError("duplicate ZIP member")
        rows = _read_response(z)
        split = {s: _indices(z, f"{s}_split_{fold}_train.txt", rows["rows"])
                 | _indices(z, f"{s}_split_{fold}_val.txt", rows["rows"])
                 | _indices(z, f"{s}_split_{fold}_test.txt", rows["rows"])
                 for s in STUDIES}
        tables = {}
        for table in TABLES:
            isdrug = table in DRUG_NUMERIC
            codes = rows["drugs"] if isdrug else rows["samples"]
            id_map = rows["drug_code_by_id"] if isdrug else rows["sample_code_by_id"]
            bystudy = {study: {key for key,c in id_map.items()
                               if c in {codes[i] for i in indices}}
                       for study, indices in split.items()}
            tables[table] = scan_geometry(z, table, set(id_map), bystudy)
    return {
        "schema": "tumorlab-improve-phase15b-geometry-v1",
        "fold": fold, "feature_tables": tables,
        "outcomes_read": False, "training_performed": False,
        "notes": [
            "Histogram labels are inclusive upper bounds on a fraction, not quantiles.",
            "Study entity groups overlap; counts must not be added across studies.",
            "Global min/max are diagnostic, not source-supported valid biological bounds.",
            "Methylation missing values are not zero observations; zeros are counted separately.",
            "Study summaries cover unique linked entities, not repeated response rows.",
            "All diagnostics descriptive and pooled; no test-informed imputation or column selection.",
            "Raw IDs, per-entity masks and molecular values are never exported.",
        ],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=Path(__file__).resolve().parents[2] / "data/raw/improve_cross_study/csa_data.zip")
    p.add_argument("--fold", type=int, default=0)
    a = p.parse_args(argv)
    try:
        print(json.dumps(audit_geometry(a.archive, a.fold), indent=2))
    except (OSError, ValueError, zipfile.BadZipFile, UnicodeDecodeError, csv.Error) as exc:
        p.exit(2, "IMPROVE geometry QC: " + str(exc) + "\n")


if __name__ == "__main__":
    main()
