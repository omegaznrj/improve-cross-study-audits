"""IMPROVE Phase 1.5C: targeted copy-number missing-mask concordance and Mordred tails.

Read-only ZIP, source-linked identities only. Reports aggregate counts, never
chemical identifiers, sample identifiers, individual feature names or outcomes.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
from pathlib import Path
import zipfile
from collections import Counter

from research.oncology.improve_split_audit import BASE, _read_response
from research.oncology.improve_numeric_quality_audit import _value_kind

COPY = ("cancer_copy_number.tsv", "cancer_discretized_copy_number.tsv")


def _table(z, name, ids):
    path = BASE + "x_data/" + name
    with z.open(path) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t", strict=True)
        headings = next(reader)
        if len(headings) != len(set(headings)) or len(headings) < 2:
            raise ValueError(name + ": duplicate/missing headings")
        seen = set()
        for values in reader:
            if len(values) != len(headings):
                raise ValueError(name + ": malformed row")
            key = values[0].strip()
            if key not in ids:
                continue
            if key in seen:
                raise ValueError(name + ": duplicate response-linked key")
            seen.add(key)
            yield key, headings[1:], values[1:]


def mask_comparison(z, ids):
    masks = {}
    headers = None
    for file in COPY:
        current = {}
        for key, cols, vals in _table(z, file, ids):
            if headers is None:
                headers = tuple(cols)
            elif headers != tuple(cols):
                raise ValueError("copy-number feature columns differ")
            current[key] = frozenset(i for i, val in enumerate(vals)
                                     if _value_kind(val)[0] != "valid")
        masks[file] = current
    a, b = masks[COPY[0]], masks[COPY[1]]
    shared = set(a) & set(b)
    return {
        "shared_linked_entities": len(shared),
        "missing_in_continuous_table_only": len(set(b) - set(a)),
        "missing_in_discretized_table_only": len(set(a) - set(b)),
        "entities_with_equal_missing_masks": sum(a[k] == b[k] for k in shared),
        "entities_with_different_missing_masks": sum(a[k] != b[k] for k in shared),
        "total_continuous_missing_cells_shared": sum(len(a[k]) for k in shared),
        "total_discretized_missing_cells_shared": sum(len(b[k]) for k in shared),
        "missing_cells_in_either_but_not_both": sum(len(a[k] ^ b[k]) for k in shared),
    }


def mordred_extremes(z, ids):
    counts = Counter()
    per_feature = {}
    for _, cols, vals in _table(z, "drug_mordred.tsv", ids):
        for col, raw in zip(cols, vals):
            kind, value = _value_kind(raw)
            counts[kind] += 1
            if kind != "valid":
                continue
            state = per_feature.setdefault(col, [math.inf, -math.inf])
            state[0] = min(state[0], value)
            state[1] = max(state[1], value)
            if abs(value) >= 1e6:
                counts["cells_abs_ge_1e6"] += 1
            if abs(value) >= 1e9:
                counts["cells_abs_ge_1e9"] += 1
    return {
        "numeric_counts": dict(sorted(counts.items())),
        "features_with_abs_value_ge_1e6": sum(max(abs(lo),abs(hi)) >= 1e6 for lo,hi in per_feature.values()),
        "features_with_abs_value_ge_1e9": sum(max(abs(lo),abs(hi)) >= 1e9 for lo,hi in per_feature.values()),
        "features_with_negative_values": sum(lo < 0 for lo,_ in per_feature.values()),
        "interpretation": "Descriptor magnitudes have different units/scales; large magnitudes alone do not establish an error.",
    }


def audit(path):
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate ZIP member names")
        response = _read_response(z)
        return {
            "schema": "tumorlab-improve-phase15c-v1",
            "copy_number_missing_mask_comparison": mask_comparison(z, set(response["sample_code_by_id"])),
            "mordred_magnitude_diagnostics": mordred_extremes(z, set(response["drug_code_by_id"])),
            "response_outcomes_read": False, "training_performed": False,
            "caveats": [
                "Reports only source-linked aggregate diagnostics.",
                "Magnitude checks are not domain verification; source descriptor units still require documentation.",
                "No feature removal, normalization, imputation, or split-informed fitting performed.",
            ],
        }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=Path(__file__).resolve().parents[2]/"data/raw/improve_cross_study/csa_data.zip")
    args = p.parse_args(argv)
    try:
        print(json.dumps(audit(args.archive), indent=2))
    except (OSError, ValueError, zipfile.BadZipFile, csv.Error, UnicodeDecodeError) as ex:
        p.exit(2, "IMPROVE Phase 1.5C: " + str(ex) + "\n")


if __name__ == "__main__":
    main()
