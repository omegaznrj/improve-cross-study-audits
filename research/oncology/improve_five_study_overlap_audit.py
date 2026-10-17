"""All-five-study IMPROVE source overlap and pooled-split risk audit (NO FIT).

Read only source identifier fields from original response.tsv, original zero-based
index lists and no outcome/omics columns. Report aggregate biological key counts,
not the keys themselves. Shared harmonized identifiers are a diagnostic rather
than proof of identical biological materials or an official protocol violation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import zipfile

from research.oncology.improve_archive_audit import DEFAULT_ARCHIVE
from research.oncology.improve_split_audit import (
    BASE, STUDIES, _read_response, _indices,
)

PARTS = ("train", "val", "test")


def _symbols(indices: set[int], response: dict) -> dict:
    samples, drugs = response["samples"], response["drugs"]
    missing_sample = response["sample_code_by_id"].get("")
    missing_drug = response["drug_code_by_id"].get("")
    valid = [
        (samples[i], drugs[i]) for i in indices
        if samples[i] != missing_sample and drugs[i] != missing_drug
    ]
    pair_freq = Counter(valid)
    return {
        "rows": len(indices),
        "keyed_rows": len(valid),
        "missing_key_rows": len(indices) - len(valid),
        "sample_keys": {s for s, _ in valid},
        "drug_keys": {d for _, d in valid},
        "pair_keys": set(pair_freq),
        "repeated_pair_excess_rows": sum(n - 1 for n in pair_freq.values()),
        "pairs_with_multiple_rows": sum(n > 1 for n in pair_freq.values()),
    }


def _intersection(a: dict, b: dict) -> dict:
    return {
        "sample_ids": len(a["sample_keys"] & b["sample_keys"]),
        "drug_ids": len(a["drug_keys"] & b["drug_keys"]),
        "sample_drug_pairs": len(a["pair_keys"] & b["pair_keys"]),
    }


def audit_five_study_overlap(path: Path = DEFAULT_ARCHIVE, fold: int = 0) -> dict:
    if type(fold) is not int or not 0 <= fold <= 9:
        raise ValueError("fold must be an integer from 0 through 9")
    if not path.is_file():
        raise ValueError("archive not found")
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicated ZIP paths")
        needed = [BASE + "y_data/response.tsv"]
        for study in STUDIES:
            needed.extend([BASE + f"splits/{study}_all.txt"])
            needed.extend(BASE + f"splits/{study}_split_{fold}_{part}.txt"
                          for part in PARTS)
        for name in needed:
            if name not in names:
                raise ValueError("missing required source manifest")
        response = _read_response(z)
        membership = {}
        views = {}
        by_study = {}
        for study in STUDIES:
            partitions = {"all": _indices(z, f"{study}_all.txt", response["rows"])}
            partitions.update({
                p: _indices(z, f"{study}_split_{fold}_{p}.txt", response["rows"])
                for p in PARTS
            })
            membership[study] = partitions
            v = {k: _symbols(indices, response)
                 for k, indices in partitions.items()}
            views[study] = v
            source_code = response["source_ids"].get(study)
            by_study[study] = {
                "all_rows": len(partitions["all"]),
                "unique_cell_models": len(v["all"]["sample_keys"]),
                "unique_drugs": len(v["all"]["drug_keys"]),
                "unique_model_drug_pairs": len(v["all"]["pair_keys"]),
                "repeated_assay_pair_excess_rows": v["all"]["repeated_pair_excess_rows"],
                "pairs_with_multiple_rows": v["all"]["pairs_with_multiple_rows"],
                "missing_biological_keys": v["all"]["missing_key_rows"],
                "source_label_mismatch": (
                    None if source_code is None else
                    sum(response["source_codes"][i] != source_code
                        for i in partitions["all"])
                ),
                "fold_partition_rows": {p: len(partitions[p]) for p in PARTS},
                "split_row_index_overlap": {
                    "train_val": len(partitions["train"] & partitions["val"]),
                    "train_test": len(partitions["train"] & partitions["test"]),
                    "val_test": len(partitions["val"] & partitions["test"]),
                },
                "partition_indices_not_in_study_all": len(
                    set.union(*(partitions[p] for p in PARTS)) - partitions["all"]
                ),
                "study_all_rows_not_covered": len(
                    partitions["all"] - set.union(*(partitions[p] for p in PARTS))
                ),
                "within_source_pair_intersections": {
                    "train_val": _intersection(v["train"], v["val"]),
                    "train_test": _intersection(v["train"], v["test"]),
                    "val_test": _intersection(v["val"], v["test"]),
                },
            }
        comparisons = {}
        for left_pos, left in enumerate(STUDIES):
            for right in STUDIES[left_pos+1:]:
                lv, rv = views[left], views[right]
                pair_name = left + "_vs_" + right
                comparisons[pair_name] = {
                    "all_vs_all": _intersection(lv["all"], rv["all"]),
                    "train_vs_foreign_validation": _intersection(lv["train"], rv["val"]),
                    "train_vs_foreign_test": _intersection(lv["train"], rv["test"]),
                    "foreign_train_vs_validation": _intersection(rv["train"], lv["val"]),
                    "foreign_train_vs_test": _intersection(rv["train"], lv["test"]),
                    "shared_raw_response_row_indices": len(
                        membership[left]["all"] & membership[right]["all"]),
                }

        all_train_pairs = set.union(*(views[s]["train"]["pair_keys"] for s in STUDIES))
        all_valid_pairs = set.union(*(views[s]["val"]["pair_keys"] for s in STUDIES))
        all_test_pairs = set.union(*(views[s]["test"]["pair_keys"] for s in STUDIES))
        cross_study_train_to_heldout = set()
        for a in STUDIES:
            for b in STUDIES:
                if a != b:
                    cross_study_train_to_heldout |= (
                        views[a]["train"]["pair_keys"] &
                        (views[b]["val"]["pair_keys"] |
                         views[b]["test"]["pair_keys"])
                    )
        return {
            "schema": "tumorlab-improve-five-study-overlap-v1",
            "fold": fold,
            "response_total_rows": response["rows"],
            "studies": by_study,
            "study_pair_overlaps": comparisons,
            "study_pair_count": len(comparisons),
            "pooled_partition_risks": {
                "pooled_train_unique_pairs": len(all_train_pairs),
                "pooled_validation_unique_pairs": len(all_valid_pairs),
                "pooled_test_unique_pairs": len(all_test_pairs),
                "train_vs_validation_shared_pairs": len(all_train_pairs & all_valid_pairs),
                "train_vs_test_shared_pairs": len(all_train_pairs & all_test_pairs),
                "cross_source_train_to_other_source_heldout_pairs": len(
                    cross_study_train_to_heldout
                ),
            },
            "no_response_outcome_values_used": True,
            "no_raw_cell_or_drug_identifiers_reported": True,
            "no_training_performed": True,
            "interpretation": [
                "Harmonized cell/drug identifiers are not a confirmed physical sample or structure crosswalk.",
                "Repeated within-source model-drug pairs may have legitimate different assays or replicates.",
                "Original fold row-index validity and grouping are separate from pooled-source leakage.",
                "Cross-source train-to-heldout intersections prevent a naive independent pooled-pair claim.",
                "No individual response targets, molecular values or processed .pt files are used.",
            ],
        }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--fold", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        result = audit_five_study_overlap(args.archive, args.fold)
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        parser.exit(2, "IMPROVE five-study overlap audit: " + str(exc) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
