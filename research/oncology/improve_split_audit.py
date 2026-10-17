"""Read-only split-index and feature-key audit for IMPROVE CSA Zenodo v2.

Locally scans only response identifiers/source values and text split manifests.
Does not interpret or retain response outcomes, ingest omics matrices, or train models.
The file remains unextracted in the local ZIP. This is a structural audit,
not authorization to use the data for scientific learning.
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter
import csv
import io
from pathlib import Path
import zipfile

BASE = "csa_data/raw_data/"
STUDIES = ("CCLE", "CTRPv2", "GDSCv1", "GDSCv2", "gCSI")
MAX_ROWS = 2_000_000
MAX_INDICES = 2_000_000
REQUIRED_RESPONSE_FIELDS = ("source", "improve_sample_id", "improve_chem_id", "study")


def _read_response(zf: zipfile.ZipFile) -> dict:
    filename = BASE + "y_data/response.tsv"
    sample_codes = array("I")
    drug_codes = array("I")
    sources = array("I")
    samples: dict[str, int] = {}
    drugs: dict[str, int] = {}
    source_values: dict[str, int] = {}
    source_counts: Counter[str] = Counter()
    study_counts: Counter[str] = Counter()
    missing_sample_or_drug = 0
    with zf.open(filename) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("response.tsv has missing or duplicate heading")
        for name in REQUIRED_RESPONSE_FIELDS:
            if name not in reader.fieldnames:
                raise ValueError("response.tsv missing key column: " + name)
        for row_index, row in enumerate(reader):
            if row_index >= MAX_ROWS:
                raise ValueError("response table exceeds row-safety cap")
            if None in row:
                raise ValueError("response row has excess TSV fields: " + str(row_index))
            sample = (row["improve_sample_id"] or "").strip()
            drug = (row["improve_chem_id"] or "").strip()
            source = (row["source"] or "").strip()
            study = (row["study"] or "").strip()
            if not sample or not drug:
                missing_sample_or_drug += 1
            sample_codes.append(samples.setdefault(sample, len(samples)))
            drug_codes.append(drugs.setdefault(drug, len(drugs)))
            sources.append(source_values.setdefault(source, len(source_values)))
            source_counts[source] += 1
            study_counts[study] += 1
    return {
        "rows": len(sample_codes), "samples": sample_codes,
        "drugs": drug_codes, "source_codes": sources,
        "source_ids": source_values, "source_counts": dict(sorted(source_counts.items())),
        "study_counts": dict(sorted(study_counts.items())),
        "unique_sample_ids": len(samples), "unique_drug_ids": len(drugs),
        "missing_join_keys": missing_sample_or_drug,
        # Internal only: not included in the public/CLI summary.
        "sample_code_by_id": samples, "drug_code_by_id": drugs,
    }


def _indices(zf: zipfile.ZipFile, name: str, row_count: int) -> set[int]:
    full = BASE + "splits/" + name
    result: set[int] = set()
    total = 0
    with zf.open(full) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as stream:
        for line_number, line in enumerate(stream, 1):
            value = line.strip()
            if not value:
                continue
            total += 1
            if total > MAX_INDICES:
                raise ValueError(f"{name}: too many row indices")
            if not value.isascii() or not value.isdecimal():
                raise ValueError(f"{name}: invalid integer at line {line_number}")
            index = int(value)
            if not 0 <= index < row_count:
                raise ValueError(f"{name}: index outside response rows at line {line_number}")
            if index in result:
                raise ValueError(f"{name}: duplicate index at line {line_number}")
            result.add(index)
    return result


def _overlap_counts(a: set[int], b: set[int], rows: dict) -> dict:
    # A row split may share cell lines, drugs, or (cell, drug) pairs across groups.
    sa = {rows["samples"][i] for i in a}
    sb = {rows["samples"][i] for i in b}
    da = {rows["drugs"][i] for i in a}
    db = {rows["drugs"][i] for i in b}
    pa = {(rows["samples"][i], rows["drugs"][i]) for i in a}
    pb = {(rows["samples"][i], rows["drugs"][i]) for i in b}
    return {
        "overlapping_rows": len(a & b),
        "overlapping_sample_ids": len(sa & sb),
        "overlapping_drug_ids": len(da & db),
        "overlapping_sample_drug_pairs": len(pa & pb),
    }


def audit_split_indices(path: Path, fold: int = 0) -> dict:
    if type(fold) is not int or fold not in range(10):
        raise ValueError("fold must be between 0 and 9")
    if not path.is_file():
        raise ValueError("ZIP archive not found: " + str(path))
    with zipfile.ZipFile(path) as archive:
        members = archive.namelist()
        if len(members) != len(set(members)):
            raise ValueError("duplicate ZIP member paths")
        names = set(members)
        required = [BASE + "y_data/response.tsv"]
        for study in STUDIES:
            required += [BASE + f"splits/{study}_all.txt"]
            required += [BASE + f"splits/{study}_split_{fold}_{kind}.txt"
                         for kind in ("train", "val", "test")]
        missing = [name for name in required if name not in names]
        if missing:
            raise ValueError("missing response/split file: " + missing[0])
        rows = _read_response(archive)
        audit = {}
        for study in STUDIES:
            all_rows = _indices(archive, f"{study}_all.txt", rows["rows"])
            train = _indices(archive, f"{study}_split_{fold}_train.txt", rows["rows"])
            val = _indices(archive, f"{study}_split_{fold}_val.txt", rows["rows"])
            test = _indices(archive, f"{study}_split_{fold}_test.txt", rows["rows"])
            total_parts = train | val | test
            outside_study = total_parts - all_rows
            source_id = rows["source_ids"].get(study)
            mismatching_source_rows = (sum(rows["source_codes"][i] != source_id
                                           for i in all_rows) if source_id is not None else None)
            audit[study] = {
                "all_rows": len(all_rows),
                "train_rows": len(train), "validation_rows": len(val),
                "test_rows": len(test),
                "row_overlap_train_val": len(train & val),
                "row_overlap_train_test": len(train & test),
                "row_overlap_val_test": len(val & test),
                "partition_indices_outside_study_all": len(outside_study),
                "all_rows_not_covered_by_partitions": len(all_rows - total_parts),
                "study_all_indices_with_different_source_label": mismatching_source_rows,
                "train_test_entity_overlap": _overlap_counts(train, test, rows),
                "interpretation": "Row-disjoint is not sample/drug-disjoint; verify study protocol.",
            }
    return {
        "scope": "local source IDs and one official split fold only; no response metrics retained",
        "fold": fold,
        "response_row_count": rows["rows"],
        "response_source_counts": rows["source_counts"],
        "response_study_counts": rows["study_counts"],
        "unique_sample_ids": rows["unique_sample_ids"],
        "unique_drug_ids": rows["unique_drug_ids"],
        "rows_missing_sample_or_drug_key": rows["missing_join_keys"],
        "split_audit": audit,
        "notes": [
            "All split values treated as zero-based row positions in full response.tsv.",
            "No response labels, drug names, molecular values or entity IDs displayed.",
            "Entity overlap is a risk diagnostic, not necessarily a violation of the published protocol.",
            "Archive reuse/license terms remain unverified; no model training permitted by this audit.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path(__file__).resolve().parents[2] /
                        "data/raw/improve_cross_study/csa_data.zip")
    parser.add_argument("--fold", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        report = audit_split_indices(args.archive, args.fold)
    except (OSError, ValueError, zipfile.BadZipFile, UnicodeDecodeError, csv.Error) as exc:
        parser.exit(2, f"IMPROVE split audit: {exc}\n")
    print("TUMORLAB / READ-ONLY IMPROVE RESPONSE & SPLIT AUDIT")
    print("Response rows:", report["response_row_count"])
    print("Source counts:", report["response_source_counts"])
    print("Study field counts:", report["response_study_counts"])
    print("Unique sample/drug keys:", report["unique_sample_ids"], report["unique_drug_ids"])
    print("Response rows missing join keys:", report["rows_missing_sample_or_drug_key"])
    print("Official fold:", report["fold"])
    for study, result in report["split_audit"].items():
        print(f"  {study}: all={result['all_rows']}, train={result['train_rows']}, "
              f"val={result['validation_rows']}, test={result['test_rows']}")
        print("    Partition row overlaps:",
              result["row_overlap_train_val"], result["row_overlap_train_test"],
              result["row_overlap_val_test"])
        print("    Out of study / not covered:",
              result["partition_indices_outside_study_all"],
              result["all_rows_not_covered_by_partitions"])
        print("    Study/source mismatch:", result["study_all_indices_with_different_source_label"])
        print("    Train-test biological key overlaps:", result["train_test_entity_overlap"])
    print("All data remained inside local ZIP; response values not retained; .pt objects NOT loaded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
