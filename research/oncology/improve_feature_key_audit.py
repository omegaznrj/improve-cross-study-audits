"""IMPROVE Zenodo v2 ZIP: read-only exact identifier/feature-availability audit.

Scans only the first identifier field of molecular feature table rows and the
explicit improve_chem_id column in drug_info. Skips all scientific feature
values; reports only aggregate join coverage and duplicate/blank keys.
No full matrix, .pt, target value, or patient-level data is loaded.
"""
from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path
import zipfile

from research.oncology.improve_split_audit import (
    BASE, STUDIES, _indices, _read_response,
)

MAX_LINE_BYTES = 16 * 1024 * 1024
MAX_KEY_BYTES = 2048

# Molecular tables exported from Pandas have one unnamed first index heading.
# This is a HYPOTHESIS about identifier semantics until actual values match.
SAMPLE_FEATURES = (
    "cancer_gene_expression.tsv",
    "cancer_copy_number.tsv",
    "cancer_discretized_copy_number.tsv",
    "cancer_DNA_methylation.tsv",
    "cancer_RPPA.tsv",
    "cancer_miRNA_expression.tsv",
    "cancer_mutation_count.tsv",
)
DRUG_FEATURES = (
    "drug_SMILES.tsv",
    "drug_ecfp4_nbits512.tsv",
    "drug_mordred.tsv",
    "drug_info.tsv",
)


def _key_table(archive: zipfile.ZipFile, table: str, entity: str,
               expected_ids: dict[str, int]) -> dict:
    name = BASE + "x_data/" + table
    if name not in archive.namelist():
        return {
            "status": "missing_file", "file": table, "entity": entity,
            "data_rows": 0, "unique_nonblank_ids": 0,
            "unique_matching_response_ids": 0, "duplicate_key_rows": 0,
            "empty_id_rows": 0, "matched_codes": set(),
            "header_index_name": None, "key_column_position": None,
        }
    rows = 0
    matches: set[int] = set()
    unique: set[str] = set()
    duplicate = 0
    empty = 0
    with archive.open(name) as stream:
        first = stream.readline(MAX_LINE_BYTES + 1)
        if not first or len(first) > MAX_LINE_BYTES:
            raise ValueError("missing/oversized header: " + table)
        header = first.decode("utf-8-sig").rstrip("\r\n").split("\t")
        column = 0
        if table == "drug_info.tsv":
            if len(header) != len(set(header)) or header.count("improve_chem_id") != 1:
                raise ValueError("drug_info missing unique improve_chem_id heading")
            column = header.index("improve_chem_id")
        else:
            if entity == "drug" and (not header or header[0] != "improve_chem_id"):
                raise ValueError("unexpected drug feature key header: " + table)
        # Only the selected key cell is interpreted. Every scientific feature
        # value after that remains uninterpreted and never enters the output.
        for row in stream:
            rows += 1
            if rows > 2_000_000:
                raise ValueError("feature row safety cap exceeded: " + table)
            if len(row) > MAX_LINE_BYTES:
                raise ValueError("feature line safety cap exceeded: " + table)
            if column == 0:
                raw = row.partition(b"\t")[0]
            else:
                # drug_info.tsv is 7 small columns; no wide omics parsed here.
                fields = next(csv.reader([row.decode("utf-8-sig").rstrip("\r\n")],
                                         delimiter="\t", strict=True))
                if len(fields) != len(header):
                    raise ValueError("drug_info inconsistent column count")
                raw = fields[column].encode("utf-8")
            if len(raw) > MAX_KEY_BYTES:
                raise ValueError("identifier too long: " + table)
            key = raw.decode("utf-8-sig").strip()
            if not key:
                empty += 1
                continue
            if key in unique:
                duplicate += 1
            unique.add(key)
            code = expected_ids.get(key)
            if code is not None:
                matches.add(code)
    return {
        "status": "keys_scanned", "file": table, "entity": entity,
        "data_rows": rows, "unique_nonblank_ids": len(unique),
        "unique_matching_response_ids": len(matches),
        "duplicate_key_rows": duplicate, "empty_id_rows": empty,
        "header_index_name": header[column] or None,
        "key_column_position": column,
        "matched_codes": matches,
    }


def _counts_per_study(rows: dict, splits: dict, matched_codes: set[int],
                      entity: str) -> dict:
    encoded = rows["samples"] if entity == "sample" else rows["drugs"]
    result = {}
    for study, partitions in splits.items():
        result[study] = {}
        for part, positions in partitions.items():
            matched = sum(encoded[i] in matched_codes for i in positions)
            result[study][part] = {
                "response_rows": len(positions), "exact_key_matched_rows": matched,
                "unmatched_rows": len(positions) - matched,
            }
    return result


def audit_feature_keys(path: Path, fold: int = 0) -> dict:
    """Inspect exact join-key coverage; NO join or matrix values are emitted."""
    if type(fold) is not int or fold not in range(10):
        raise ValueError("fold must be 0 to 9")
    if not path.is_file():
        raise ValueError("ZIP not found: " + str(path))
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate ZIP member paths")
        rows = _read_response(archive)
        splits = {
            study: {
                part: _indices(archive, f"{study}_split_{fold}_{part}.txt", rows["rows"])
                for part in ("train", "val", "test")
            } for study in STUDIES
        }
        output = {}
        for entity, nameset in (("sample", SAMPLE_FEATURES), ("drug", DRUG_FEATURES)):
            identity_map = rows["sample_code_by_id"] if entity == "sample" else rows["drug_code_by_id"]
            for table in nameset:
                item = _key_table(archive, table, entity, identity_map)
                codes = item.pop("matched_codes")
                item["per_study"] = _counts_per_study(rows, splits, codes, entity)
                item["join_status"] = (
                    "file_not_present" if item["status"] == "missing_file"
                    else ("no_exact_key_matches_needs_crosswalk" if not codes
                          else "exact_key_coverage_only_not_value_quality")
                )
                output[table] = item
    return {
        "source": "IMPROVE CSA Zenodo v2 DOI 10.5281/zenodo.15258883",
        "fold": fold,
        "scope": "exact key matches only; training/validation/test values not parsed",
        "response_rows": rows["rows"], "feature_tables": output,
        "deliberately_unread": ["x_data/cancer_mutation.parquet",
                               "processed/*.pt", "raw scientific feature values",
                               "response target values"],
        "limitations": [
            "Matching identifier strings does not establish biological specimen identity.",
            "Key presence does not establish nonmissing or valid molecular measurements.",
            "Unavailable or differently named IDs require an evidence-backed crosswalk.",
            "Duplicate keys may indicate repeated rows and preclude a one-to-one join.",
            "Drug response in cultured cell lines is not patient treatment outcome.",
            "Data-use rights and target freeze remain separate pre-training checks.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path(__file__).resolve().parents[2] /
                        "data/raw/improve_cross_study/csa_data.zip")
    parser.add_argument("--fold", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        report = audit_feature_keys(args.archive, args.fold)
    except (OSError, ValueError, zipfile.BadZipFile, UnicodeDecodeError, csv.Error) as exc:
        parser.exit(2, f"IMPROVE feature key audit: {exc}\n")
    print("TUMORLAB / IMPROVE FEATURE IDENTIFIER COVERAGE -- READ ONLY")
    print("Response rows:", report["response_rows"], "fold:", report["fold"])
    for name, item in report["feature_tables"].items():
        print(f"  {name}: {item['join_status']}; unique file keys="
              f"{item['unique_nonblank_ids']}; matched response IDs="
              f"{item['unique_matching_response_ids']}; duplicate key rows="
              f"{item['duplicate_key_rows']}; empty keys={item['empty_id_rows']}")
        for study in ("CCLE", "gCSI"):
            counts = item["per_study"][study]
            summary = ", ".join(
                f"{part} {counts[part]['exact_key_matched_rows']}/{counts[part]['response_rows']}"
                for part in ("train", "val", "test")
            )
            print(f"    {study} exact-key coverage: {summary}")
    print("NOTE: No scientific feature values or response labels loaded, no training.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
