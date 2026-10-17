"""Read-only IMPROVE Zenodo v2 ZIP inventory/header audit; no data extraction.

Does not deserialize .pt/.gz, import patient records, analyze response values,
change official splits, or authorize reuse for training.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import zipfile

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ARCHIVE = ROOT / "data" / "raw" / "improve_cross_study" / "csa_data.zip"
PUBLISHED_MD5 = "984eabec912914e62d56d8e0073475f8"
BASE = "csa_data/raw_data/"
STUDIES = ("CCLE", "CTRPv2", "GDSCv1", "GDSCv2", "gCSI")
REQUIRED_TABLES = (
    "y_data/response.tsv", "x_data/cancer_gene_expression.tsv",
    "x_data/drug_SMILES.tsv", "x_data/drug_info.tsv",
)
OTHER_OMICS = (
    "x_data/cancer_DNA_methylation.tsv", "x_data/cancer_RPPA.tsv",
    "x_data/cancer_copy_number.tsv", "x_data/cancer_discretized_copy_number.tsv",
    "x_data/cancer_miRNA_expression.tsv", "x_data/cancer_mutation.parquet",
    "x_data/cancer_mutation_count.tsv", "x_data/cancer_mutation_long_format.tsv",
    "x_data/drug_ecfp4_nbits512.tsv", "x_data/drug_mordred.tsv",
)
SPLIT_PATTERN = re.compile(
    r"^(CCLE|CTRPv2|GDSCv1|GDSCv2|gCSI)_(?:(all)|split_(\d+)_(train|val|test))\.txt$"
)
MAX_MEMBERS = 10000
MAX_TOTAL_DECLARED_BYTES = 25 * 1024 ** 3
MAX_HEADER_BYTES = 32 * 1024 * 1024  # bound memory for very wide tables


def verify_published_md5(path: Path) -> dict:
    """Streams local ZIP bytes; not a cryptographic publisher identity check."""
    digest = hashlib.md5()  # Published archive checksum, NOT used for security.
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    return {"published_md5": PUBLISHED_MD5, "actual_md5": actual,
            "matches_published_md5": actual == PUBLISHED_MD5}


def _check_zip_metadata(items: list[zipfile.ZipInfo]) -> None:
    if len(items) > MAX_MEMBERS:
        raise ValueError("too many ZIP entries")
    seen = set()
    total = 0
    for item in items:
        name = item.filename
        path = PurePosixPath(name)
        if (not name or "\\" in name or path.is_absolute() or ".." in path.parts
                or name.startswith("/") or ":" in path.parts[0]):
            raise ValueError("unsafe ZIP member path")
        if name in seen:
            raise ValueError("duplicate ZIP member")
        seen.add(name)
        if item.flag_bits & 1:
            raise ValueError("encrypted ZIP member is unsupported")
        total += item.file_size
        if total > MAX_TOTAL_DECLARED_BYTES:
            raise ValueError("ZIP declared uncompressed size exceeds audit limit")


def _header(archive: zipfile.ZipFile, name: str) -> dict:
    """Inspect only the first TSV physical line, capped at 32 MiB.

    Wide molecular tables can have 100,000+ columns in their FIRST line.
    Over-cap headers are reported rather than aborting every other table.
    A final header line without a newline is also valid.
    """
    result = {"table": name.removeprefix(BASE)}
    with archive.open(name) as stream:
        first = stream.readline(MAX_HEADER_BYTES + 1)
    if not first:
        raise ValueError("empty TSV header: " + name)
    if len(first) > MAX_HEADER_BYTES:
        return {**result, "status": "header_exceeds_safety_cap",
                "header_bytes_examined": len(first),
                "header_byte_limit": MAX_HEADER_BYTES,
                "column_count": None, "columns_preview": [],
                "columns_preview_truncated": None,
                "note": "Header exceeds bounded inspection; other tables are still audited."}
    fields = next(csv.reader([first.decode("utf-8-sig").rstrip("\r\n")],
                             delimiter="\t", strict=True))
    if not fields:
        raise ValueError("empty TSV header: " + name)
    # Some TSV exporters leave an index column's heading blank. Preserve its
    # original zero-based position; do not invent a feature name or shift IDs.
    # Duplicate named columns also require review, but cannot block inspection
    # of unrelated tables during this metadata-only stage.
    unnamed_positions = []
    duplicated_positions = []
    duplicate_examples = []
    seen_named = set()
    for position, field in enumerate(fields):
        if not field.strip():
            unnamed_positions.append(position)
        elif field in seen_named:
            duplicated_positions.append(position)
            if len(duplicate_examples) < 8:
                duplicate_examples.append(field)
        else:
            seen_named.add(field)
    flags = []
    if unnamed_positions:
        flags.append("unnamed_columns")
    if duplicated_positions:
        flags.append("duplicate_named_columns")
    return {
        **result,
        "status": "needs_review" if flags else "ok",
        "header_quality_flags": flags,
        "column_count": len(fields),
        "columns_preview": [f if f.strip() else None for f in fields[:16]],
        "columns_preview_truncated": len(fields) > 16,
        "unnamed_column_count": len(unnamed_positions),
        "unnamed_column_positions_preview": unnamed_positions[:16],
        "duplicate_named_column_count": len(duplicated_positions),
        "duplicate_named_columns_preview": duplicate_examples,
    }


def audit_archive(path: Path, verify_md5: bool = False) -> dict:
    """Inspect archive metadata and table headers; split membership NOT read."""
    if not path.is_file():
        raise ValueError("archive does not exist: " + str(path))
    with zipfile.ZipFile(path) as archive:
        members = archive.infolist()
        _check_zip_metadata(members)
        names = {item.filename for item in members if not item.is_dir()}
        missing_required = [name for name in REQUIRED_TABLES if BASE + name not in names]
        missing_other = [name for name in OTHER_OMICS if BASE + name not in names]
        splits = {}
        for study in STUDIES:
            expected_names = [f"{study}_all.txt"] + [
                f"{study}_split_{fold}_{kind}.txt"
                for fold in range(10) for kind in ("train", "val", "test")
            ]
            missing = [name for name in expected_names
                       if BASE + "splits/" + name not in names]
            splits[study] = {"expected": len(expected_names),
                             "present": len(expected_names) - len(missing),
                             "missing": missing}
        unknown_splits = []
        for name in names:
            if name.startswith(BASE + "splits/"):
                short = name[len(BASE + "splits/"):]
                if not SPLIT_PATTERN.fullmatch(short):
                    unknown_splits.append(short)
        tsvs = sorted(name for name in names if name.startswith(BASE)
                      and name.endswith(".tsv")
                      and ("/x_data/" in name or "/y_data/" in name))
        headers = [_header(archive, name) for name in tsvs]
        not_loaded = sorted(name.removeprefix(BASE) for name in names
                            if name.endswith((".pt", ".pkl", ".pickle", ".joblib")))
    report = {
        "source": "IMPROVE CSA Zenodo v2 / 10.5281/zenodo.15258883",
        "archive": str(path),
        "verification_scope": "ZIP metadata and first TSV header lines only",
        "file_entries": len(names),
        "response_table_present": BASE + "y_data/response.tsv" in names,
        "missing_required_tables": missing_required,
        "missing_other_feature_tables": missing_other,
        "split_manifests": splits,
        "unrecognized_split_filenames": sorted(unknown_splits),
        "table_headers": headers,
        "headers_exceeding_safety_cap": [h["table"] for h in headers
                                         if h["status"] == "header_exceeds_safety_cap"],
        "headers_needing_review": [h["table"] for h in headers
                                   if h["status"] == "needs_review"],
        "serialized_python_artifacts_not_loaded": not_loaded,
        "critical_limitations": [
            "Split-file memberships have not been read; no train/test overlap check yet.",
            "Header inspection does not establish response metric scale or missing-value semantics.",
            "Source/version cited; exact archive redistribution/training rights still pending.",
            "No archive contents extracted; do not deserialize .pt or scaler objects.",
        ],
    }
    if verify_md5:
        report["integrity"] = verify_published_md5(path)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--verify-md5", action="store_true")
    args = parser.parse_args(argv)
    try:
        report = audit_archive(args.archive, verify_md5=args.verify_md5)
    except (OSError, ValueError, zipfile.BadZipFile, UnicodeDecodeError) as exc:
        parser.exit(2, f"IMPROVE archive audit: {exc}\n")
    print("TUMORLAB / READ-ONLY IMPROVE ZENODO V2 AUDIT")
    print("Archive members:", report["file_entries"])
    print("Required tables missing:", report["missing_required_tables"])
    print("Other omics/drug features missing:", report["missing_other_feature_tables"])
    for study, summary in report["split_manifests"].items():
        print(f"  {study}: split filenames {summary['present']}/{summary['expected']}")
        if summary["missing"]:
            print("    missing:", summary["missing"])
    print("Unexpected split filenames:", report["unrecognized_split_filenames"])
    print("TSV headers (metadata only):")
    for item in report["table_headers"]:
        if item["status"] in ("ok", "needs_review"):
            print(f"  {item['table']}: {item['column_count']} columns; first names={item['columns_preview']}")
            if item["status"] == "needs_review":
                print(f"    REVIEW: {item['header_quality_flags']}; "
                      f"unnamed count={item['unnamed_column_count']} at first positions "
                      f"{item['unnamed_column_positions_preview']}; "
                      f"duplicate named count={item['duplicate_named_column_count']}; "
                      f"duplicate name examples={item['duplicate_named_columns_preview']}")
        else:
            print(f"  {item['table']}: HEADER EXCEEDS 32 MiB SAFETY CAP; columns not counted")
    if report["headers_needing_review"]:
        print("Table headers requiring source-specific review:", report["headers_needing_review"])
    if report["headers_exceeding_safety_cap"]:
        print("Wide headers needing follow-up:", report["headers_exceeding_safety_cap"])
    if args.verify_md5:
        print("Published checksum match:", report["integrity"]["matches_published_md5"])
    print("NOTE: Split contents, sample rows, licensing and patient outcomes NOT inspected; .pt never loaded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
