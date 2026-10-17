"""Read-only IMPROVE chemical-identity concordance and cross-study overlap audit.

Compares text records for training-associated CCLE/gCSI compounds, without
asserting that differing SMILES imply differing chemical structures.
Quantifies possible cross-study duplicated biological pairs and identical
ECFP fingerprints. Uses no response values, model fitting, or data export.
"""
from __future__ import annotations

import argparse
import csv
import io
from collections import Counter, defaultdict
from pathlib import Path
import zipfile

from research.oncology.improve_split_audit import BASE, _indices, _read_response

STUDIES = ("CCLE", "gCSI")
PARTITIONS = ("train", "val", "test")
MAX_ROWS = 2_000_000


def _read_tsv(archive: zipfile.ZipFile, filename: str, required: tuple[str, ...]):
    path = BASE + "x_data/" + filename
    if path not in archive.namelist():
        raise ValueError("missing chemical table: " + filename)
    with archive.open(path) as raw, io.TextIOWrapper(
        raw, encoding="utf-8-sig", newline=""
    ) as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if (reader.fieldnames is None or
                len(set(reader.fieldnames)) != len(reader.fieldnames) or
                not all(name in reader.fieldnames for name in required)):
            raise ValueError("invalid chemical table headings: " + filename)
        if filename == "drug_ecfp4_nbits512.tsv":
            expected = ["improve_chem_id"] + [f"ecfp4.{i}" for i in range(512)]
            if reader.fieldnames != expected:
                raise ValueError("fingerprint table does not have 512 ordered ECFP4 bits")
        for index, row in enumerate(reader):
            if index >= MAX_ROWS or None in row or any(x is None for x in row.values()):
                raise ValueError("unexpected chemical table structure: " + filename)
            yield row


def _study_sets(archive: zipfile.ZipFile, rows: dict, fold: int) -> dict:
    result = {}
    for study in STUDIES:
        parts = {}
        for part in PARTITIONS:
            indices = _indices(archive, f"{study}_split_{fold}_{part}.txt", rows["rows"])
            parts[part] = {
                "samples": {rows["samples"][i] for i in indices},
                "drugs": {rows["drugs"][i] for i in indices},
                "pairs": {(rows["samples"][i], rows["drugs"][i]) for i in indices},
                "rows": len(indices),
            }
        parts["all"] = {
            "samples": set().union(*(parts[x]["samples"] for x in PARTITIONS)),
            "drugs": set().union(*(parts[x]["drugs"] for x in PARTITIONS)),
            "pairs": set().union(*(parts[x]["pairs"] for x in PARTITIONS)),
            "rows": sum(parts[x]["rows"] for x in PARTITIONS),
        }
        result[study] = parts
    return result


def _overlap(one: dict, other: dict) -> dict:
    return {key: len(one[key] & other[key])
            for key in ("samples", "drugs", "pairs")}


def _chemical_records(archive: zipfile.ZipFile,
                      training_drug_keys: set[str]) -> dict:
    reference: dict[str, set[str]] = defaultdict(set)
    primary_counts: Counter[str] = Counter()
    blank_reference = 0
    for row in _read_tsv(archive, "drug_SMILES.tsv",
                         ("improve_chem_id", "canSMILES")):
        key = (row["improve_chem_id"] or "").strip()
        if key not in training_drug_keys:
            continue
        primary_counts[key] += 1
        value = (row["canSMILES"] or "").strip()
        if value:
            reference[key].add(value)
        else:
            blank_reference += 1

    aliases: dict[str, set[str]] = defaultdict(set)
    pubchem: dict[str, set[str]] = defaultdict(set)
    alias_counts: Counter[str] = Counter()
    for row in _read_tsv(archive, "drug_info.tsv",
                         ("improve_chem_id", "canSMILES", "PUBCHEM_ID")):
        key = (row["improve_chem_id"] or "").strip()
        if key not in training_drug_keys:
            continue
        alias_counts[key] += 1
        val = (row["canSMILES"] or "").strip()
        if val:
            aliases[key].add(val)
        pub = (row["PUBCHEM_ID"] or "").strip()
        if pub:
            pubchem[key].add(pub)

    fingerprint_counts: Counter[str] = Counter()
    vectors: dict[str, tuple[str, ...]] = {}
    vector_duplicate_key_rows = 0
    for row in _read_tsv(archive, "drug_ecfp4_nbits512.tsv",
                         ("improve_chem_id", "ecfp4.0")):
        key = (row["improve_chem_id"] or "").strip()
        if key not in training_drug_keys:
            continue
        fingerprint_counts[key] += 1
        if key in vectors:
            vector_duplicate_key_rows += 1
            continue
        vectors[key] = tuple(value.strip() for field, value in row.items()
                             if field != "improve_chem_id")
    fingerprint_groups: dict[tuple[str, ...], int] = Counter(vectors.values())
    reference_groups: dict[str, set[str]] = defaultdict(set)
    for key in training_drug_keys:
        if len(reference[key]) == 1:
            reference_groups[next(iter(reference[key]))].add(key)

    result = {
        "distinct_training_chemical_ids": len(training_drug_keys),
        "drug_smiles": {
            "ids_with_rows": len(primary_counts),
            "ids_with_nonblank_canonical_smiles": sum(bool(reference[k])
                                                       for k in training_drug_keys),
            "ids_with_more_than_one_reference_smiles": sum(len(reference[k]) > 1
                                                           for k in training_drug_keys),
            "duplicate_reference_id_rows": sum(n - 1 for n in primary_counts.values()),
            "blank_reference_smiles_rows": blank_reference,
            "canonical_text_shared_by_multiple_ids": sum(
                len(group) > 1 for group in reference_groups.values()
            ),
        },
        "drug_info_comparison": {
            "ids_with_metadata_rows": len(alias_counts),
            "ids_with_repeated_metadata_rows": sum(n > 1 for n in alias_counts.values()),
            "ids_with_multiple_distinct_nonblank_alias_smiles": sum(
                len(aliases[key]) > 1 for key in training_drug_keys
            ),
            "ids_with_multiple_nonblank_pubchem_ids": sum(
                len(pubchem[key]) > 1 for key in training_drug_keys
            ),
            "ids_with_any_exact_alias_reference_match": sum(
                bool(reference[key] & aliases[key])
                for key in training_drug_keys
            ),
            "ids_with_alias_but_no_exact_reference_text_match": sum(
                bool(aliases[key]) and bool(reference[key]) and
                not bool(reference[key] & aliases[key])
                for key in training_drug_keys
            ),
            "ids_without_nonblank_alias_canonical_text": sum(
                not bool(aliases[key]) for key in training_drug_keys
            ),
        },
        "fingerprints": {
            "ids_with_rows": len(fingerprint_counts),
            "duplicate_key_rows": vector_duplicate_key_rows,
            "identical_512bit_fingerprint_groups_across_distinct_ids": sum(
                n > 1 for n in fingerprint_groups.values()
            ),
            "distinct_ids_in_identical_fingerprint_groups": sum(
                n for n in fingerprint_groups.values() if n > 1
            ),
        },
        "interpretation": [
            "Exact SMILES-string differences are not proof of different molecules.",
            "Duplicate compound IDs and identical fingerprints require chemical crosswalk review.",
            "No RDKit normalization, chemical-identity verification or model training performed.",
        ],
    }
    return result


def audit_chemical_and_cross_study(path: Path, fold: int = 0) -> dict:
    if type(fold) is not int or fold not in range(10):
        raise ValueError("fold must be 0..9")
    if not path.is_file():
        raise ValueError("ZIP archive missing: " + str(path))
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate ZIP member paths")
        rows = _read_response(archive)
        studies = _study_sets(archive, rows, fold)
        training_codes = (
            studies["CCLE"]["train"]["drugs"] |
            studies["gCSI"]["train"]["drugs"]
        )
        selected_ids = {
            key for key, code in rows["drug_code_by_id"].items()
            if code in training_codes
        }
        chemicals = _chemical_records(archive, selected_ids)
        overlap = {
            "ccle_all_vs_gcsi_all": _overlap(
                studies["CCLE"]["all"], studies["gCSI"]["all"]),
            "ccle_train_vs_gcsi_test": _overlap(
                studies["CCLE"]["train"], studies["gCSI"]["test"]),
            "gcsi_train_vs_ccle_test": _overlap(
                studies["gCSI"]["train"], studies["CCLE"]["test"]),
            "ccle_train_vs_gcsi_val": _overlap(
                studies["CCLE"]["train"], studies["gCSI"]["val"]),
            "gcsi_train_vs_ccle_val": _overlap(
                studies["gCSI"]["train"], studies["CCLE"]["val"]),
        }
        training_by_study = {
            study: {
                "rows": studies[study]["train"]["rows"],
                "unique_samples": len(studies[study]["train"]["samples"]),
                "unique_drugs": len(studies[study]["train"]["drugs"]),
                "unique_pairs": len(studies[study]["train"]["pairs"]),
            }
            for study in STUDIES
        }
    return {
        "source": "IMPROVE Zenodo v2, cell-line screens",
        "fold": fold, "training_by_study": training_by_study,
        "cross_study_overlap": overlap,
        "chemical_text_and_fingerprint_diagnostics": chemicals,
        "no_response_target_values_used": True,
        "no_model_training": True,
        "limitations": [
            "Cross-study identical harmonized IDs do not guarantee identical biological samples.",
            "Cross-study overlap complicates independent-study and unseen-pair claims.",
            "String matches/mismatches do not establish true chemical equivalence or differences.",
            "ECFP collisions can occur; no inferred chemical resolution or filtering is done.",
            "No data mutation, raw compound names/strings output, or model fitting.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path,
                        default=Path(__file__).resolve().parents[2] /
                        "data/raw/improve_cross_study/csa_data.zip")
    parser.add_argument("--fold", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        r = audit_chemical_and_cross_study(args.archive, args.fold)
    except (OSError, ValueError, UnicodeDecodeError, csv.Error, zipfile.BadZipFile) as exc:
        parser.exit(2, "IMPROVE chemical audit: " + str(exc) + "\n")
    print("TUMORLAB / IMPROVE CHEMICAL AND CROSS-STUDY AUDIT (NO FIT)")
    print("Fold:", r["fold"])
    print("Training study entity counts:", r["training_by_study"])
    print("Cross-study sample/drug/pair overlaps:")
    for key, val in r["cross_study_overlap"].items():
        print("  " + key + ":", val)
    print("Chemical source text and ECFP concordance:")
    for name, values in r["chemical_text_and_fingerprint_diagnostics"].items():
        if isinstance(values, dict):
            print("  " + name + ":", values)
    print("No outcome measurements read for modeling; no chemical identities adjudicated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
