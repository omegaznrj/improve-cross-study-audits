"""Read-only IMPROVE repeated-pair provenance and metadata availability audit.

No outcome numbers, identifiers, or source row content in output. Counts are
descriptive: shared cell-drug pairs do not establish shared assay occasions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import io
import json
from pathlib import Path
import zipfile

from research.oncology.improve_archive_audit import DEFAULT_ARCHIVE
from research.oncology.improve_split_audit import BASE, STUDIES, _indices

PARTS = ("train", "val", "test")
META_CANDIDATES = (
    "source", "study", "experiment_id", "assay_id", "screen_id",
    "replicate", "replicate_id", "batch", "batch_id", "dose",
    "concentration", "time", "duration", "endpoint", "response_metric",
    "improve_sample_id", "improve_chem_id",
)
LIMIT = 2_000_000


def audit_repeat_provenance(path: Path = DEFAULT_ARCHIVE, fold: int = 0) -> dict:
    if type(fold) is not int or not 0 <= fold <= 9:
        raise ValueError("fold must be 0..9")
    if not path.is_file():
        raise ValueError("archive not found")
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate archive member paths")
        member = BASE + "y_data/response.tsv"
        if member not in names:
            raise ValueError("response table missing")
        split_sets = {}
        for study in STUDIES:
            all_idx = _indices(z, f"{study}_all.txt", LIMIT)
            parts = {p: _indices(z, f"{study}_split_{fold}_{p}.txt", LIMIT) for p in PARTS}
            if set.union(*parts.values()) != all_idx:
                raise ValueError(study + ": split union mismatch")
            if any(parts[a] & parts[b] for a in PARTS for b in PARTS if a < b):
                raise ValueError(study + ": split row overlap")
            split_sets[study] = (all_idx, parts)
        # The original row index remains the sole split locator.
        header = []
        counts = {s: Counter() for s in STUDIES}
        pairs = {s: defaultdict(lambda: Counter()) for s in STUDIES}
        source_study = Counter()
        source_mismatch = Counter()
        study_values = {s: Counter() for s in STUDIES}
        metadata_fields = []
        source_vocab = {s: set() for s in STUDIES}
        processed = 0
        with z.open(member) as raw, io.TextIOWrapper(
            raw, encoding="utf-8-sig", newline=""
        ) as stream:
            reader = csv.DictReader(stream, delimiter="\t", strict=True)
            if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ValueError("invalid response header")
            header = reader.fieldnames
            required = {"source", "study", "improve_sample_id", "improve_chem_id"}
            if not required.issubset(header):
                raise ValueError("missing required identifier fields")
            # Use field names to establish metadata *availability*, not interpret
            # numeric assay targets or invent absent protocol descriptors.
            metadata_fields = [c for c in META_CANDIDATES if c in header]
            membership = {}
            for study, (all_idx, parts) in split_sets.items():
                for i in all_idx:
                    if i in membership:
                        raise ValueError("source all sets share original row index")
                    membership[i] = (study, next(
                        (p for p in PARTS if i in parts[p]), None))
            for i, row in enumerate(reader):
                if i >= LIMIT or None in row or any(v is None for v in row.values()):
                    raise ValueError("invalid response row")
                processed += 1
                if i not in membership:
                    continue
                study, part = membership[i]
                source = (row["source"] or "").strip()
                row_study = (row["study"] or "").strip()
                source_vocab[study].add(source)
                study_values[study][row_study] += 1
                if source != study:
                    source_mismatch[study] += 1
                if row_study != study:
                    source_study[study] += 1
                sample = (row["improve_sample_id"] or "").strip()
                drug = (row["improve_chem_id"] or "").strip()
                if not sample or not drug:
                    counts[study]["blank_join_key_rows"] += 1
                    continue
                key = (sample, drug)
                pairs[study][key][part] += 1
                counts[study]["rows"] += 1
        if processed != sum(len(x[0]) for x in split_sets.values()):
            raise ValueError("response row total differs from study manifest union")
        result = {}
        for study in STUDIES:
            frequencies = pairs[study]
            repeated = {k: v for k, v in frequencies.items() if sum(v.values()) > 1}
            hist = Counter()
            for v in repeated.values():
                n = sum(v.values())
                band = "2" if n == 2 else "3-5" if n <= 5 else "6-10" if n <= 10 else "11+"
                hist[band] += 1
            result[study] = {
                "rows": counts[study]["rows"],
                "source_label_mismatches": source_mismatch[study],
                "study_label_not_equal_to_source_name": source_study[study],
                "distinct_study_field_values": len(study_values[study]),
                "study_field_value_multiplicity_histogram": dict(sorted(
                    Counter(study_values[study].values()).items()
                )),
                "blank_study_field_rows": study_values[study].get("", 0),
                "distinct_source_labels": len(source_vocab[study]),
                "blank_join_key_rows": counts[study]["blank_join_key_rows"],
                "unique_pairs": len(frequencies),
                "repeated_pair_groups": len(repeated),
                "repeated_pair_excess_rows": sum(sum(v.values()) - 1 for v in repeated.values()),
                "largest_pair_multiplicity": max((sum(v.values()) for v in frequencies.values()), default=0),
                "repeated_group_size_histogram": dict(sorted(hist.items())),
                "repeat_groups_spanning_partitions": sum(
                    sum(bool(v[p]) for p in PARTS) >= 2 for v in repeated.values()
                ),
                "repeat_groups_entirely_in_one_partition": sum(
                    sum(bool(v[p]) for p in PARTS) == 1 for v in repeated.values()
                ),
                "possible_assay_metadata_columns": [
                    k for k in metadata_fields if k not in required
                ],
                "confirmed_assay_occasion_identity": False,
            }
        return {
            "schema": "tumorlab-improve-repeat-provenance-audit-v1",
            "fold": fold,
            "response_rows": processed,
            "present_response_columns": header,
            "studies": result,
            "training_performed": False,
            "response_values_interpreted": False,
            "raw_biological_identifiers_exported": False,
            "limitations": [
                "Repeated cell-drug pair is not proof of same original experiment or replicate.",
                "The study field may identify experiments or other study-level groupings; definition remains unverified.",
                "Study field values and biological identifiers are not included in the output; only frequencies of their multiplicities.",
                "Available columns do not establish assay context semantics; consult originating studies.",
                "No response values or experimental occasion identifiers are adjudicated.",
                "Original partition assignments preserved; no models trained.",
            ],
        }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument("--fold", type=int, default=0)
    args = ap.parse_args(argv)
    try:
        print(json.dumps(audit_repeat_provenance(args.archive, args.fold), indent=2))
    except (OSError, ValueError, zipfile.BadZipFile, csv.Error, UnicodeError) as exc:
        ap.exit(2, "IMPROVE repeat provenance audit: " + str(exc) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
