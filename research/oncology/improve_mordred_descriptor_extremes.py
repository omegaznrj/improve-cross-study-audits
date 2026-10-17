"""Read-only Phase 1.5D: identify Mordred feature columns with extreme magnitude.

Exports descriptor NAMES and aggregate ranges only. No chemical IDs or row
values, response labels, fit or source modifications. Names are source-provided,
not proof of descriptor units or correctness.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
from pathlib import Path
import zipfile
from research.oncology.improve_split_audit import BASE, _read_response
from research.oncology.improve_numeric_quality_audit import _value_kind


def extreme_descriptors(z: zipfile.ZipFile, expected: set[str], threshold=1e6):
    if threshold <= 0 or not math.isfinite(threshold):
        raise ValueError("threshold must be positive finite")
    member = BASE+"x_data/drug_mordred.tsv"
    with z.open(member) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t", strict=True)
        header = next(reader)
        if header[0] != "improve_chem_id" or len(header) != len(set(header)):
            raise ValueError("invalid Mordred column headings")
        descriptors = header[1:]
        mins = [math.inf]*len(descriptors)
        maxs = [-math.inf]*len(descriptors)
        count = [0]*len(descriptors)
        above = [0]*len(descriptors)
        seen=set()
        for row in reader:
            if len(row) != len(header):
                raise ValueError("invalid Mordred row width")
            key = row[0].strip()
            if key not in expected:
                continue
            if key in seen:
                raise ValueError("duplicate linked Mordred key")
            seen.add(key)
            for i, rawval in enumerate(row[1:]):
                kind, value = _value_kind(rawval)
                if kind != "valid":
                    continue
                count[i] += 1
                mins[i] = min(mins[i], value)
                maxs[i] = max(maxs[i], value)
                if abs(value) >= threshold:
                    above[i] += 1
    return {
        "linked_drugs": len(seen),
        "threshold_absolute": threshold,
        "extreme_descriptor_columns": [
            {"descriptor": name, "minimum": mins[i], "maximum": maxs[i],
             "valid_linked_drugs": count[i], "values_at_or_above_threshold": above[i]}
            for i,name in enumerate(descriptors) if above[i]
        ],
    }


def audit(path: Path, threshold=1e6):
    with zipfile.ZipFile(path) as z:
        if len(z.namelist()) != len(set(z.namelist())):
            raise ValueError("duplicate ZIP members")
        drug_ids=set(_read_response(z)["drug_code_by_id"])
        return {
            "schema":"tumorlab-improve-mordred-extremes-v1",
            **extreme_descriptors(z,drug_ids,threshold),
            "outcomes_read":False,"training_performed":False,
            "notes":[
                "Column names and observed ranges are descriptive, not chemically validated.",
                "Different Mordred descriptor definitions have different scales.",
                "Source-level feature definitions and generator versions must be verified before preprocessing.",
                "No compound identifiers or compound-specific values displayed.",
            ]
        }


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive",type=Path,default=Path(__file__).resolve().parents[2]/"data/raw/improve_cross_study/csa_data.zip")
    parser.add_argument("--threshold",type=float,default=1e6)
    args=parser.parse_args(argv)
    try:
        print(json.dumps(audit(args.archive,args.threshold),indent=2))
    except (OSError,ValueError,zipfile.BadZipFile,csv.Error) as exc:
        parser.exit(2,"IMPROVE Mordred extreme names: "+str(exc)+"\n")


if __name__=="__main__":
    main()
