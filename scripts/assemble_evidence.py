"""Validate and assemble supplied evidence; no predictions, simulation or ranking.

The bundled example is synthetic. A file digest verifies unchanged bytes, not
chemical identity, scientific validity, authorship or a license to publish.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "medchem_evidence_v1"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def required_text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be nonempty text")
    return value


def number(value, label, minimum=None, integer=False):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number or null")
    if minimum is not None and value < minimum:
        raise ValueError(f"{label} must be at least {minimum}")
    if integer and not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")
    return value


def local_file(root, relative):
    relative = Path(required_text(relative, "structure_file"))
    if relative.is_absolute():
        raise ValueError("Structure paths must be relative to the input file")
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("Structure path escapes input directory")
    if not path.is_file():
        raise ValueError("Structure file missing")
    return path


def assemble(input_path):
    input_path = Path(input_path).resolve()
    def reject_nonfinite(value):
        raise ValueError(f"Invalid JSON numeric constant: {value}")
    data = json.loads(input_path.read_text(encoding="utf-8"), parse_constant=reject_nonfinite)
    if data.get("schema") != SCHEMA:
        raise ValueError(f"Expected schema {SCHEMA}")
    if data.get("scope") not in ("synthetic_example", "existing_evidence_only"):
        raise ValueError("Unsupported evidence scope")
    required_text(data.get("source_note"), "source_note")
    limitations = data.get("limitations")
    if not isinstance(limitations, list) or not limitations:
        raise ValueError("At least one limitation is required")
    for limitation in limitations:
        required_text(limitation, "limitation")
    targets = data.get("targets", [])
    if not isinstance(targets, list) or not targets:
        raise ValueError("At least one target is required")
    for target in targets:
        required_text(target, "target")
    if len(set(targets)) != len(targets):
        raise ValueError("Duplicate target")
    candidates = {}
    for candidate in data.get("candidates", []):
        name = required_text(candidate.get("id"), "candidate id")
        if name in candidates:
            raise ValueError("Duplicate candidate id")
        path = local_file(input_path.parent, candidate.get("structure_file"))
        if sha(path) != candidate.get("structure_sha256"):
            raise ValueError(f"Structure file SHA mismatch: {name}")
        required_text(candidate.get("smiles"), "smiles")
        required_text(candidate.get("source"), "candidate source")
        candidates[name] = candidate
    if not candidates:
        raise ValueError("At least one candidate is required")
    evidence = {}
    for item in data.get("evidence", []):
        key = (item.get("candidate"), item.get("target"))
        if key[0] not in candidates or key[1] not in targets:
            raise ValueError("Evidence references unknown candidate or target")
        if key in evidence:
            raise ValueError("Duplicate candidate-target evidence")
        required_text(item.get("source"), "evidence source")
        evidence[key] = item
    expected = {(name, target) for name in candidates for target in targets}
    if set(evidence) != expected:
        raise ValueError("Supply one explicit record for every candidate-target pair, including missing evidence")
    rows = []
    for name, candidate in candidates.items():
        props = candidate.get("properties", {})
        for target in targets:
            item = evidence[(name, target)]
            md = item.get("md", {})
            status = md.get("status")
            if status not in ("not_available", "incomplete_snapshot", "completed_reviewed"):
                raise ValueError("MD status must explicitly identify missing, provisional or completed reviewed evidence")
            md_ns = number(md.get("sampling_ns"), "sampling_ns", minimum=0)
            md_seed = number(md.get("seed"), "seed", minimum=1, integer=True)
            md_rmsd = number(md.get("ligand_rmsd_mean_A"), "ligand_rmsd_mean_A", minimum=0)
            if status == "not_available" and any(v is not None for v in (md_ns, md_seed, md_rmsd)):
                raise ValueError("Missing MD must not contain apparent measured metrics")
            if status != "not_available":
                required_text(md.get("source"), "MD source")
                if md_ns is None or md_ns <= 0 or md_seed is None or md_rmsd is None:
                    raise ValueError("Present MD requires a positive duration, seed and finite RMSD")
            rows.append({
                "candidate_id": name, "target": target, "SMILES": candidate["smiles"],
                "scope": data["scope"], "source": item["source"],
                "supplied_docking_score_kcal_mol": number(item.get("docking_score_kcal_mol"), "docking score"),
                "supplied_MW": number(props.get("MW"), "MW", minimum=0),
                "supplied_cLogP": number(props.get("cLogP"), "cLogP"),
                "supplied_TPSA_A2": number(props.get("TPSA_A2"), "TPSA", minimum=0),
                "supplied_ADMET_raw_score": number(props.get("ADMET_raw_score"), "ADMET raw score"),
                "ADMET_endpoint": props.get("ADMET_endpoint"),
                "MD_evidence_status": status, "MD_seed": md_seed,
                "MD_sampling_ns": md_ns, "MD_ligand_RMSD_mean_A": md_rmsd,
                "measured_affinity": None, "human_decision": "not_recorded",
            })
    return data, rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "examples" / "synthetic_evidence.json")
    parser.add_argument("--out", type=Path, default=ROOT / "results" / "demo")
    args = parser.parse_args(argv)
    data, rows = assemble(args.input)
    if args.out.exists() and any(args.out.iterdir()):
        raise ValueError("Output directory must be new or empty; preserve earlier receipts")
    args.out.mkdir(parents=True, exist_ok=True)
    table = args.out / "candidates.csv"
    with table.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "schema": SCHEMA, "scope": data["scope"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_note": data["source_note"], "candidate_target_rows": rows,
        "limitations": data["limitations"], "new_predictions": False, "new_MD": False,
        "automatic_selection": False,
    }
    output = args.out / "results.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    receipt = {
        "scope": data["scope"], "input_sha256": sha(args.input),
        "script_sha256": sha(__file__), "rows": len(rows),
        "structure_files_verified": len(data["candidates"]),
        "file_hash_verification_is_not_chemical_identity_validation": True,
        "output_sha256": {p.name: sha(p) for p in (table, output)},
        "new_MD": False, "new_AI_inference": False,
    }
    (args.out / "execution_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return receipt


if __name__ == "__main__":
    main()
