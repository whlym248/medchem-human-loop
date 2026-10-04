"""Freeze supplied evidence into a review packet and record human decisions.

This coordinates evidence, local case search and transparent route suggestions.
It never launches MD, docking, inference, a generator, or a remote job. A decision
is a local attestation, not authenticated identity or scientific validation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import shutil
import tempfile

from assemble_evidence import assemble, local_file, required_text
from run_md_v2 import run_lock

SCHEMA = "medchem_review_round_v1"
ACTIONS = ("approve_next_round", "reject", "request_evidence", "defer")
METRICS = ("supplied_docking_score_kcal_mol", "supplied_MW", "supplied_cLogP", "supplied_TPSA_A2", "supplied_ADMET_raw_score")
LIMITATIONS = [
    "Supplied evidence and review labels are not independently scientifically authenticated.",
    "Displayed differences are arithmetic only; protocols and applicability must be reviewed before interpretation.",
    "Missing values remain missing. No affinity gain, clinical risk, convergence or synthesis feasibility is inferred.",
    "A human next-round decision is not permission to bypass MD pilot review or to launch paid resources.",
    "No candidate is automatically ranked, selected or submitted. Diffusion and a learned optimization policy are not implemented.",
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    def reject(value):
        raise ValueError(f"Nonfinite JSON value: {value}")
    return json.loads(Path(path).read_text(encoding="utf-8-sig"), parse_constant=reject)


def write_new(path, value):
    with Path(path).open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def render(packet):
    esc = lambda x: html.escape("missing" if x is None else str(x), quote=True)
    columns = ["candidate_id", "target", "source", "ADMET_endpoint", "MD_evidence_status",
               "MD_seed", "MD_sampling_ns", "MD_ligand_RMSD_mean_A", *METRICS]
    header = "".join(f"<th>{esc(key)}</th>" for key in columns)
    rows = "".join("<tr>" + "".join(f"<td>{esc(row.get(key))}</td>" for key in columns) + "</tr>" for row in packet["rows"])
    limits = "".join(f"<li>{esc(item)}</li>" for item in packet["limitations"])
    diffs = esc(json.dumps(packet["comparisons"], ensure_ascii=False, indent=2))
    context = esc(json.dumps(packet["planning_context"], ensure_ascii=False, indent=2))
    return f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Review packet</title>
<style>body{{font:16px system-ui;margin:32px;line-height:1.5}}table{{border-collapse:collapse;display:block;overflow:auto}}th,td{{border:1px solid #bbb;padding:8px}}pre{{white-space:pre-wrap}}.scope{{background:#fff1bf;padding:12px}}</style>
<h1>Human review packet</h1><p class="scope">Scope: {esc(packet['scope'])}. Supplied evidence, not new scientific validation.</p>
<p>Parent: {esc(packet['parent_id'])}. Created UTC: {esc(packet['created_utc'])}.</p>
<p><a href="review_packet.json">Full review packet</a> · <a href="evidence/input.json">Input, structures and source declarations</a></p>
<ul>{limits}</ul><table><tr>{header}</tr>{rows}</table>
<h2>Arithmetic differences from supplied parent values</h2><pre>{diffs}</pre>
<h2>Case search and route suggestion</h2><pre>{context}</pre>
<p>Record decisions with workflow_round.py decide. This HTML is static and does not execute jobs.</p></html>'''


def create_round(input_path, out_dir, parent_id, *, case_path=None, query="", route_path=None, previous_round=None):
    source = Path(input_path).resolve()
    out = Path(out_dir).resolve()
    if out.exists():
        raise ValueError("Round output must not already exist")
    source_digest = sha(source)
    data, rows = assemble(source)
    if sha(source) != source_digest:
        raise ValueError("Evidence changed while being read")
    ids = [item["id"] for item in data["candidates"]]
    if parent_id not in ids:
        raise ValueError("Parent must name a supplied candidate")
    ancestor = None
    if previous_round is not None:
        previous = inspect_round(previous_round)
        if previous["scope"] != data["scope"]:
            raise ValueError("Cannot change synthetic/real evidence scope between rounds")
        if not any(item["action"] == "approve_next_round" for item in previous["latest_decisions"].values()):
            raise ValueError("Previous round has no human approve_next_round decision")
        approved = {key for key, item in previous["latest_decisions"].items() if item["action"] == "approve_next_round"}
        if parent_id not in approved:
            raise ValueError("New parent must be explicitly approved in the previous round")
        old_data, _ = assemble(Path(previous_round) / "evidence" / "input.json")
        old_candidate = next(item for item in old_data["candidates"] if item["id"] == parent_id)
        new_candidate = next(item for item in data["candidates"] if item["id"] == parent_id)
        if old_candidate["smiles"] != new_candidate["smiles"] or old_candidate["structure_sha256"] != new_candidate["structure_sha256"]:
            raise ValueError("Approved parent identity changed between rounds")
        ancestor = {"manifest_sha256": previous["manifest_sha256"], "decision_chain_tip": previous["decision_chain_tip"], "parent_id": parent_id}
    parent_rows = {r["target"]: r for r in rows if r["candidate_id"] == parent_id}
    comparisons = []
    for row in rows:
        if row["candidate_id"] == parent_id:
            continue
        base = parent_rows[row["target"]]
        differences = {m: None if row[m] is None or base[m] is None else row[m] - base[m] for m in METRICS}
        same_endpoint = bool(row.get("ADMET_endpoint")) and row.get("ADMET_endpoint") == base.get("ADMET_endpoint")
        if not same_endpoint:
            differences["supplied_ADMET_raw_score"] = None
        comparisons.append({"candidate_id": row["candidate_id"], "target": row["target"],
            "delta_candidate_minus_parent": differences, "ADMET_endpoint_labels_match": same_endpoint,
            "comparable_protocols_verified": False, "benefit_or_risk_classification": "requires_human_interpretation"})
    planning = {"case_search": None, "route_suggestion": None}
    if case_path is not None:
        from case_library import search_cases
        planning["case_search"] = search_cases(case_path, query=query, allow_synthetic=data["scope"] == "synthetic_example")
    if route_path is not None:
        from route_proposal import recommend
        request = read(route_path)
        matches = (planning["case_search"] or {}).get("matches", [])
        # Synthetic fixtures must never become real-world case support.
        request["case_support"] = any(item.get("scope") == "literature_curated" for item in matches)
        planning["route_suggestion"] = recommend(request)
    packet = {"schema": SCHEMA, "scope": data["scope"], "created_utc": now(), "parent_id": parent_id,
              "candidate_ids": ids, "rows": rows, "comparisons": comparisons,
              "planning_context": planning, "previous_round": ancestor,
              "limitations": data["limitations"] + LIMITATIONS,
              "new_MD": False, "new_inference": False, "automatic_selection": False}
    out.parent.mkdir(parents=True, exist_ok=True)
    # A failed preparation leaves no apparently complete round at the final path.
    with tempfile.TemporaryDirectory(prefix=".round-preparing-", dir=out.parent) as temporary:
        stage = Path(temporary) / "round"
        evidence = stage / "evidence"
        evidence.mkdir(parents=True)
        for candidate in data["candidates"]:
            original = local_file(source.parent, candidate["structure_file"])
            relative = Path(candidate["structure_file"])
            if relative.is_absolute() or relative.drive or ".." in relative.parts:
                raise ValueError("Snapshot structure paths must stay inside evidence without parent traversal")
            target = (evidence / relative).resolve()
            if not target.is_relative_to(evidence.resolve()) or target == evidence.resolve():
                raise ValueError("Snapshot structure path escapes evidence directory")
            if target == evidence / "input.json":
                raise ValueError("Structure filename input.json is reserved for the snapshot")
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                shutil.copyfile(original, target)
        write_new(evidence / "input.json", data)
        assemble(evidence / "input.json")  # Recheck the frozen bytes, not the original filenames.
        if case_path is not None:
            shutil.copyfile(case_path, stage / "case_library_snapshot.jsonl")
            # Reject case files changed while building the packet.
            if sha(stage / "case_library_snapshot.jsonl") != planning["case_search"]["library_sha256"]:
                raise ValueError("Case library bytes changed during snapshot")
            from case_library import search_cases
            if search_cases(stage / "case_library_snapshot.jsonl", query=query, allow_synthetic=data["scope"] == "synthetic_example")["matches"] != planning["case_search"]["matches"]:
                raise ValueError("Case library changed during snapshot")
        if route_path is not None:
            write_new(stage / "route_request.json", request)
        write_new(stage / "review_packet.json", packet)
        (stage / "report.html").write_text(render(packet), encoding="utf-8", newline="\n")
        files = {p.relative_to(stage).as_posix(): sha(p) for p in sorted(stage.rglob("*")) if p.is_file()}
        write_new(stage / "round_manifest.json", {"schema": SCHEMA, "created_utc": packet["created_utc"],
                  "input_sha256": source_digest, "snapshot_files": files, "script_sha256": sha(__file__)})
        stage.rename(out)
    return inspect_round(out)


def verify_snapshot(directory):
    root = Path(directory).resolve()
    manifest = read(root / "round_manifest.json")
    if manifest.get("schema") != SCHEMA:
        raise ValueError("Unknown round schema")
    files = manifest.get("snapshot_files", {})
    if not all(name in files for name in ("review_packet.json", "report.html", "evidence/input.json")):
        raise ValueError("Snapshot manifest lacks required artifacts")
    for relative, expected in files.items():
        path = local_file(root, relative)
        if sha(path) != expected:
            raise ValueError(f"Frozen evidence changed: {relative}")
    packet = read(root / "review_packet.json")
    if packet.get("schema") != SCHEMA:
        raise ValueError("Unknown packet schema")
    assemble(root / "evidence" / "input.json")
    return root, packet, sha(root / "round_manifest.json")


def inspect_round(directory):
    root, packet, manifest_sha = verify_snapshot(directory)
    events = sorted((root / "decisions").glob("*.json"))
    latest, previous = {}, None
    for sequence, path in enumerate(events, 1):
        if path.name != f"decision_{sequence:04d}.json":
            raise ValueError("Decision sequence has a gap or unrecognized file")
        event = read(path)
        if event.get("schema") != SCHEMA or event.get("scope") != packet["scope"] or event.get("limitations_acknowledged") is not True:
            raise ValueError("Decision schema, scope or acknowledgement mismatch")
        if (event.get("round_manifest_sha256") != manifest_sha or event.get("previous_event_sha256") != previous
                or event.get("sequence") != sequence):
            raise ValueError("Decision history identity or chain mismatch")
        if event.get("action") not in ACTIONS or event.get("candidate_id") not in packet["candidate_ids"]:
            raise ValueError("Invalid recorded decision")
        required_text(event.get("reviewer"), "reviewer")
        required_text(event.get("reason"), "reason")
        latest[event["candidate_id"]] = event
        previous = sha(path)
    return {"schema": SCHEMA, "scope": packet["scope"], "manifest_sha256": manifest_sha,
            "decision_chain_tip": previous, "events": len(events), "latest_decisions": latest,
            "awaiting_human_decision": [key for key in packet["candidate_ids"] if key not in latest],
            "jobs_launched": 0, "reviewer_identity_authenticated": False}


def record_decision(directory, candidate_id, action, reviewer, reason, *, acknowledge_limitations=False):
    if action not in ACTIONS:
        raise ValueError("Unknown decision action")
    required_text(reviewer, "reviewer")
    required_text(reason, "reason")
    if acknowledge_limitations is not True:
        raise ValueError("Explicit acknowledgement of evidence limitations is required")
    root, packet, _ = verify_snapshot(directory)
    if candidate_id not in packet["candidate_ids"]:
        raise ValueError("Unknown candidate")
    with run_lock(root):
        state = inspect_round(root)
        event = {"schema": SCHEMA, "sequence": state["events"] + 1, "recorded_utc": now(),
                 "round_manifest_sha256": state["manifest_sha256"], "previous_event_sha256": state["decision_chain_tip"],
                 "scope": packet["scope"], "candidate_id": candidate_id, "action": action,
                 "reviewer": reviewer, "reason": reason, "limitations_acknowledged": True,
                 "reviewer_identity_authenticated": False, "scientific_validation": False, "jobs_launched": 0}
        folder = root / "decisions"
        folder.mkdir(exist_ok=True)
        write_new(folder / f"decision_{event['sequence']:04d}.json", event)
    return event


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare")
    prepare.add_argument("--input", type=Path, required=True)
    prepare.add_argument("--out", type=Path, required=True)
    prepare.add_argument("--parent", required=True)
    prepare.add_argument("--cases", type=Path)
    prepare.add_argument("--query", default="")
    prepare.add_argument("--route-request", type=Path)
    prepare.add_argument("--previous-round", type=Path)
    status = sub.add_parser("status")
    status.add_argument("--round", type=Path, required=True)
    decision = sub.add_parser("decide")
    decision.add_argument("--round", type=Path, required=True)
    decision.add_argument("--candidate", required=True)
    decision.add_argument("--action", choices=ACTIONS, required=True)
    decision.add_argument("--reviewer", required=True)
    decision.add_argument("--reason", required=True)
    decision.add_argument("--acknowledge-limitations", action="store_true")
    args = parser.parse_args(argv)
    if args.command == "prepare":
        result = create_round(args.input, args.out, args.parent, case_path=args.cases, query=args.query,
                              route_path=args.route_request, previous_round=args.previous_round)
    elif args.command == "status":
        result = inspect_round(args.round)
    else:
        result = record_decision(args.round, args.candidate, args.action, args.reviewer, args.reason,
                                  acknowledge_limitations=args.acknowledge_limitations)
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return result


if __name__ == "__main__":
    main()
