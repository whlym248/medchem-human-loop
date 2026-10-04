"""Export graph proposals as missing-evidence records; no independent chemistry check."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

RULES = {"aromatic_ch_to_n", "aromatic_ch_add_cn", "aromatic_ch_add_f", "aromatic_ch_add_cl"}
NOTE = ("Graph proposals only, not predictions or reviewed scientific results. "
        "This exporter checks supplied hashes and provenance fields, not chemistry or generation authenticity.")


def _sha(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _smiles(value, digest):
    if (not isinstance(value, str) or not value or len(value) > 4096
            or any(char.isspace() or ord(char) < 32 for char in value)):
        raise ValueError("Expected a single plain SMILES string; chemistry is not verified")
    if _sha(value) != digest:
        raise ValueError("Structure SHA mismatch")
    return value


def _unique(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _invalid_number(value):
    raise ValueError(f"Invalid JSON number: {value}")


def export_proposals(proposals, targets, scope, out):
    """Create new .smi files and input.json accepted by assemble_evidence.py."""
    out = Path(out)
    if out.exists() or out.is_symlink():
        raise ValueError("Output must be a new directory; refusing to overwrite")
    if scope not in {"synthetic_example", "existing_evidence_only"}:
        raise ValueError("scope must be synthetic_example or existing_evidence_only")
    if (not isinstance(targets, (list, tuple)) or not 1 <= len(targets) <= 50
            or any(not isinstance(t, str) or not t.strip() or len(t) > 128
                   or any(ord(c) < 32 for c in t) for t in targets)
            or len(set(targets)) != len(targets)):
        raise ValueError("Provide 1 to 50 unique nonempty target names")
    with Path(proposals).open("rb") as handle:
        raw = handle.read(2 * 1024 * 1024 + 1)
    if len(raw) > 2 * 1024 * 1024:
        raise ValueError("Proposal report exceeds 2 MiB")
    report = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=_unique,
                        parse_constant=_invalid_number)
    if not isinstance(report, dict) or report.get("schema_version") != "local-edits/v1":
        raise ValueError("Expected local-edits/v1 proposal report")
    parent, request = report.get("parent"), report.get("request")
    if not isinstance(parent, dict) or not isinstance(request, dict):
        raise ValueError("Parent and request provenance are required")
    _smiles(parent.get("input_smiles"), parent.get("input_smiles_sha256"))
    canonical = _smiles(parent.get("canonical_smiles"), parent.get("canonical_smiles_sha256"))
    site, rules = request.get("site_index"), request.get("rules")
    if (type(site) is not int or not 0 <= site < 256 or not isinstance(rules, list)
            or not 1 <= len(rules) <= 64 or any(not isinstance(r, str) or r not in RULES for r in rules)):
        raise ValueError("Invalid requested site or rules")
    proposals_list = report.get("candidates")
    if not isinstance(proposals_list, list) or len(proposals_list) > 20:
        raise ValueError("Expected at most 20 candidate proposals")
    # Scope is provenance, not something a CLI flag may upgrade. Earlier
    # generator reports lack it, so their required CLI scope is only an
    # unauthenticated caller declaration (made explicit in the exported input).
    scope_sources = [report, parent, request, report.get("provenance", {}), *proposals_list]
    declared_scopes = [item["scope"] for item in scope_sources
                       if isinstance(item, dict) and "scope" in item]
    if any(value not in ("synthetic_example", "existing_evidence_only") for value in declared_scopes):
        raise ValueError("Unsupported source scope declaration")
    if any(value != scope for value in declared_scopes):
        raise ValueError("Cannot change declared source scope; synthetic evidence cannot be upgraded")
    structures = [("parent", canonical, "Supplied parent structure; " + NOTE)]
    ids, seen_smiles = {"parent"}, {canonical}
    for item in proposals_list:
        if not isinstance(item, dict):
            raise ValueError("Each proposal must be an object")
        name = item.get("candidate_id")
        if not isinstance(name, str) or not re.fullmatch(r"edit_[0-9a-f]{12}", name):
            raise ValueError("Unsafe or unsupported candidate ID")
        if name in ids:
            raise ValueError("Duplicate candidate ID")
        smiles = _smiles(item.get("smiles"), item.get("smiles_sha256"))
        if name != "edit_" + _sha(smiles)[:12] or smiles in seen_smiles:
            raise ValueError("Candidate identity mismatch or duplicate structure")
        origin = item.get("origin", {})
        if (not isinstance(origin, dict)
                or origin.get("parent_input_smiles_sha256") != parent["input_smiles_sha256"]
                or type(origin.get("site_atom_map")) is not int or origin["site_atom_map"] != site + 1
                or type(origin.get("transformation_count")) is not int or origin["transformation_count"] != 1):
            raise ValueError("Candidate parent origin or site provenance mismatch")
        if (item.get("rule") not in rules or item.get("rule") not in RULES
                or type(item.get("site_index")) is not int or item["site_index"] != site):
            raise ValueError("Candidate rule or site differs from request")
        ids.add(name)
        seen_smiles.add(smiles)
        structures.append((name, smiles, f"Supplied graph proposal from {item['rule']}; " + NOTE))
    payload = {"schema": "medchem_evidence_v1", "scope": scope, "source_note": NOTE,
               "targets": list(targets), "candidates": [], "evidence": [],
               "limitations": [NOTE, "All property, docking and MD data are missing, not zero.",
                               "No synthesis, affinity, safety or optimization benefit is established."],
               "proposal_provenance": {"generation_file_sha256": hashlib.sha256(raw).hexdigest(),
                                       "generation_schema": "local-edits/v1",
                                       "source_scope_declared": bool(declared_scopes),
                                       "scope_basis": "source_declaration" if declared_scopes else "caller_declaration_only",
                                       "scope_authenticity_verified": False,
                                       "independent_chemical_validation": False,
                                       "generation_authenticity_verified": False}}
    if not declared_scopes:
        payload["limitations"].append("Source report has no scope; the required CLI scope is an unauthenticated caller declaration.")
    files = {}
    for name, smiles, source in structures:
        filename, contents = name + ".smi", smiles + "\n"
        files[filename] = contents
        payload["candidates"].append({"id": name, "smiles": smiles, "source": source,
                                      "structure_file": filename, "structure_sha256": _sha(contents),
                                      "properties": {}})
        for target in targets:
            payload["evidence"].append({"candidate": name, "target": target, "source": NOTE,
                                        "docking_score_kcal_mol": None, "md": {"status": "not_available"}})
    files["input.json"] = json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    out.mkdir(parents=True, exist_ok=False)
    for name, contents in files.items():
        (out / name).write_text(contents, encoding="utf-8", newline="\n")
    return payload


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proposals", required=True, type=Path)
    parser.add_argument("--target", required=True, action="append")
    parser.add_argument("--scope", required=True, choices=["synthetic_example", "existing_evidence_only"])
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    result = export_proposals(args.proposals, args.target, args.scope, args.out)
    print(json.dumps({"input": str(args.out / "input.json"), "parent_id": "parent",
                      "candidate_target_rows": len(result["evidence"])}))
    return result


if __name__ == "__main__":
    main()
