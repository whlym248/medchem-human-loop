"""Bounded, human-selected aromatic C-H edits; no docking, MD or model inference.

RDKit is optional for the rest of this repository. Importing this module does not
require it; calling enumerate_edits does. Atom indices refer to the parsed input
SMILES before canonicalization, and atom-map numbers equal those indices + 1.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem, rdBase
except ImportError:
    Chem = None
    rdBase = None


RULES = (
    "aromatic_ch_to_n",
    "aromatic_ch_add_cn",
    "aromatic_ch_add_f",
    "aromatic_ch_add_cl",
)
SCHEMA_VERSION = "local-edits/v1"
MAX_ATOMS = 256


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _index(value: Any, count: int, label: str) -> int:
    if type(value) is not int or not 0 <= value < count:
        raise ValueError(f"{label} must be an integer in [0, {count - 1}]")
    return value


def _atom_label(atom: Any) -> tuple:
    return (
        atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetIsotope(),
        atom.GetIsAromatic(), str(atom.GetChiralTag()),
        atom.GetNumRadicalElectrons(), atom.GetAtomMapNum(),
    )


def _bond_label(bond: Any) -> tuple:
    return (str(bond.GetBondType()), bond.GetIsAromatic(), str(bond.GetStereo()))


def _edges(mol: Any, limit: int) -> dict:
    return {
        tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))): _bond_label(bond)
        for bond in mol.GetBonds()
        if bond.GetBeginAtomIdx() < limit and bond.GetEndAtomIdx() < limit
    }


def _fixed_neighborhood(mol: Any, idx: int) -> tuple:
    atom = mol.GetAtomWithIdx(idx)
    neighbors = sorted(
        (neighbor.GetIdx(), _atom_label(neighbor),
         _bond_label(mol.GetBondBetweenAtoms(idx, neighbor.GetIdx())))
        for neighbor in atom.GetNeighbors()
    )
    return (_atom_label(atom), atom.GetTotalNumHs(), tuple(neighbors))


def _canonical(mol: Any) -> str:
    clean = Chem.Mol(mol)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(clean, canonical=True, isomericSmiles=True)


def _parse(request: dict) -> tuple:
    if Chem is None:
        raise RuntimeError("RDKit is required for local edits; see requirements-chem.txt")
    if not isinstance(request, dict):
        raise ValueError("request must be a JSON object")
    allowed = {"smiles", "site_index", "rules", "fixed_atom_indices", "max_candidates", "scope"}
    unknown = set(request) - allowed
    if unknown:
        raise ValueError(f"Unsupported request fields: {sorted(unknown)}")
    if "scope" in request and (not isinstance(request["scope"], str)
                               or request["scope"] not in ("synthetic_example", "existing_evidence_only")):
        raise ValueError("scope must be synthetic_example or existing_evidence_only when supplied")
    smiles = request.get("smiles")
    if not isinstance(smiles, str) or not smiles or len(smiles) > 4096:
        raise ValueError("smiles must be a nonempty string of at most 4096 characters")
    if any(char.isspace() for char in smiles) or "|" in smiles:
        raise ValueError("Use plain SMILES only, without names, whitespace or CXSMILES")
    params = Chem.SmilesParserParams()
    params.removeHs = False
    params.parseName = False
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        raise ValueError("SMILES could not be parsed and sanitized by RDKit")
    if not 1 <= mol.GetNumAtoms() <= MAX_ATOMS:
        raise ValueError(f"Input must have 1 to {MAX_ATOMS} atoms")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Multiple fragments and salts are unsupported")
    if any(atom.GetFormalCharge() for atom in mol.GetAtoms()):
        raise ValueError("Charged atoms are unsupported, including net-neutral zwitterions")
    if any(atom.GetAtomMapNum() for atom in mol.GetAtoms()):
        raise ValueError("Input atom maps are unsupported; this tool assigns index + 1")
    if any(atom.GetIsotope() or atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        raise ValueError("Isotopes and radicals are unsupported")
    if any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        raise ValueError("Explicit hydrogen atoms are unsupported; use implicit hydrogens")
    if any(atom.GetAtomicNum() not in {6, 7, 8, 9, 15, 16, 17, 35, 53}
           for atom in mol.GetAtoms()):
        raise ValueError("Only neutral C/N/O/F/P/S/Cl/Br/I inputs are supported")
    if (any(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            for atom in mol.GetAtoms()) or Chem.FindPotentialStereo(mol)):
        raise ValueError("Specified or potential stereochemistry is unsupported in this version")
    site = _index(request.get("site_index"), mol.GetNumAtoms(), "site_index")
    rules = request.get("rules")
    if (not isinstance(rules, list) or not 1 <= len(rules) <= 64
            or any(not isinstance(rule, str) or rule not in RULES for rule in rules)):
        raise ValueError(f"rules must contain 1 to 64 names from {list(RULES)}")
    fixed = request.get("fixed_atom_indices", [])
    if not isinstance(fixed, list):
        raise ValueError("fixed_atom_indices must be a list")
    fixed = [_index(idx, mol.GetNumAtoms(), "fixed_atom_indices item") for idx in fixed]
    if len(set(fixed)) != len(fixed):
        raise ValueError("fixed_atom_indices must not contain duplicates")
    maximum = request.get("max_candidates", 4)
    if type(maximum) is not int or not 1 <= maximum <= 20:
        raise ValueError("max_candidates must be an integer from 1 to 20")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    return mol, site, rules, sorted(fixed), maximum


def _edit(parent: Any, site: int, rule: str) -> Any:
    editable = Chem.RWMol(parent)
    atom = editable.GetAtomWithIdx(site)
    atom.SetNumExplicitHs(0)
    if rule == "aromatic_ch_to_n":
        atom.SetAtomicNum(7)
        atom.SetNoImplicit(True)
    else:
        atom.SetNoImplicit(False)
        element = {"aromatic_ch_add_cn": 6, "aromatic_ch_add_f": 9,
                   "aromatic_ch_add_cl": 17}[rule]
        new_idx = editable.AddAtom(Chem.Atom(element))
        editable.GetAtomWithIdx(new_idx).SetAtomMapNum(new_idx + 1)
        editable.AddBond(site, new_idx, Chem.BondType.SINGLE)
        if rule == "aromatic_ch_add_cn":
            nitrogen = editable.AddAtom(Chem.Atom(7))
            editable.GetAtomWithIdx(nitrogen).SetAtomMapNum(nitrogen + 1)
            editable.AddBond(new_idx, nitrogen, Chem.BondType.TRIPLE)
    candidate = editable.GetMol()
    Chem.SanitizeMol(candidate)
    Chem.AssignStereochemistry(candidate, cleanIt=True, force=True)
    return candidate


def _audit(parent: Any, candidate: Any, site: int, fixed: list[int]) -> dict:
    count = parent.GetNumAtoms()
    if candidate.GetNumAtoms() < count or _edges(parent, count) != _edges(candidate, count):
        raise ValueError("An existing scaffold edge changed")
    for idx in range(count):
        before, after = parent.GetAtomWithIdx(idx), candidate.GetAtomWithIdx(idx)
        if after.GetAtomMapNum() != before.GetAtomMapNum():
            raise ValueError(f"Parent atom map {idx + 1} changed")
        if idx != site and (_atom_label(before) != _atom_label(after)
                            or before.GetTotalNumHs() != after.GetTotalNumHs()):
            raise ValueError(f"Unselected scaffold atom {idx} changed")
    for idx in fixed:
        if _fixed_neighborhood(parent, idx) != _fixed_neighborhood(candidate, idx):
            raise ValueError(f"Protected atom {idx} or its immediate neighborhood changed")
    if (len(Chem.GetMolFrags(candidate)) != 1
            or any(atom.GetFormalCharge() or atom.GetNumRadicalElectrons()
                   for atom in candidate.GetAtoms())):
        raise ValueError("Candidate is fragmented, charged or radical")
    if Chem.FindPotentialStereo(candidate):
        raise ValueError("Candidate introduces unsupported potential stereochemistry")
    for bond in candidate.GetBonds():
        ends = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if min(ends) < count <= max(ends) and min(ends) != site:
            raise ValueError("Added fragment attaches outside selected site")
    if candidate.GetAtomWithIdx(site).GetTotalNumHs() != 0:
        raise ValueError("Selected C-H hydrogen was not replaced")
    if not candidate.GetAtomWithIdx(site).GetIsAromatic():
        raise ValueError("Selected site lost aromaticity")
    return {
        "edited_parent_atom_indices": [site],
        "retained_parent_atom_maps": list(range(1, count + 1)),
        "added_atom_maps": list(range(count + 1, candidate.GetNumAtoms() + 1)),
        "existing_scaffold_edges_unchanged": True,
        "unselected_parent_atom_labels_unchanged": True,
        "fixed_neighborhoods_unchanged": True,
        "sanitized_by_rdkit": True,
        "net_formal_charge": 0,
    }


def enumerate_edits(request: dict) -> dict:
    """Return deterministic candidates/rejections; invalid requests raise ValueError.

    These are one-site graph proposals only. No biological ranking is assigned.
    The caller chooses the site, rules and fixed atoms explicitly.
    """
    parent, site, rules, fixed, maximum = _parse(request)
    original = request["smiles"]
    canonical = _canonical(parent)
    parent_sha = _sha(original)
    output = {
        "schema_version": SCHEMA_VERSION,
        "parent": {
            "input_smiles": original,
            "input_smiles_sha256": parent_sha,
            "canonical_smiles": canonical,
            "canonical_smiles_sha256": _sha(canonical),
            "mapped_smiles": Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True),
            "atoms": [{"index": atom.GetIdx(), "atom_map": atom.GetAtomMapNum(),
                       "element": atom.GetSymbol(), "aromatic": atom.GetIsAromatic(),
                       "hydrogen_count": atom.GetTotalNumHs()}
                      for atom in parent.GetAtoms()],
        },
        "request": {"site_index": site, "rules": list(rules),
                    "fixed_atom_indices": fixed, "max_candidates": maximum},
        "candidates": [], "rejected": [],
        "provenance": {
            "generator": "enumerate_local_edits.py",
            "method": "finite human-selected aromatic C-H graph transformations",
            "rdkit_version": rdBase.rdkitVersion,
            "new_MD": False, "new_docking": False, "new_AI_inference": False,
            "generation_ranking": "none; requested rule order only",
        },
        "warnings": [
            "Graph validity does not establish synthesizability, stability, safety or efficacy.",
            "These candidates are hypotheses, not recommendations or affinity predictions.",
            "No 3D coordinates, protonation selection, docking, ADMET or MD were computed.",
            "Input atom indices are RDKit input order, not canonical SMILES order.",
        ],
    }
    if "scope" in request:
        output["scope"] = request["scope"]
        output["request"]["scope"] = request["scope"]
    atom = parent.GetAtomWithIdx(site)
    site_problem = None
    if site in fixed:
        site_problem = "Selected site is protected by fixed_atom_indices"
    elif not (atom.GetAtomicNum() == 6 and atom.GetIsAromatic()
              and atom.GetTotalNumHs() == 1 and atom.GetDegree() == 2):
        site_problem = "Selected site must be an aromatic carbon with one H and two neighbors"
    seen = {canonical}
    seen_rules = set()
    for rule in rules:
        reason = None
        if rule in seen_rules:
            reason = "Duplicate rule; candidate already considered"
        elif site_problem:
            reason = site_problem
        elif len(output["candidates"]) >= maximum:
            reason = "Candidate limit reached"
        seen_rules.add(rule)
        if reason is None:
            try:
                candidate = _edit(parent, site, rule)
                audit = _audit(parent, candidate, site, fixed)
                identity = _canonical(candidate)
                if identity in seen:
                    reason = "Duplicate canonical structure or unchanged parent"
                else:
                    seen.add(identity)
                    output["candidates"].append({
                        "candidate_id": "edit_" + _sha(identity)[:12],
                        "rule": rule, "site_index": site,
                        "smiles": identity,
                        "smiles_sha256": _sha(identity),
                        "mapped_smiles": Chem.MolToSmiles(candidate, canonical=True,
                                                              isomericSmiles=True),
                        "origin": {"parent_input_smiles_sha256": parent_sha,
                                   "site_atom_map": site + 1,
                                   "transformation_count": 1},
                        "audit": audit,
                    })
            except (ValueError, RuntimeError) as error:
                reason = f"Chemical/constraint validation failed: {error}"
        if reason is not None:
            output["rejected"].append({"rule": rule, "site_index": site, "reason": reason})
    return output


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="JSON edit request")
    parser.add_argument("--out", required=True, type=Path, help="New, nonexistent output directory")
    args = parser.parse_args(argv)
    if args.out.exists():
        raise ValueError("Output must be a new, nonexistent directory; refusing to overwrite")
    request = json.loads(args.input.read_text(encoding="utf-8-sig"))
    result = enumerate_edits(request)
    report = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    candidates = "".join(f"{item['smiles']}\t{item['candidate_id']}\n"
                         for item in result["candidates"])
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "local_edits.json").write_text(report, encoding="utf-8")
    (args.out / "candidates.smi").write_text(candidates, encoding="utf-8")
    print(json.dumps({"candidate_count": len(result["candidates"]),
                      "rejected_count": len(result["rejected"]),
                      "report": str(args.out / "local_edits.json")}, ensure_ascii=False))
    return result


if __name__ == "__main__":
    main()
