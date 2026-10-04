"""Small graph-only tests. No coordinates, GPU work or MD execution."""
import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import enumerate_local_edits as edits


@unittest.skipIf(edits.Chem is None, "RDKit unavailable; optional graph tests skipped")
class LocalEditTests(unittest.TestCase):
    def setUp(self):
        self.request = json.loads((ROOT / "examples" / "local_edit_request.json")
                                  .read_text(encoding="utf-8"))

    def test_four_graphs_are_valid_and_preserve_parent_maps(self):
        report = edits.enumerate_edits(self.request)
        self.assertEqual(len(report["candidates"]), 4)
        self.assertEqual(report["rejected"], [])
        parent_count = len(report["parent"]["atoms"])
        self.assertEqual(parent_count, 7)
        for item in report["candidates"]:
            molecule = edits.Chem.MolFromSmiles(item["mapped_smiles"])
            self.assertIsNotNone(molecule)
            maps = {atom.GetAtomMapNum(): atom for atom in molecule.GetAtoms()}
            self.assertEqual(sorted(maps), list(range(1, molecule.GetNumAtoms() + 1)))
            self.assertTrue(all(maps[idx + 1].GetSymbol() == "C"
                                for idx in range(parent_count) if idx != 3))
            self.assertEqual(item["audit"]["edited_parent_atom_indices"], [3])
            self.assertEqual(item["origin"]["parent_input_smiles_sha256"],
                             report["parent"]["input_smiles_sha256"])
            self.assertTrue(item["audit"]["existing_scaffold_edges_unchanged"])
            self.assertEqual(sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()), 0)
            site = maps[4]
            self.assertEqual(site.GetTotalNumHs(), 0)
            self.assertTrue(site.GetIsAromatic())
        identities = {item["smiles"] for item in report["candidates"]}
        self.assertEqual(len(identities), 4)
        nitrogen_edit = report["candidates"][0]
        self.assertEqual(nitrogen_edit["audit"]["added_atom_maps"], [])
        self.assertEqual(report["candidates"][1]["audit"]["added_atom_maps"], [8, 9])

    def test_output_is_deterministic_and_input_unchanged(self):
        original = copy.deepcopy(self.request)
        first = edits.enumerate_edits(self.request)
        second = edits.enumerate_edits(self.request)
        self.assertEqual(first, second)
        self.assertEqual(self.request, original)

    def test_scope_propagates_and_legacy_scope_stays_undeclared(self):
        report = edits.enumerate_edits(self.request)
        self.assertEqual(report["scope"], "synthetic_example")
        self.assertEqual(report["request"]["scope"], "synthetic_example")
        self.request["scope"] = "existing_evidence_only"
        self.assertEqual(edits.enumerate_edits(self.request)["scope"], "existing_evidence_only")
        self.request.pop("scope")
        legacy = edits.enumerate_edits(self.request)
        self.assertNotIn("scope", legacy)
        self.assertNotIn("scope", legacy["request"])

    def test_invalid_scope_is_rejected(self):
        for scope in [None, True, 1, [], {}, "validated", "", "synthetic_example "]:
            with self.subTest(scope=scope):
                self.request["scope"] = scope
                with self.assertRaisesRegex(ValueError, "scope must be"):
                    edits.enumerate_edits(self.request)

    def test_fixed_site_rejects_addition_and_substitution(self):
        self.request["fixed_atom_indices"] = [3]
        result = edits.enumerate_edits(self.request)
        self.assertEqual(result["candidates"], [])
        self.assertEqual(len(result["rejected"]), 4)
        self.assertTrue(all("protected" in row["reason"] for row in result["rejected"]))

    def test_fixed_neighbor_element_changes_are_rejected(self):
        self.request["site_index"] = 2
        self.request["fixed_atom_indices"] = [1]
        result = edits.enumerate_edits(self.request)
        self.assertEqual(result["rejected"][0]["rule"], "aromatic_ch_to_n")
        self.assertIn("neighborhood", result["rejected"][0]["reason"])
        self.assertEqual(len(result["candidates"]), 3)

    def test_duplicate_rules_and_limit_are_reported(self):
        self.request["rules"] = ["aromatic_ch_add_f", "aromatic_ch_add_f",
                                 "aromatic_ch_add_cn"]
        self.request["max_candidates"] = 1
        result = edits.enumerate_edits(self.request)
        self.assertEqual(len(result["candidates"]), 1)
        self.assertIn("Duplicate", result["rejected"][0]["reason"])
        self.assertIn("limit", result["rejected"][1]["reason"])

    def test_duplicate_canonical_structure_is_not_emitted(self):
        self.request["rules"] = ["aromatic_ch_add_f", "aromatic_ch_add_cl"]
        original_edit = edits._edit

        def same_product(parent, site, rule):
            return original_edit(parent, site, "aromatic_ch_add_f")

        with patch.object(edits, "_edit", side_effect=same_product):
            result = edits.enumerate_edits(self.request)
        self.assertEqual(len(result["candidates"]), 1)
        self.assertIn("Duplicate canonical", result["rejected"][0]["reason"])

    def test_non_ch_aromatic_and_aliphatic_sites_are_rejected(self):
        for site in [0, 1]:
            with self.subTest(site=site):
                self.request["site_index"] = site
                self.request["fixed_atom_indices"] = []
                result = edits.enumerate_edits(self.request)
                self.assertEqual(result["candidates"], [])
                self.assertIn("aromatic carbon with one H", result["rejected"][0]["reason"])

    def test_bracket_ch_hydrogen_is_replaced(self):
        self.request.update(smiles="[cH]1ccccc1", site_index=0, fixed_atom_indices=[])
        result = edits.enumerate_edits(self.request)
        self.assertEqual(len(result["candidates"]), 4)

    def test_symmetry_breaking_that_creates_stereocenter_is_rejected(self):
        # The parent has two equivalent phenyl groups; editing only one makes
        # the central carbon stereogenic even though its bonds are unchanged.
        self.request.update(smiles="CC(c1ccccc1)c1ccccc1", site_index=3,
                            fixed_atom_indices=[])
        result = edits.enumerate_edits(self.request)
        self.assertEqual(result["candidates"], [])
        self.assertEqual(len(result["rejected"]), 4)
        self.assertTrue(all("stereochemistry" in item["reason"]
                            for item in result["rejected"]))

    def test_input_rejects_charge_fragments_stereo_and_radicals(self):
        bad_inputs = {
            "C[NH2+]Cc1ccccc1": "Charged",
            "[O-]C(=O)C[NH3+]": "Charged",
            "Cc1ccccc1.O": "Multiple fragments",
            "C[C@H](F)c1ccccc1": "stereochemistry",
            "CC(F)c1ccccc1": "stereochemistry",
            "C/C=C/c1ccccc1": "stereochemistry",
            "[CH2]c1ccccc1": "radicals",
            "[13CH3]c1ccccc1": "Isotopes",
            "[CH3:1]c1ccccc1": "atom maps",
            "[H]c1ccccc1": "Explicit hydrogen",
            "CC invalid": "plain SMILES",
            "C(C)(C)(C)(C)C": "parsed and sanitized",
        }
        for smiles, error in bad_inputs.items():
            with self.subTest(smiles=smiles):
                self.request.update(smiles=smiles, site_index=0, fixed_atom_indices=[])
                with self.assertRaisesRegex(ValueError, error):
                    edits.enumerate_edits(self.request)

    def test_invalid_indices_rules_and_bounds_rejected(self):
        invalid_fields = [
            {"site_index": -1}, {"site_index": 7}, {"site_index": True},
            {"fixed_atom_indices": [7]}, {"fixed_atom_indices": [0, 0]},
            {"max_candidates": 0}, {"max_candidates": 21}, {"max_candidates": True},
            {"rules": []}, {"rules": ["[c:1]>>[n:1]"]},
            {"rules": ["aromatic_ch_add_f"] * 65}, {"script": "print('no')"},
        ]
        for change in invalid_fields:
            with self.subTest(change=change):
                request = dict(self.request, **change)
                with self.assertRaises(ValueError):
                    edits.enumerate_edits(request)

    def test_new_directory_only_and_saved_report_matches(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "request.json"
            source.write_text(json.dumps(self.request), encoding="utf-8")
            output = Path(temporary) / "new"
            with contextlib.redirect_stdout(io.StringIO()):
                report = edits.main(["--input", str(source), "--out", str(output)])
            self.assertEqual(json.loads((output / "local_edits.json").read_text()), report)
            self.assertEqual(len((output / "candidates.smi").read_text().splitlines()), 4)
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                edits.main(["--input", str(source), "--out", str(output)])

    def test_existing_bond_change_fails_audit(self):
        parent, site, _, fixed, _ = edits._parse(self.request)
        altered = edits.Chem.RWMol(parent)
        altered.RemoveBond(0, 1)
        with self.assertRaisesRegex(ValueError, "scaffold edge"):
            edits._audit(parent, altered.GetMol(), site, fixed)


class OptionalDependencyTests(unittest.TestCase):
    def test_missing_rdkit_has_actionable_error(self):
        with patch.object(edits, "Chem", None):
            with self.assertRaisesRegex(RuntimeError, "requirements-chem.txt"):
                edits.enumerate_edits({})


if __name__ == "__main__":
    unittest.main()
