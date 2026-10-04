"""Standard-library hand-written proposal fixtures; no chemical or MD calculations."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import assemble_evidence
import proposals_to_evidence as bridge
import workflow_round


def sha(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


class ProposalBridgeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source, self.out = self.root / "proposals.json", self.root / "new"
        parent, candidate = "c1ccccc1", "c1ccncc1"
        self.fixture = {
            "schema_version": "local-edits/v1",
            "parent": {"input_smiles": parent, "input_smiles_sha256": sha(parent),
                       "canonical_smiles": parent, "canonical_smiles_sha256": sha(parent)},
            "request": {"site_index": 3, "rules": ["aromatic_ch_to_n"]},
            "candidates": [{"candidate_id": "edit_" + sha(candidate)[:12],
                            "smiles": candidate, "smiles_sha256": sha(candidate),
                            "site_index": 3, "rule": "aromatic_ch_to_n",
                            "origin": {"parent_input_smiles_sha256": sha(parent),
                                       "site_atom_map": 4, "transformation_count": 1}}],
        }

    def write(self):
        self.source.write_text(json.dumps(self.fixture), encoding="utf-8")

    def export(self):
        self.write()
        return bridge.export_proposals(self.source, ["toy_A", "toy_B"], "synthetic_example", self.out)

    def test_export_assembles_without_inventing_missing_evidence(self):
        result = self.export()
        parsed, rows = assemble_evidence.assemble(self.out / "input.json")
        self.assertEqual(result, parsed)
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(candidate["properties"] == {} for candidate in result["candidates"]))
        for row in rows:
            self.assertIsNone(row["supplied_docking_score_kcal_mol"])
            self.assertIsNone(row["supplied_ADMET_raw_score"])
            self.assertIsNone(row["MD_sampling_ns"])
            self.assertEqual(row["MD_evidence_status"], "not_available")
            self.assertEqual(row["human_decision"], "not_recorded")
        provenance = result["proposal_provenance"]
        self.assertEqual(provenance["generation_file_sha256"],
                         hashlib.sha256(self.source.read_bytes()).hexdigest())
        self.assertFalse(provenance["independent_chemical_validation"])

    def test_hash_tampering_fails_before_output_creation(self):
        self.fixture["candidates"][0]["smiles"] = "CCO"
        with self.assertRaisesRegex(ValueError, "SHA mismatch"):
            self.export()
        self.assertFalse(self.out.exists())

    def test_duplicate_ids_fail(self):
        self.fixture["candidates"].append(dict(self.fixture["candidates"][0]))
        with self.assertRaisesRegex(ValueError, "Duplicate candidate ID"):
            self.export()

    def test_unsafe_path_ids_fail(self):
        for name in ["../escape", "..\\escape", "D:/escape", "/absolute", "parent", "CON"]:
            with self.subTest(name=name):
                self.fixture["candidates"][0]["candidate_id"] = name
                with self.assertRaisesRegex(ValueError, "candidate ID"):
                    self.export()

    def test_existing_output_is_never_overwritten(self):
        self.out.mkdir()
        with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
            self.export()

    def test_origin_and_rule_mismatches_fail(self):
        candidate = self.fixture["candidates"][0]
        candidate["origin"]["parent_input_smiles_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "origin"):
            self.export()
        candidate["origin"]["parent_input_smiles_sha256"] = self.fixture["parent"]["input_smiles_sha256"]
        candidate["rule"] = "aromatic_ch_add_f"
        with self.assertRaisesRegex(ValueError, "rule"):
            self.export()

    def test_supplied_scores_are_not_carried_into_evidence(self):
        candidate = self.fixture["candidates"][0]
        candidate.update(properties={"ADMET_raw_score": 0.1}, docking_score_kcal_mol=-100,
                         md={"status": "completed_reviewed"})
        result = self.export()
        self.assertTrue(all(item["md"]["status"] == "not_available" for item in result["evidence"]))
        self.assertTrue(all(item["docking_score_kcal_mol"] is None for item in result["evidence"]))
        self.assertTrue(all(item["properties"] == {} for item in result["candidates"]))

    def test_scope_targets_and_parent_hash_validated(self):
        self.write()
        for targets, scope in [(["a", "a"], "synthetic_example"), (["a"], "validated")]:
            with self.subTest(targets=targets, scope=scope):
                with self.assertRaises(ValueError):
                    bridge.export_proposals(self.source, targets, scope, self.out)
        self.fixture["parent"]["canonical_smiles_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "SHA mismatch"):
            self.export()

    def test_cli_creates_input_json_and_parent_smi(self):
        self.write()
        with contextlib.redirect_stdout(io.StringIO()) as output:
            bridge.main(["--proposals", str(self.source), "--target", "toy_A",
                         "--scope", "existing_evidence_only", "--out", str(self.out)])
        receipt = json.loads(output.getvalue())
        self.assertEqual(receipt["parent_id"], "parent")
        self.assertTrue((self.out / "input.json").is_file())
        self.assertTrue((self.out / "parent.smi").is_file())
        self.assertEqual(assemble_evidence.assemble(self.out / "input.json")[0]["scope"],
                         "existing_evidence_only")

    def test_declared_synthetic_scope_cannot_be_upgraded(self):
        for location in [self.fixture, self.fixture["parent"], self.fixture["candidates"][0]]:
            with self.subTest(location=list(location)):
                location["scope"] = "synthetic_example"
                self.write()
                with self.assertRaisesRegex(ValueError, "Cannot change declared source scope"):
                    bridge.export_proposals(self.source, ["toy"], "existing_evidence_only", self.out)
                location.pop("scope")
        self.assertFalse(self.out.exists())

    def test_missing_source_scope_is_explicitly_unauthenticated(self):
        result = self.export()
        provenance = result["proposal_provenance"]
        self.assertFalse(provenance["source_scope_declared"])
        self.assertEqual(provenance["scope_basis"], "caller_declaration_only")
        self.assertFalse(provenance["scope_authenticity_verified"])

    def test_round_keeps_scope_missing_metrics_and_integrity(self):
        self.fixture["scope"] = "synthetic_example"
        self.export()
        round_dir = self.root / "round"
        state = workflow_round.create_round(self.out / "input.json", round_dir, "parent")
        self.assertEqual(state["scope"], "synthetic_example")
        self.assertEqual(state["jobs_launched"], 0)
        packet = json.loads((round_dir / "review_packet.json").read_text(encoding="utf-8"))
        self.assertTrue(all(value is None for comparison in packet["comparisons"]
                            for value in comparison["delta_candidate_minus_parent"].values()))
        self.assertEqual(packet["parent_id"], "parent")
        self.assertEqual(set(packet["candidate_ids"]), {"parent", self.fixture["candidates"][0]["candidate_id"]})
        frozen_input = json.loads((round_dir / "evidence" / "input.json").read_text(encoding="utf-8"))
        self.assertTrue(frozen_input["proposal_provenance"]["source_scope_declared"])
        (round_dir / "evidence" / "parent.smi").write_text("CCO\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Frozen evidence changed"):
            workflow_round.inspect_round(round_dir)


if __name__ == "__main__":
    unittest.main()
