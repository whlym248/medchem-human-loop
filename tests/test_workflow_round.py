"""No MD or predictions: exercise local review packets and decision provenance."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import workflow_round as workflow


class WorkflowRoundTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.inputs = self.work / "inputs"
        self.inputs.mkdir()
        for name in ("synthetic_evidence.json", "toy_parent.smi", "toy_variant.smi"):
            shutil.copyfile(ROOT / "examples" / name, self.inputs / name)
        self.input = self.inputs / "synthetic_evidence.json"
        self.round = self.work / "round1"
        workflow.create_round(self.input, self.round, "toy_parent")

    def decide(self, candidate="toy_variant", action="approve_next_round"):
        return workflow.record_decision(self.round, candidate, action, "test reviewer", "Synthetic test only", acknowledge_limitations=True)

    def test_snapshot_is_independent_and_preserves_missing_data(self):
        self.input.write_text("modified after snapshot", encoding="utf-8")
        state = workflow.inspect_round(self.round)
        self.assertEqual(state["jobs_launched"], 0)
        self.assertEqual(state["events"], 0)
        packet = workflow.read(self.round / "review_packet.json")
        comparison = next(row for row in packet["comparisons"] if row["target"] == "toy_target_B")
        self.assertIsNone(comparison["delta_candidate_minus_parent"]["supplied_docking_score_kcal_mol"])
        self.assertFalse(comparison["comparable_protocols_verified"])

    def test_modified_snapshot_blocks_decision(self):
        (self.round / "evidence" / "toy_parent.smi").write_text("changed", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Frozen evidence changed"):
            self.decide()

    def test_required_acknowledgement_reviewer_reason_and_candidate(self):
        with self.assertRaisesRegex(ValueError, "acknowledgement"):
            workflow.record_decision(self.round, "toy_variant", "reject", "Reviewer", "Reason")
        with self.assertRaisesRegex(ValueError, "reviewer"):
            workflow.record_decision(self.round, "toy_variant", "reject", "", "Reason", acknowledge_limitations=True)
        with self.assertRaisesRegex(ValueError, "Unknown candidate"):
            self.decide(candidate="not in evidence")
        with self.assertRaisesRegex(ValueError, "acknowledgement"):
            workflow.record_decision(self.round, "toy_variant", "reject", "Reviewer", "Reason", acknowledge_limitations="false")

    def test_revised_decision_appends_chain_and_detects_prior_change(self):
        first = self.decide(action="defer")
        second = self.decide()
        self.assertEqual(first["sequence"], 1)
        self.assertEqual(second["sequence"], 2)
        self.assertEqual(workflow.inspect_round(self.round)["events"], 2)
        path = self.round / "decisions" / "decision_0001.json"
        item = workflow.read(path)
        item["reason"] = "changed after later event"
        path.write_text(json.dumps(item), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "chain mismatch"):
            workflow.inspect_round(self.round)

    def test_cannot_overwrite_or_continue_without_human_approval(self):
        with self.assertRaisesRegex(ValueError, "already exist"):
            workflow.create_round(self.input, self.round, "toy_parent")
        with self.assertRaisesRegex(ValueError, "no human"):
            workflow.create_round(self.input, self.work / "next", "toy_variant", previous_round=self.round)
        self.decide()
        state = workflow.create_round(self.input, self.work / "next", "toy_variant", previous_round=self.round)
        self.assertEqual(state["events"], 0)
        packet = workflow.read(self.work / "next" / "review_packet.json")
        self.assertEqual(packet["previous_round"]["parent_id"], "toy_variant")

    def test_approved_parent_cannot_change_identity(self):
        self.decide()
        data = workflow.read(self.input)
        data["candidates"][1]["smiles"] = "CCCCO"
        self.input.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "identity changed"):
            workflow.create_round(self.input, self.work / "next", "toy_variant", previous_round=self.round)

    def test_synthetic_cannot_become_real_by_renaming_scope(self):
        self.decide()
        data = workflow.read(self.input)
        data["scope"] = "existing_evidence_only"
        self.input.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "scope"):
                workflow.create_round(self.input, self.work / "next", "toy_variant", previous_round=self.round)

    def test_snapshot_rejects_parent_traversal_even_when_source_resolves_inside(self):
        data = workflow.read(self.input)
        data["candidates"][0]["structure_file"] = "../inputs/toy_parent.smi"
        self.input.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "parent traversal"):
            workflow.create_round(self.input, self.work / "escape", "toy_parent")
        self.assertFalse((self.work / "escape").exists())

    def test_different_admet_endpoints_are_not_subtracted(self):
        data = workflow.read(self.input)
        data["candidates"][1]["properties"]["ADMET_endpoint"] = "different_endpoint"
        self.input.write_text(json.dumps(data), encoding="utf-8")
        directory = self.work / "different_admet"
        workflow.create_round(self.input, directory, "toy_parent")
        packet = workflow.read(directory / "review_packet.json")
        self.assertTrue(all(row["delta_candidate_minus_parent"]["supplied_ADMET_raw_score"] is None for row in packet["comparisons"]))

    def test_report_escapes_supplied_html(self):
        packet = workflow.read(self.round / "review_packet.json")
        packet["limitations"].append('<script>alert("bad")</script>')
        rendered = workflow.render(packet)
        self.assertNotIn('<script>', rendered)
        self.assertIn('&lt;script&gt;', rendered)
        self.assertIn('ADMET_endpoint', rendered)
        self.assertIn('evidence/input.json', rendered)

    def test_synthetic_case_match_never_becomes_real_support(self):
        directory = self.work / "with_context"
        workflow.create_round(self.input, directory, "toy_parent",
                              case_path=ROOT / "examples" / "cases" / "synthetic_cases.jsonl",
                              query="aromatic", route_path=ROOT / "examples" / "route_request.json")
        packet = workflow.read(directory / "review_packet.json")
        route = packet["planning_context"]["route_suggestion"]
        self.assertFalse(route["request"]["case_support"])
        self.assertEqual(route["recommendation"], "collect_evidence")

    def test_changed_library_bytes_rejected_even_if_matches_unchanged(self):
        import case_library
        library = self.work / "cases.jsonl"
        shutil.copyfile(ROOT / "examples" / "cases" / "synthetic_cases.jsonl", library)
        original_search = case_library.search_cases
        def mutate_after_read(path, *args, **kwargs):
            result = original_search(path, *args, **kwargs)
            with library.open("a", encoding="utf-8") as stream:
                stream.write("\n")
            return result
        with patch.object(case_library, "search_cases", side_effect=mutate_after_read):
            with self.assertRaisesRegex(ValueError, "bytes changed"):
                workflow.create_round(self.input, self.work / "changing", "toy_parent", case_path=library, query="aromatic")


if __name__ == "__main__":
    unittest.main()
