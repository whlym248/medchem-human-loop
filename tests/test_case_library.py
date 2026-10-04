"""Traceability and retrieval boundaries; no web, inference, chemistry engine or MD."""
import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import case_library as cases


class CaseLibraryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "cases.jsonl"
        self.records = cases.load_library(ROOT / "examples/cases/synthetic_cases.jsonl")
        self.write()

    def write(self):
        self.path.write_text("\n".join(json.dumps(r) for r in self.records), encoding="utf-8")

    def literature(self, record):
        record["scope"] = "literature_curated"
        # Deliberately fictitious in-memory metadata: verifies shape, not a real case.
        record["source"] = {"url": "https://example.org/not-a-real-paper",
                            "citation": "TEST ONLY: a fictional citation",
                            "locator": "TEST ONLY: Table 1, compounds 1 and 2",
                            "license_note": "Test fixture, no third-party source copied"}
        record["review"] = {"status": "reviewed", "reviewer": "Test reviewer",
                            "note": "Test-only attestation; do not publish as a real reviewed case"}

    def test_fixture_counts_and_real_library_empty(self):
        summary = cases.validate_library(self.path)
        self.assertEqual(summary["scope_counts"]["synthetic_example"], 3)
        self.assertEqual(summary["review_status_counts"]["reviewed"], 2)
        real = cases.validate_library(ROOT / "examples/cases/literature_curated.jsonl")
        self.assertEqual(real["record_count"], 0)
        self.assertFalse(summary["scientific_claims_verified"])

    def test_synthetic_isolation_by_default_and_no_false_support(self):
        self.assertEqual(cases.search_cases(self.path, tags=["small_local"])["matches"], [])
        found = cases.search_cases(self.path, tags=["small_local"], allow_synthetic=True)
        self.assertEqual([r["id"] for r in found["matches"]], ["TOY-001", "TOY-002"])
        self.assertFalse(found["real_case_support_available"])

    def test_exact_query_tags_and_target_all_required(self):
        found = cases.search_cases(self.path, "AZA polarity", tags=["SMALL_LOCAL"],
                                   target="toy_target_a", allow_synthetic=True)
        self.assertEqual([r["id"] for r in found["matches"]], ["TOY-002"])
        self.assertFalse(cases.search_cases(self.path, "polarity unknownword", allow_synthetic=True)["matches"])
        self.assertFalse(cases.search_cases(self.path, "polarity", target="TOY_TARGET_B", allow_synthetic=True)["matches"])
        self.assertFalse(cases.search_cases(self.path, "fixture author", allow_synthetic=True)["matches"])

    def test_no_hit_is_not_a_recommendation(self):
        result = cases.search_cases(self.path, "nonexistent", allow_synthetic=True)
        self.assertEqual(result["returned_count"], 0)
        self.assertFalse(result["real_case_support_available"])
        self.assertNotIn("recommendation", result)

    def test_incomplete_or_boolean_review_rejected(self):
        self.records[0]["review"]["reviewer"] = ""
        self.write()
        with self.assertRaisesRegex(ValueError, "review.reviewer"):
            cases.load_library(self.path)
        self.records[0]["review"]["reviewer"] = "Some asserted name"
        self.records[0]["review"]["status"] = True
        self.write()
        with self.assertRaisesRegex(ValueError, "review.status"):
            cases.load_library(self.path)

    def test_literature_provenance_and_no_authenticity_claim(self):
        self.literature(self.records[0])
        self.write()
        found = cases.search_cases(self.path, "fluorination")
        self.assertTrue(found["real_case_support_available"])
        self.assertFalse(found["review_attestations_authenticated"])
        self.assertFalse(found["scientific_claims_verified"])
        self.records[0]["source"]["locator"] = ""
        self.write()
        with self.assertRaisesRegex(ValueError, "source.locator"):
            cases.load_library(self.path)

    def test_draft_and_rejected_literature_not_retrieved(self):
        for record in self.records[:2]:
            self.literature(record)
        self.records[0]["review"]["status"] = "draft"
        self.records[1]["review"]["status"] = "rejected"
        self.write()
        self.assertFalse(cases.search_cases(self.path, tags=["small_local"])["matches"])

    def test_protocols_are_preserved_not_pooled(self):
        result = cases.search_cases(self.path, tags=["small_local"], allow_synthetic=True)
        self.assertNotEqual(result["matches"][0]["endpoints"][0]["protocol"],
                            result["matches"][1]["endpoints"][0]["protocol"])
        self.assertNotIn("mean_effect", result)

    def test_duplicate_ids_keys_and_nonfinite_values_rejected(self):
        original = copy.deepcopy(self.records)
        self.records.append(copy.deepcopy(self.records[0]))
        self.write()
        with self.assertRaisesRegex(ValueError, "Duplicate case id"):
            cases.load_library(self.path)
        self.records = original
        self.records[0]["endpoints"][0]["value"] = float("nan")
        self.write()
        with self.assertRaisesRegex(ValueError, "numeric constant"):
            cases.load_library(self.path)
        self.path.write_text('{"id":"a", "id":"b"}', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Duplicate JSON key"):
            cases.load_library(self.path)

    def test_invalid_limits_blank_query_and_truthy_flag_rejected(self):
        for kwargs in ({}, {"query": "!!!"}, {"query": "aza", "limit": True},
                       {"query": "aza", "allow_synthetic": "yes"}, {"query": "aza", "tags": "aza"}):
            with self.assertRaises(ValueError):
                cases.search_cases(self.path, **kwargs)

    def test_urls_not_executed_and_credentials_rejected(self):
        self.records[0]["source"]["url"] = "https://user:secret@example.org/paper"
        self.write()
        with self.assertRaisesRegex(ValueError, "without credentials"):
            cases.load_library(self.path)
        self.records[0]["source"]["url"] = "file:///private/path"
        self.write()
        with self.assertRaisesRegex(ValueError, "HTTP"):
            cases.load_library(self.path)

    def test_return_limits_and_cli(self):
        result = cases.search_cases(self.path, tags=["small_local"], allow_synthetic=True, limit=1)
        self.assertEqual(result["match_count"], 2)
        self.assertEqual(result["returned_count"], 1)
        self.assertTrue(result["truncated"])
        with contextlib.redirect_stdout(io.StringIO()) as output:
            cases.main(["validate", "--library", str(self.path)])
        self.assertTrue(json.loads(output.getvalue())["valid"])


if __name__ == "__main__":
    unittest.main()
