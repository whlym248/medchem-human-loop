"""Decision boundaries of the transparent route proposal; never runs a method."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import case_library
import route_proposal as routes


class RouteProposalTests(unittest.TestCase):
    def request(self, **overrides):
        return {"edit_scope": "small_local", "fixed_region_required": True,
                "case_support": True, **overrides}

    def test_supported_local_edit_requires_human_and_core_check(self):
        result = routes.recommend(self.request())
        self.assertEqual(result["recommendation"], "limited_rule_edits")
        self.assertTrue(result["human_approval_required"])
        self.assertTrue(result["not_learned_policy"])
        self.assertFalse(result["executes_generation"])
        self.assertFalse(result["executes_MD"])
        self.assertTrue(any("Protected atoms" in s for s in result["reasons"]))

    def test_missing_evidence_or_unclear_scope_collects_evidence(self):
        self.assertEqual(routes.recommend(self.request(case_support=False))["recommendation"], "collect_evidence")
        self.assertEqual(routes.recommend(self.request(edit_scope="unclear"))["recommendation"], "collect_evidence")

    def test_fragment_rebuild_never_claims_diffusion_available(self):
        for availability in (False, True):
            result = routes.recommend(self.request(edit_scope="fragment_rebuild", diffusion_backend_available=availability))
            self.assertEqual(result["recommendation"], "defer_diffusion")
            self.assertEqual(result["blocked_methods"][0]["status"], "not_implemented")
            self.assertTrue(any("ignored" in s for s in result["reasons"]))

    def test_invalid_fields_and_truthy_booleans_rejected(self):
        for request in ({}, self.request(edit_scope="whatever"), self.request(case_support=1),
                        self.request(fixed_region_required="true"), self.request(diffusion_backend_available="yes"),
                        self.request(reward=100), []):
            with self.assertRaises(ValueError):
                routes.recommend(request)

    def test_no_hits_and_synthetic_hits_cannot_enable_real_support(self):
        path = ROOT / "examples/cases/synthetic_cases.jsonl"
        for query in ("aza", "nonexistent"):
            result = case_library.search_cases(path, query, allow_synthetic=True)
            route = routes.recommend(self.request(case_support=result["real_case_support_available"]))
            self.assertEqual(route["recommendation"], "collect_evidence")

    def test_input_not_mutated_and_cli_has_source_hash(self):
        request = self.request()
        before = dict(request)
        result = routes.recommend(request)
        result["request"]["case_support"] = False
        self.assertEqual(request, before)
        with contextlib.redirect_stdout(io.StringIO()) as output:
            routes.main(["--request", str(ROOT / "examples/route_request.json")])
        printed = json.loads(output.getvalue())
        self.assertEqual(printed["recommendation"], "collect_evidence")
        self.assertEqual(len(printed["request_sha256"]), 64)

    def test_cli_duplicate_keys_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "request.json"
            path.write_text('{"case_support":false,"case_support":true}', encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    routes.main(["--request", str(path)])
            self.assertEqual(caught.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
