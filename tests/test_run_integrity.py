"""Test restart bookkeeping with fake bytes only; never call run() or OpenMM."""
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_md_v2 as runner


class RunIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.identity = {"settings": {"timestep_ps": 0.002}, "equilibration_steps": 100,
                         "target_absolute_step": 600, "case_id": "toy"}
        self.digest = runner.digest(self.identity)

    def segment(self, index, start, end):
        path = self.root / "attempts" / f"attempt_{index:04d}" / "segment.json"
        runner.atomic_json(path, {"attempt": f"attempt_{index:04d}",
                                 "run_identity_sha256": self.digest,
                                 "start_after_absolute_step": start,
                                 "last_output_step": end, "report_interval_steps": 100})

    def test_superseded_tail_excluded(self):
        self.segment(1, 100, 600)
        self.segment(2, 300, 500)
        manifest = runner.rebuild_segment_manifest(self.root, self.identity)
        self.assertEqual([s["include_frame_count"] for s in manifest["segments"]], [2, 2])
        self.assertEqual(manifest["segments"][0]["superseded_tail_after_step"], 300)

    def test_segment_identity_mismatch_rejected(self):
        self.segment(1, 100, 500)
        changed = dict(self.identity, case_id="different")
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            runner.rebuild_segment_manifest(self.root, changed)

    def test_corrupt_checkpoint_falls_back_to_other_slot(self):
        folder = self.root / "attempts" / "attempt_0001" / "checkpoints"
        for slot, step in ((0, 300), (1, 500)):
            blob = f"fake checkpoint slot {slot}".encode()
            runner.atomic_bytes(folder / f"slot_{slot}.chk", blob)
            runner.atomic_json(folder / f"slot_{slot}.json", {
                "run_identity_sha256": self.digest, "absolute_step": step,
                "saved_utc": "2026-01-01T00:00:00Z",
                "checkpoint_sha256": hashlib.sha256(blob).hexdigest()})
        (folder / "slot_1.chk").write_bytes(b"corrupted")
        checkpoint, rejected = runner.find_checkpoint(self.root, self.identity)
        self.assertEqual(checkpoint["metadata"]["absolute_step"], 300)
        self.assertEqual(len(rejected), 1)
        self.assertIn("binary/metadata mismatch", rejected[0]["reason"])

    def test_production_without_review_rejected(self):
        with self.assertRaisesRegex(ValueError, "explicit pilot review"):
            runner.verify_review(self.root, self.identity)

    def test_duration_must_be_integral_and_positive(self):
        self.assertEqual(runner.exact_steps(1), 500)
        for bad in (0, -1, 0.001, float("inf")):
            with self.assertRaises(ValueError):
                runner.exact_steps(bad)

    def test_case_path_cannot_escape(self):
        for bad in ("..", "../case", "case/other", "case\\other"):
            with self.assertRaises(ValueError):
                runner.case_path(self.root, bad)


if __name__ == "__main__":
    unittest.main()
