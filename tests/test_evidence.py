"""Synthetic, standard-library tests; these never execute an MD engine."""
import contextlib
import copy
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import assemble_evidence as evidence


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name in ("synthetic_evidence.json", "toy_parent.smi", "toy_variant.smi"):
            shutil.copyfile(ROOT / "examples" / name, self.root / name)
        self.path = self.root / "synthetic_evidence.json"
        self.data = json.loads(self.path.read_text(encoding="utf-8"))

    def write(self):
        self.path.write_text(json.dumps(self.data), encoding="utf-8")

    def test_demo_retains_missing_and_provisional_evidence(self):
        data, rows = evidence.assemble(self.path)
        self.assertEqual(len(rows), 4)
        self.assertEqual(rows[1]["MD_evidence_status"], "incomplete_snapshot")
        self.assertIsNone(rows[2]["MD_sampling_ns"])
        self.assertIsNone(rows[3]["supplied_docking_score_kcal_mol"])
        self.assertTrue(all(r["human_decision"] == "not_recorded" for r in rows))
        output = self.root / "output"
        with contextlib.redirect_stdout(io.StringIO()):
            receipt = evidence.main(["--input", str(self.path), "--out", str(output)])
        self.assertFalse(receipt["new_MD"])
        self.assertFalse(receipt["new_AI_inference"])
        for name, digest in receipt["output_sha256"].items():
            self.assertEqual(evidence.sha(output / name), digest)
        with self.assertRaisesRegex(ValueError, "new or empty"):
            evidence.main(["--input", str(self.path), "--out", str(output)])

    def test_tampered_structure_rejected(self):
        (self.root / "toy_parent.smi").write_text("changed", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "SHA mismatch"):
            evidence.assemble(self.path)

    def test_duplicate_or_missing_pair_rejected(self):
        original = copy.deepcopy(self.data)
        self.data["evidence"].append(self.data["evidence"][0])
        self.write()
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            evidence.assemble(self.path)
        self.data = original
        self.data["evidence"].pop()
        self.write()
        with self.assertRaisesRegex(ValueError, "every candidate-target"):
            evidence.assemble(self.path)

    def test_nonfinite_metrics_rejected(self):
        self.data["evidence"][0]["docking_score_kcal_mol"] = float("nan")
        self.write()
        with self.assertRaisesRegex(ValueError, "Invalid JSON numeric"):
            evidence.assemble(self.path)

    def test_missing_md_cannot_be_presented_as_results(self):
        self.data["evidence"][2]["md"]["sampling_ns"] = 20
        self.write()
        with self.assertRaisesRegex(ValueError, "Missing MD"):
            evidence.assemble(self.path)

    def test_structure_path_escape_rejected(self):
        self.data["candidates"][0]["structure_file"] = "../outside.smi"
        self.write()
        with self.assertRaisesRegex(ValueError, "escapes"):
            evidence.assemble(self.path)


if __name__ == "__main__":
    unittest.main()
