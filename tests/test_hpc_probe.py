"""Probe tests stub GPU subprocess calls; no remote access and no MD."""
import contextlib
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import hpc_probe as probe


class ProbeTests(unittest.TestCase):
    def setUp(self):
        # Windows platform queries can invoke subprocess internally. Isolate
        # them so GPU subprocess stubs cannot intercept host discovery.
        node_patch = patch.object(probe.platform, "node", return_value="toy-host")
        platform_patch = patch.object(probe.platform, "platform", return_value="toy-platform")
        node_patch.start()
        platform_patch.start()
        self.addCleanup(node_patch.stop)
        self.addCleanup(platform_patch.stop)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.state = self.root / "state.json"
        self.state.write_text(json.dumps({
            "state": "running", "worker_pid": 42,
            "current_progress": {"production_ns": 3.1, "secret": "do-not-copy"},
            "completed": [{"case_id": "toy", "seed": 11, "run": "/private/path"}],
            "token": "do-not-copy", "commands": ["never execute me"],
        }), encoding="utf-8")

    def test_snapshot_whitelists_state_and_queries_only_gpu(self):
        before = self.state.read_bytes()
        with patch.object(probe.shutil, "which", return_value="/trusted/nvidia-smi"), patch.object(
            probe.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "Toy GPU, 1, 40, 0, 0, 8\n", "")
        ) as run:
            result = probe.collect(self.root, "state.json", timeout=2)
        call_args, call_kwargs = run.call_args
        self.assertEqual(call_args[0][0], "/trusted/nvidia-smi")
        self.assertEqual(call_args[0][1], "--query-gpu=name,driver_version,temperature.gpu,utilization.gpu,memory.used,memory.total")
        self.assertFalse(call_kwargs["shell"])
        self.assertEqual(call_kwargs["timeout"], 2)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(result["snapshot_status"], "collected")
        self.assertEqual(result["state"]["reported"]["current_progress"], {"production_ns": 3.1})
        self.assertNotIn("do-not-copy", json.dumps(result))
        self.assertNotIn("/private/path", json.dumps(result))
        self.assertEqual(before, self.state.read_bytes())
        self.assertEqual(list(self.root.iterdir()), [self.state])

    def test_missing_gpu_and_invalid_state_are_partial_not_md_failure(self):
        self.state.write_text("{partial", encoding="utf-8")
        with patch.object(probe.shutil, "which", return_value=None), patch.object(probe.subprocess, "run") as run:
            result = probe.collect(self.root, self.state)
        run.assert_not_called()
        self.assertEqual(result["snapshot_status"], "partial")
        self.assertEqual(result["state"]["status"], "unavailable")
        self.assertNotIn("reported", result["state"])
        self.assertTrue(any("does not establish MD failure" in x for x in result["limitations"]))

    def test_gpu_timeout_is_bounded_and_reported(self):
        with patch.object(probe.shutil, "which", return_value="nvidia-smi"), patch.object(
            probe.subprocess, "run", side_effect=subprocess.TimeoutExpired("nvidia-smi", 0.25)
        ):
            result = probe.query_gpu(0.25)
        self.assertEqual(result, {"status": "timed_out", "timeout_seconds": 0.25})

    def test_gpu_nonzero_exit_is_not_success(self):
        with patch.object(probe.shutil, "which", return_value="nvidia-smi"), patch.object(
            probe.subprocess, "run", return_value=subprocess.CompletedProcess([], 9, "", "driver error")
        ):
            result = probe.query_gpu(1)
        self.assertEqual(result["status"], "query_failed")
        self.assertEqual(result["exit_code"], 9)

    def test_refuse_overwrite_and_only_new_receipt_is_written(self):
        out = self.root / "receipt.json"
        with patch.object(probe.shutil, "which", return_value=None), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(probe.main(["--workdir", str(self.root), "--out", str(out)]), 1)
        original = out.read_bytes()
        with patch.object(probe, "collect") as collect, contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(probe.main(["--workdir", str(self.root), "--out", str(out)]), 2)
        collect.assert_not_called()
        self.assertEqual(out.read_bytes(), original)
        self.assertEqual(set(self.root.iterdir()), {self.state, out})

    def test_size_limit_nonfinite_and_nonscalar_state(self):
        for content in (" " * (probe.MAX_STATE_BYTES + 1), '{"state": NaN}', "[]"):
            self.state.write_text(content, encoding="utf-8")
            result = probe.read_state(self.state, datetime.now(timezone.utc))
            self.assertEqual(result["status"], "unavailable")
            self.assertNotIn("reported", result)

    def test_state_rewritten_during_read_is_unknown(self):
        actual = self.state.stat()
        before = SimpleNamespace(st_size=actual.st_size, st_mtime_ns=100,
                                 st_ino=actual.st_ino, st_mode=actual.st_mode)
        after = SimpleNamespace(st_size=actual.st_size, st_mtime_ns=101,
                                st_ino=actual.st_ino, st_mode=actual.st_mode)
        # read_state stat / is_file stat / post-read stat.
        with patch.object(Path, "stat", side_effect=[before, before, after]):
            result = probe.read_state(self.state, datetime.now(timezone.utc))
        self.assertEqual(result["status"], "changed_during_read")
        self.assertNotIn("reported", result)
        self.assertNotIn("file_sha256", result)

    def test_reject_invalid_timeouts_before_subprocess(self):
        with patch.object(probe.subprocess, "run") as run:
            for timeout in (0, -1, 31, float("nan"), float("inf")):
                with self.assertRaises(ValueError):
                    probe.collect(self.root, timeout=timeout)
        run.assert_not_called()

    def test_missing_workdir_is_not_created(self):
        missing = self.root / "absent"
        with self.assertRaises(FileNotFoundError):
            probe.collect(missing)
        self.assertFalse(missing.exists())


if __name__ == "__main__":
    unittest.main()
