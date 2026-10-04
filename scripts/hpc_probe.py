"""Bounded, standard-library HPC telemetry; never starts or restores MD.

Reads one explicitly supplied state JSON and queries nvidia-smi. Apart from an
optional NEW output file, no files are written. Output may contain host/path
information and is for private operations, not automatic public upload.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import sys

MAX_STATE_BYTES = 1024 * 1024
MAX_OUTPUT_CHARS = 4096
STATE_KEYS = (
    "state", "status", "observed_utc", "updated_at", "observed_at",
    "current_case", "current_seed", "current_index", "worker_pid",
    "active_runner_pid", "automatic_failure_restart", "scientific_analysis",
)
PROGRESS_KEYS = ("production_ns", "absolute_step", "attempt_status", "disk_free_GiB")
COMPLETED_KEYS = ("case_id", "seed", "planned_production_ns", "ns")
LIMITATIONS = [
    "State and completed records are reported claims, not independent MD validation.",
    "PID values are not process-liveness or process-identity checks.",
    "No checkpoint, trajectory, manifest, environment, or scientific-quality audit is performed.",
    "A GPU query failure or stale state does not establish MD failure.",
    "No MD, OpenMM import/Context, installation, restart, lock removal, or queue action is performed.",
]


def utc_now():
    return datetime.now(timezone.utc)


def scalar_subset(data, keys):
    """Return a bounded whitelist; never copy commands, paths or nested secrets."""
    result = {}
    if not isinstance(data, dict):
        return result
    for key in keys:
        value = data.get(key)
        if key not in data:
            continue
        if value is None or isinstance(value, (bool, int)):
            result[key] = value
        elif isinstance(value, float) and math.isfinite(value):
            result[key] = value
        elif isinstance(value, str):
            result[key] = value[:500]
    return result


def reject_json_constant(value):
    raise ValueError("Non-finite JSON constant: " + value)


def read_state(path, observed_at):
    """Read one bounded file. A concurrent rewrite is unknown, never completed."""
    result = {"path": str(path), "status": "unavailable"}
    try:
        before = path.stat()
        if not path.is_file():
            raise ValueError("State must be a regular file")
        if before.st_size > MAX_STATE_BYTES:
            raise ValueError("State exceeds 1 MiB size limit")
        with path.open("rb") as stream:
            raw = stream.read(MAX_STATE_BYTES + 1)
        after = path.stat()
        if len(raw) > MAX_STATE_BYTES:
            raise ValueError("State exceeds 1 MiB size limit")
        if (before.st_size, before.st_mtime_ns, before.st_ino) != (
            after.st_size, after.st_mtime_ns, after.st_ino
        ):
            result.update(status="changed_during_read", retry_policy="Inspect later; no automatic loop")
            return result
        data = json.loads(raw.decode("utf-8-sig"), parse_constant=reject_json_constant)
        if not isinstance(data, dict):
            raise ValueError("State JSON must be an object")
        reported = scalar_subset(data, STATE_KEYS)
        if isinstance(data.get("current_progress"), dict):
            reported["current_progress"] = scalar_subset(data["current_progress"], PROGRESS_KEYS)
        completed = data.get("completed")
        if isinstance(completed, list):
            reported["completed_record_count"] = len(completed)
            reported["completed_records"] = [scalar_subset(row, COMPLETED_KEYS) for row in completed[:100]]
            reported["completed_records_truncated"] = len(completed) > 100
        result.update(
            status="read",
            file_sha256=hashlib.sha256(raw).hexdigest(),
            bytes=len(raw),
            file_modified_utc=datetime.fromtimestamp(after.st_mtime, timezone.utc).isoformat(),
            file_age_seconds=round(observed_at.timestamp() - after.st_mtime, 3),
            reported=reported,
            consistency_scope="File metadata unchanged around read; not a transaction or checkpoint audit",
        )
    except (OSError, UnicodeError, ValueError, RecursionError) as exc:
        result["error"] = str(exc)[:MAX_OUTPUT_CHARS]
    return result


def query_gpu(timeout):
    executable = shutil.which("nvidia-smi")
    if not executable:
        return {"status": "unavailable", "reason": "nvidia-smi not found on PATH"}
    command = [executable,
               "--query-gpu=name,driver_version,temperature.gpu,utilization.gpu,memory.used,memory.total",
               "--format=csv,noheader,nounits"]
    try:
        completed = subprocess.run(
            command, stdin=subprocess.DEVNULL, capture_output=True,
            text=True, encoding="utf-8", errors="replace", timeout=timeout,
            check=False, shell=False,
        )
        return {
            "status": "queried" if completed.returncode == 0 else "query_failed",
            "exit_code": completed.returncode,
            "columns": ["name", "driver_version", "temperature_C", "utilization_percent", "memory_used_MiB", "memory_total_MiB"],
            "csv": completed.stdout[:MAX_OUTPUT_CHARS].strip(),
            "stderr": completed.stderr[:MAX_OUTPUT_CHARS].strip(),
            "output_truncated": max(len(completed.stdout), len(completed.stderr)) > MAX_OUTPUT_CHARS,
        }
    except subprocess.TimeoutExpired:
        return {"status": "timed_out", "timeout_seconds": timeout}
    except OSError as exc:
        return {"status": "unavailable", "reason": str(exc)[:MAX_OUTPUT_CHARS]}


def collect(workdir, state_path=None, timeout=8.0):
    if not math.isfinite(timeout) or not 0 < timeout <= 30:
        raise ValueError("GPU timeout must be greater than 0 and at most 30 seconds")
    workdir = Path(workdir).resolve(strict=True)
    if not workdir.is_dir():
        raise ValueError("Workdir must be an existing directory")
    started = utc_now()
    disk = shutil.disk_usage(workdir)
    result = {
        "schema_version": 1,
        "observed_utc": started.isoformat(),
        "host": platform.node(),
        "platform": platform.platform(),
        "workdir": str(workdir),
        "scope": "Read-only telemetry snapshot; optional new receipt is the only write",
        "disk": {"total_bytes": disk.total, "used_bytes": disk.used, "free_bytes": disk.free},
        "gpu": query_gpu(timeout),
        "state": {"status": "not_requested"},
        "limitations": LIMITATIONS,
        "public_upload_safe": False,
    }
    if state_path is not None:
        selected = Path(state_path)
        if not selected.is_absolute():
            selected = workdir / selected
        result["state"] = read_state(selected, started)
    result["finished_utc"] = utc_now().isoformat()
    result["snapshot_status"] = (
        "collected" if result["gpu"]["status"] == "queried"
        and result["state"]["status"] in ("read", "not_requested") else "partial"
    )
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True, type=Path, help="Existing computation directory; never created")
    parser.add_argument("--state", type=Path, help="Optional JSON; relative paths are resolved under --workdir")
    parser.add_argument("--gpu-timeout", type=float, default=8.0, help="GPU query timeout in seconds (0, 30]")
    parser.add_argument("--out", type=Path, help="Optional NEW local receipt file; existing paths are never overwritten")
    args = parser.parse_args(argv)
    try:
        if args.out is not None and (args.out.exists() or args.out.is_symlink()):
            raise FileExistsError("Output must be a new file; refusing overwrite")
        result = collect(args.workdir, args.state, args.gpu_timeout)
        payload = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        if args.out is not None:
            # Exclusive creation also protects a race after the preflight check.
            # No parent mkdir or atomic rename on a potentially remote mount.
            with args.out.open("x", encoding="utf-8") as stream:
                stream.write(payload)
        print(payload, end="")
        return 0 if result["snapshot_status"] == "collected" else 1
    except (OSError, ValueError) as exc:
        print("hpc_probe: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
