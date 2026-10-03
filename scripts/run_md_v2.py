"""Run prebuilt OpenMM systems with immutable attempts and auditable restarts.

No parameterization is performed. Each restart creates a new trajectory segment.
Read segment_manifest.json before combining trajectories; superseded tails must
be excluded. A reviewed, completed pilot is required before production.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import platform as pyplatform
import time
import uuid


ROOT = Path(__file__).resolve().parents[1]
SETTINGS = {
    "timestep_ps": 0.002,
    "temperature_K": 310.0,
    "friction_per_ps": 1.0,
    "pressure_bar": 1.0,
    "barostat_interval_steps": 25,
    "minimization_iterations": 2500,
    "minimization_tolerance_kj_mol_nm": 10.0,
    "pilot_equilibration_ps_k": [[2, 1000], [4, 100], [4, 10]],
    "production_equilibration_ps_k": [[100, 1000], [100, 500], [100, 100], [100, 10], [100, 0]],
    "pilot_production_ns": 0.01,
    "pilot_report_steps": 500,
    "production_report_steps": 5000,
    "checkpoint_steps": 50000,
}
INPUT_FILES = ("system.xml", "complex.cif", "positions_nm.npy", "metadata.json")
REVIEW_CHECKS = ("finite_energy", "temperature_density", "geometry", "throughput")


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def atomic_bytes(path, content):
    """Replace only this output file; a killed writer leaves its prior version."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".tmp-" + uuid.uuid4().hex[:12])
    try:
        with temporary.open("xb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_json(path, value):
    atomic_bytes(path, (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8"))


@contextmanager
def run_lock(directory):
    """OS advisory locks release automatically if the process is killed."""
    directory.mkdir(parents=True, exist_ok=True)
    stream = (directory / ".runner.lock").open("a+b")
    stream.seek(0, 2)
    if stream.tell() == 0:
        stream.write(b"0")
        stream.flush()
    stream.seek(0)
    locked = False
    try:
        if os.name == "nt":
            import msvcrt
            try:
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError as exc:
                raise RuntimeError("Another process owns this run") from exc
        else:
            import fcntl
            try:
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as exc:
                raise RuntimeError("Another process owns this run") from exc
        locked = True
        yield
    finally:
        if locked:
            stream.seek(0)
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        stream.close()


def case_path(root, case):
    if not case or Path(case).name != case or "/" in case or "\\" in case or case in (".", ".."):
        raise ValueError("Invalid case identifier")
    return root / "systems" / case


def input_hashes(root, case):
    source = case_path(root, case)
    hashes = {name: sha(source / name) for name in INPUT_FILES}
    meta = read_json(source / "metadata.json")
    if meta.get("case_id") != case:
        raise ValueError("Metadata case identity differs from requested case")
    if meta.get("system_sha256") != hashes["system.xml"]:
        raise ValueError("System hash differs from preparation metadata")
    # Verify the frozen package as well when its manifest is available.
    package = root / "package_manifest.json"
    if package.exists():
        entries = read_json(package)
        if isinstance(entries, dict):
            entries = entries.get("files", [])
        expected = {e["path"].replace("\\", "/"): e["sha256"] for e in entries}
        for filename, actual in hashes.items():
            relative = f"systems/{case}/{filename}"
            if relative not in expected or expected[relative].lower() != actual:
                raise ValueError(f"Frozen package input mismatch: {relative}")
    return hashes


def exact_steps(duration_ps):
    value = Decimal(str(duration_ps)) / Decimal(str(SETTINGS["timestep_ps"]))
    if not value.is_finite() or value != value.to_integral_value() or value <= 0:
        raise ValueError("Duration must be a positive integral number of steps")
    return int(value)


def make_identity(root, case, mode, seed, platform, ns, openmm_version):
    if mode not in ("pilot", "production") or seed <= 0 or seed > 2_000_000_000:
        raise ValueError("Invalid mode or random seed")
    requested_ns = SETTINGS["pilot_production_ns"] if mode == "pilot" else ns
    production_steps = exact_steps(Decimal(str(requested_ns)) * 1000)
    schedule = SETTINGS[f"{mode}_equilibration_ps_k"]
    eq_steps = sum(exact_steps(ps) for ps, _ in schedule)
    props = {"Precision": "mixed"} if platform in ("CUDA", "OpenCL") else {"Threads": "4"} if platform == "CPU" else {}
    spec = {
        "schema": 2, "case_id": case, "mode": mode, "seed": seed,
        "input_sha256": input_hashes(root, case),
        "production_ns": str(Decimal(str(requested_ns)).normalize()),
        "production_steps": production_steps, "equilibration_steps": eq_steps,
        "target_absolute_step": eq_steps + production_steps,
        "settings": json.loads(json.dumps(SETTINGS)),
        "platform": platform, "platform_properties": props,
        "openmm_version": openmm_version,
        "python_major_minor": list(__import__("sys").version_info[:2]),
    }
    return spec


def attempt_directories(run_dir):
    return sorted(p for p in (run_dir / "attempts").glob("attempt_*") if p.is_dir())


def rebuild_segment_manifest(run_dir, identity):
    """A later segment replaces every earlier frame beyond its restart step."""
    segments = []
    for folder in attempt_directories(run_dir):
        path = folder / "segment.json"
        if not path.exists():
            continue
        segment = read_json(path)
        if segment["run_identity_sha256"] != digest(identity):
            raise ValueError("Segment identity mismatch")
        segments.append(segment)
    for index, segment in enumerate(segments):
        last = segment.get("last_output_step")
        future_starts = [later["start_after_absolute_step"] for later in segments[index + 1:]]
        cap = min([last] + future_starts) if last is not None else None
        segment["include_steps_gt"] = segment["start_after_absolute_step"]
        segment["include_steps_le"] = cap
        segment["superseded_tail_after_step"] = cap if last is not None and cap < last else None
        interval = segment["report_interval_steps"]
        first = (segment["include_steps_gt"] // interval + 1) * interval
        segment["include_first_frame_index_zero_based"] = 0
        segment["include_frame_count"] = max(0, (cap - first) // interval + 1) if cap is not None else 0
        segment["ignore_trailing_uncommitted_frames"] = True
    payload = {
        "schema": 2, "run_identity_sha256": digest(identity), "updated_utc": utcnow(),
        "time_unit": "ps", "timestep_ps": identity["settings"]["timestep_ps"],
        "equilibration_steps": identity["equilibration_steps"],
        "target_absolute_step": identity["target_absolute_step"],
        "instructions": "Use segments in listed order. Include only the first include_frame_count paired CSV/DCD frames from each segment. Absolute CSV steps are authoritative. Ignore superseded tails and any uncommitted trailing frame; never concatenate complete DCD files blindly.",
        "segments": segments,
    }
    atomic_json(run_dir / "segment_manifest.json", payload)
    return payload


def find_checkpoint(run_dir, identity):
    candidates, rejected = [], []
    for attempt in attempt_directories(run_dir):
        for path in (attempt / "checkpoints").glob("slot_*.json"):
            try:
                metadata = read_json(path)
                binary = path.with_suffix(".chk")
                if metadata["run_identity_sha256"] != digest(identity):
                    raise ValueError("checkpoint identity mismatch")
                if metadata["checkpoint_sha256"] != sha(binary):
                    raise ValueError("checkpoint binary/metadata mismatch")
                if not identity["equilibration_steps"] <= metadata["absolute_step"] <= identity["target_absolute_step"]:
                    raise ValueError("checkpoint step outside requested run")
                candidates.append((metadata["absolute_step"], metadata["saved_utc"], binary, metadata))
            except (OSError, ValueError, KeyError, TypeError) as exc:
                rejected.append({"path": str(path.relative_to(run_dir)), "reason": str(exc)})
    if not candidates:
        return None, rejected
    _, _, binary, metadata = max(candidates, key=lambda entry: (entry[0], entry[1]))
    return {"binary": binary, "metadata": metadata}, rejected


def verify_review(root, production_identity):
    case = production_identity["case_id"]
    review_path = root / "reviews" / (case + ".json")
    if not review_path.exists():
        raise ValueError("Production requires an explicit pilot review; use the review action after inspecting the pilot")
    review = read_json(review_path)
    if review.get("approved") is not True or review.get("case_id") != case:
        raise ValueError("Pilot review is not approved for this case")
    if review.get("input_sha256") != production_identity["input_sha256"]:
        raise ValueError("Pilot review refers to different inputs")
    if not all(review.get("checks", {}).get(key) is True for key in REVIEW_CHECKS):
        raise ValueError("Required pilot review checks are incomplete")
    pilot_dir = root / "runs_v2" / "pilot" / case / f"seed_{int(review['pilot_seed'])}"
    pilot_identity = read_json(pilot_dir / "run_identity.json")
    if review.get("pilot_identity_sha256") != digest(pilot_identity):
        raise ValueError("Reviewed pilot identity changed")
    for field in ("input_sha256", "platform", "platform_properties", "openmm_version", "settings"):
        if pilot_identity[field] != production_identity[field]:
            raise ValueError(f"Pilot review does not match production {field}")
    attempt = review["pilot_attempt"]
    if Path(attempt).name != attempt or not attempt.startswith("attempt_"):
        raise ValueError("Invalid reviewed pilot attempt")
    status = read_json(pilot_dir / "attempts" / attempt / "attempt_status.json")
    if status.get("status") != "completed" or status.get("completed_absolute_step") != pilot_identity["target_absolute_step"]:
        raise ValueError("The reviewed pilot did not complete")
    if status.get("run_identity_sha256") != digest(pilot_identity):
        raise ValueError("Reviewed pilot status identity differs")
    return {"path": str(review_path.relative_to(root)), "sha256": sha(review_path)}


def approve_pilot(root, case, seed, reviewer, note, checks, attempt=None):
    """Record a review decision; this never infers scientific approval from completion."""
    if not reviewer.strip() or not note.strip() or not all(checks.get(k) is True for k in REVIEW_CHECKS):
        raise ValueError("Provide reviewer, a substantive note, and all four explicit review checks")
    pilot_dir = root / "runs_v2" / "pilot" / case / f"seed_{seed}"
    identity = read_json(pilot_dir / "run_identity.json")
    if identity["input_sha256"] != input_hashes(root, case):
        raise ValueError("Pilot inputs have changed since the run")
    completed, missing_status = [], []
    for folder in attempt_directories(pilot_dir):
        status_path = folder / "attempt_status.json"
        if not status_path.exists():
            missing_status.append(folder.name)
            continue
        status = read_json(status_path)
        if status.get("status") == "completed" and status.get("completed_absolute_step") == identity["target_absolute_step"]:
            completed.append(folder)
    selected = [p for p in completed if attempt is None or p.name == attempt]
    if not selected:
        raise ValueError("No matching completed pilot attempt")
    selected = selected[-1]
    review = {
        "schema": 2, "approved": True, "case_id": case, "pilot_seed": seed,
        "pilot_attempt": selected.name, "pilot_identity_sha256": digest(identity),
        "input_sha256": identity["input_sha256"], "checks": checks,
        "reviewer": reviewer, "note": note, "reviewed_utc": utcnow(),
        "skipped_attempts_missing_status": missing_status,
        "scope": "Technical pilot review only; not a binding, affinity, convergence or activity conclusion",
    }
    folder = root / "reviews"
    folder.mkdir(exist_ok=True)
    # Retain every decision, while the case pointer identifies the current review.
    archive = folder / f"{case}.{uuid.uuid4().hex}.json"
    atomic_json(archive, review)
    atomic_json(folder / (case + ".json"), review)
    return review


class SegmentReporter:
    def __init__(self, app, folder, run_dir, identity, start_step, checkpoint_source):
        self.folder = folder
        self.interval = identity["settings"][f"{identity['mode']}_report_steps"]
        self.csv = app.StateDataReporter(str(folder / "state.csv"), self.interval, step=True, time=True,
            potentialEnergy=True, temperature=True, density=True, speed=True, separator=",")
        self.dcd = app.DCDReporter(str(folder / "trajectory.dcd"), self.interval, enforcePeriodicBox=False)
        self.metadata = {
            "attempt": folder.name, "run_identity_sha256": digest(identity),
            "trajectory": str((folder / "trajectory.dcd").relative_to(run_dir)).replace("\\", "/"),
            "state_csv": str((folder / "state.csv").relative_to(run_dir)).replace("\\", "/"),
            "start_after_absolute_step": start_step, "last_output_step": None,
            "paired_frames_committed": 0, "report_interval_steps": self.interval,
            "resumed_from": checkpoint_source,
        }
        atomic_json(folder / "segment.json", self.metadata)

    def describeNextReport(self, simulation):
        return {"steps": self.interval - simulation.currentStep % self.interval,
                "periodic": False, "include": ["positions", "energy"]}

    def report(self, simulation, state):
        self.csv.report(simulation, state)
        self.dcd.report(simulation, state)
        for reporter in (self.csv, self.dcd):
            stream = getattr(reporter, "_out", None)
            if stream is not None and hasattr(stream, "flush"):
                stream.flush()
        # Commit the pair only after both writes. A killed writer may leave an
        # extra trailing frame, which the manifest explicitly excludes.
        self.metadata["last_output_step"] = simulation.currentStep
        self.metadata["paired_frames_committed"] += 1
        atomic_json(self.folder / "segment.json", self.metadata)

    def close(self):
        for reporter in (self.csv, self.dcd):
            stream = getattr(reporter, "_out", None)
            if stream is not None and hasattr(stream, "close") and not stream.closed:
                stream.close()


class AtomicCheckpointReporter:
    def __init__(self, folder, run_dir, identity, record):
        self.folder, self.run_dir, self.identity, self.record = folder, run_dir, identity, record
        self.interval = identity["settings"]["checkpoint_steps"]
        self.sequence = 0

    def describeNextReport(self, simulation):
        return {"steps": self.interval - simulation.currentStep % self.interval,
                "periodic": None, "include": []}

    def report(self, simulation, state=None):
        slot = self.sequence % 2
        self.sequence += 1
        target = self.folder / "checkpoints" / f"slot_{slot}.chk"
        blob = simulation.context.createCheckpoint()
        metadata = {"schema": 2, "run_identity_sha256": digest(self.identity),
            "absolute_step": simulation.currentStep, "saved_utc": utcnow(),
            "checkpoint_sha256": hashlib.sha256(blob).hexdigest(),
            "phase": "production", "attempt": self.folder.name}
        atomic_bytes(target, blob)
        # If interrupted between these replaces, this slot is rejected by hash
        # and the other independently committed slot remains usable.
        atomic_json(target.with_suffix(".json"), metadata)
        self.record["last_checkpoint_absolute_step"] = simulation.currentStep
        atomic_json(self.folder / "attempt_status.json", self.record)
        rebuild_segment_manifest(self.run_dir, self.identity)


def run(case, mode, seed=11, platform="CUDA", ns=20, resume=False, root=ROOT, restart_from_inputs=False):
    import numpy as np
    import openmm as mm
    from openmm import app, unit

    root = Path(root).resolve()
    identity = make_identity(root, case, mode, seed, platform, ns, mm.__version__)
    source = case_path(root, case)
    directory = root / "runs_v2" / mode / case / f"seed_{seed}"
    with run_lock(directory):
        identity_path = directory / "run_identity.json"
        if identity_path.exists():
            old = read_json(identity_path)
            if old != identity:
                changed = [key for key in identity if old.get(key) != identity[key]]
                raise ValueError("Refusing changed run identity: " + ", ".join(changed))
            summary_path = directory / "run_summary.json"
            if summary_path.exists():
                summary = read_json(summary_path)
                if summary.get("status") == "completed":
                    if summary.get("run_identity_sha256") != digest(identity) or summary.get("completed_absolute_step") != identity["target_absolute_step"]:
                        raise ValueError("Completed summary does not match requested identity/length")
                    completed_attempt = summary.get("completed_attempt", "")
                    if Path(completed_attempt).name != completed_attempt or not completed_attempt.startswith("attempt_"):
                        raise ValueError("Invalid completed attempt in run summary")
                    source_status_path = directory / "attempts" / completed_attempt / "attempt_status.json"
                    if not source_status_path.exists():
                        raise ValueError("Completed summary has no supporting attempt status")
                    source_status = read_json(source_status_path)
                    if (source_status.get("status") != "completed"
                        or source_status.get("run_identity_sha256") != digest(identity)
                        or source_status.get("completed_absolute_step") != identity["target_absolute_step"]):
                        raise ValueError("Completed summary disagrees with supporting attempt status; inspect before retrying")
                    return {**summary, "idempotent_skip": True}
            if not resume:
                raise FileExistsError("Existing incomplete run. Use --resume; old outputs will remain untouched")
        else:
            if resume:
                raise ValueError("Cannot resume a run that has no recorded identity")

        review = verify_review(root, identity) if mode == "production" else None
        checkpoint, rejected = find_checkpoint(directory, identity) if resume else (None, [])
        if resume and checkpoint is None and not restart_from_inputs:
            raise ValueError("No valid production checkpoint. Inspect failure, then explicitly use --resume --restart-from-inputs if a fresh attempt is appropriate")
        if restart_from_inputs and (not resume or checkpoint is not None):
            raise ValueError("--restart-from-inputs only applies to an existing run with no valid checkpoint")
        if not identity_path.exists():
            atomic_json(identity_path, identity)
        previous = attempt_directories(directory)
        index = max([int(p.name.split("_")[-1]) for p in previous] + [0]) + 1
        folder = directory / "attempts" / f"attempt_{index:04d}"
        folder.mkdir(parents=True, exist_ok=False)
        source_checkpoint = None
        if checkpoint:
            source_checkpoint = {**checkpoint["metadata"],
                "path": str(checkpoint["binary"].relative_to(directory)).replace("\\", "/")}
        record = {"schema": 2, "case_id": case, "mode": mode, "seed": seed,
            "attempt": folder.name, "run_identity_sha256": digest(identity), "status": "running",
            "started_utc": utcnow(), "openmm": mm.__version__, "python": pyplatform.python_version(),
            "platform": platform, "resumed_from": source_checkpoint, "rejected_checkpoints": rejected,
            "restart_from_inputs": bool(restart_from_inputs), "pilot_review": review,
            "is_scientific_production": mode == "production"}
        atomic_json(folder / "attempt_status.json", record)
        simulation, segment, checkpoint_reporter, equilibration = None, None, None, None
        try:
            system = mm.XmlSerializer.deserialize((source / "system.xml").read_text(encoding="utf-8"))
            barostats = [f for f in system.getForces() if isinstance(f, mm.MonteCarloBarostat)]
            if len(barostats) != 1:
                raise ValueError("Expected exactly one MonteCarloBarostat")
            barostat = barostats[0]
            barostat.setRandomNumberSeed(seed + 100000)
            barostat.setFrequency(SETTINGS["barostat_interval_steps"] if checkpoint else 0)
            integrator = mm.LangevinMiddleIntegrator(SETTINGS["temperature_K"] * unit.kelvin,
                SETTINGS["friction_per_ps"] / unit.picosecond, SETTINGS["timestep_ps"] * unit.picosecond)
            integrator.setRandomNumberSeed(seed)
            pdb = app.PDBxFile(str(source / "complex.cif"))
            selected_platform = mm.Platform.getPlatformByName(platform)
            simulation = app.Simulation(pdb.topology, system, integrator, selected_platform, identity["platform_properties"])
            if checkpoint:
                simulation.loadCheckpoint(str(checkpoint["binary"]))
                if simulation.currentStep != source_checkpoint["absolute_step"]:
                    raise ValueError("Checkpoint binary step differs from its committed metadata")
                for parameter in ("k_protein", "k_ligand"):
                    if abs(simulation.context.getParameter(parameter)) > 1e-12:
                        raise ValueError("Production checkpoint unexpectedly contains position restraints")
            else:
                positions = np.load(source / "positions_nm.npy", allow_pickle=False)
                if positions.shape != (system.getNumParticles(), 3) or not np.isfinite(positions).all():
                    raise ValueError("Invalid input coordinate array")
                simulation.context.setPositions(positions * unit.nanometer)
                simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
                simulation.minimizeEnergy(maxIterations=SETTINGS["minimization_iterations"],
                    tolerance=SETTINGS["minimization_tolerance_kj_mol_nm"] * unit.kilojoule_per_mole / unit.nanometer)
                simulation.context.setVelocitiesToTemperature(SETTINGS["temperature_K"] * unit.kelvin, seed)
                equilibration = app.StateDataReporter(str(folder / "equilibration.csv"),
                    min(500, identity["settings"][f"{mode}_report_steps"]), step=True, time=True,
                    potentialEnergy=True, temperature=True, density=True, separator=",")
                simulation.reporters.append(equilibration)
                for schedule_index, (ps, force_constant) in enumerate(SETTINGS[f"{mode}_equilibration_ps_k"]):
                    if schedule_index == 1:
                        barostat.setFrequency(SETTINGS["barostat_interval_steps"])
                        simulation.context.reinitialize(preserveState=True)
                    for parameter in ("k_protein", "k_ligand"):
                        simulation.context.setParameter(parameter, force_constant)
                    simulation.step(exact_steps(ps))
                simulation.reporters = []
                for parameter in ("k_protein", "k_ligand"):
                    simulation.context.setParameter(parameter, 0)
            start_step = simulation.currentStep
            if start_step < identity["equilibration_steps"] or start_step > identity["target_absolute_step"]:
                raise ValueError("Simulation step outside declared production interval")
            segment = SegmentReporter(app, folder, directory, identity, start_step, source_checkpoint)
            simulation.reporters.append(segment)
            checkpoint_reporter = AtomicCheckpointReporter(folder, directory, identity, record)
            simulation.reporters.append(checkpoint_reporter)
            rebuild_segment_manifest(directory, identity)
            remaining = identity["target_absolute_step"] - start_step
            timer = time.perf_counter()
            if remaining:
                simulation.step(remaining)
            elapsed = time.perf_counter() - timer
            final = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            energy = float(final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            if not np.isfinite(energy):
                raise ValueError("Nonfinite final potential energy")
            checkpoint_reporter.report(simulation)
            pdb.topology.setPeriodicBoxVectors(final.getPeriodicBoxVectors())
            with (folder / "final.cif").open("w", encoding="utf-8") as stream:
                app.PDBxFile.writeFile(pdb.topology, final.getPositions(), stream, keepIds=True)
            record.update(status="completed", completed_absolute_step=simulation.currentStep,
                segment_start_after_step=start_step, timed_steps=remaining, elapsed_s=elapsed,
                ns_per_day=remaining * SETTINGS["timestep_ps"] / 1000 / elapsed * 86400 if remaining and elapsed else None,
                final_potential_energy_kj_mol=energy,
                limitation="Technical pilot only; not affinity, activity or convergence evidence" if mode == "pilot" else "Sampling, pose behavior and independent-repeat consistency require analysis")
            summary = {"status": "completed", "run_identity_sha256": digest(identity),
                "completed_absolute_step": simulation.currentStep, "completed_attempt": folder.name,
                "production_ns": identity["production_ns"], "completed_utc": utcnow()}
            # The attempt record is committed first; the run summary is merely
            # an index and must never claim completion before its source does.
            record["finished_utc"] = utcnow()
            atomic_json(folder / "attempt_status.json", record)
            atomic_json(directory / "run_summary.json", summary)
            return {**record, "idempotent_skip": False}
        except BaseException as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            if simulation is not None:
                try:
                    record["last_observed_absolute_step"] = simulation.currentStep
                except Exception:
                    pass
            raise
        finally:
            if segment is not None:
                segment.close()
            if equilibration is not None:
                stream = getattr(equilibration, "_out", None)
                if stream is not None and hasattr(stream, "close") and not stream.closed:
                    stream.close()
            record["finished_utc"] = utcnow()
            atomic_json(folder / "attempt_status.json", record)
            rebuild_segment_manifest(directory, identity)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["pilot", "production", "review"])
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--case", required=True)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--platform", choices=["CUDA", "OpenCL", "CPU", "Reference"], default="CUDA")
    parser.add_argument("--ns", type=float, default=20)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--restart-from-inputs", action="store_true")
    parser.add_argument("--reviewer")
    parser.add_argument("--note")
    parser.add_argument("--pilot-attempt")
    for name in REVIEW_CHECKS:
        parser.add_argument("--" + name.replace("_", "-") + "-checked", action="store_true")
    args = parser.parse_args(argv)
    if args.action == "review":
        if not args.reviewer or not args.note:
            parser.error("review requires --reviewer and --note")
        result = approve_pilot(args.root.resolve(), args.case, args.seed, args.reviewer, args.note,
            {name: getattr(args, name + "_checked") for name in REVIEW_CHECKS}, args.pilot_attempt)
    else:
        result = run(args.case, args.action, args.seed, args.platform, args.ns, args.resume,
            args.root, args.restart_from_inputs)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    main()
