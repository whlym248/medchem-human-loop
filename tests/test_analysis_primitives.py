"""Optional synthetic geometry/ledger checks; no trajectories or MD engine runs.

Install requirements-analysis.txt to execute these; the standard-library CI
reports an explicit skip otherwise. OpenMM is imported as a file-format library;
no Context, Integrator or Simulation is instantiated.
"""
import itertools
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import numpy as np
    import analyze_production_md as analysis
    from check_pilot_quality import MinimumImage, csv_summary, json_safe
    OPTIONAL_AVAILABLE = True
except ModuleNotFoundError:
    OPTIONAL_AVAILABLE = False


@unittest.skipUnless(OPTIONAL_AVAILABLE, "Optional analysis dependencies not installed")
class AnalysisPrimitiveTests(unittest.TestCase):
    def test_no_hbond_candidates_return_empty_geometry_and_occupancy(self):
        atoms=[SimpleNamespace(element=SimpleNamespace(atomic_number=6),residue=SimpleNamespace(name="ALA")),
               SimpleNamespace(element=SimpleNamespace(atomic_number=6),residue=SimpleNamespace(name="LIG"))]
        groups,triplets,starts,triplet_groups=analysis.hydrogen_bonds(atoms,np.array([0,1]),{0:0,1:1},[set(),set()],[1],set(),set())
        self.assertEqual(groups,[])
        self.assertEqual(triplets.shape,(0,3))
        class FrameCountOnly:
            def __len__(self):return 3
        with patch.object(analysis.md,"compute_distances",side_effect=AssertionError("Empty pairs must not reach MDTraj")), patch.object(analysis.md,"compute_angles",side_effect=AssertionError("Empty triples must not reach MDTraj")):
            da,ha,angles,qualified,group_qualified,group_angles=analysis.hbond_geometry(FrameCountOnly(),np.empty((0,2),dtype=int),triplets,starts,triplet_groups,SimpleNamespace())
        for values in (da,ha,angles,qualified,group_qualified,group_angles):
            self.assertEqual(values.shape,(3,0))
        counts=group_qualified.sum(axis=0)
        self.assertEqual(counts.tolist(),[])
        self.assertEqual(int((counts>0).sum()),0)
        self.assertEqual(da.min(axis=0).tolist(),[])
        self.assertEqual(group_angles.max(axis=0).tolist(),[])
        for frame in range(3):
            self.assertEqual(np.flatnonzero(da[frame,triplet_groups]<=5.0).tolist(),[])
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp)/"hydrogen_bond_occupancy.csv"
            analysis.write_table(output,[],fields=["donor","acceptor","occupancy"])
            self.assertEqual(output.read_text(encoding="utf-8-sig").strip(),"donor,acceptor,occupancy")

    def test_minimum_image_matches_brute_force(self):
        for box in (np.diag([2., 3., 4.]),
                    np.array([[2., 0., 0.], [.7, 2.5, 0.], [-.5, .6, 3.]])):
            delta = np.array([[.1, .2, .3], [1.5, -1.4, 1.6], [-2.3, 3.8, 4.9]])
            shifts = np.array(list(itertools.product(range(-3, 4), repeat=3))) @ box
            brute = np.min(np.linalg.norm(delta[:, None, :] - shifts[None, :, :], axis=-1), axis=1)
            for function in (MinimumImage(box), lambda x: analysis.mic(x, box)):
                self.assertTrue(np.allclose(np.linalg.norm(function(delta), axis=1), brute))
                self.assertTrue(np.allclose(np.linalg.norm(function(delta + np.array([2, -1, 1]) @ box), axis=1), brute))

    def test_whole_molecule_reconstruction(self):
        reference = np.array([[0, 0, 0], [.1, 0, 0], [.1, .1, 0], [.3, .3, .3], [.4, .3, .3]])
        _, trees = analysis.make_components(5, [(0, 1), (1, 2), (3, 4)])
        box = np.diag([2., 2., 2.])
        shifted = reference.copy()
        shifted[1] += box[0]
        shifted[3:] -= box[1]
        restored = analysis.whole_coordinates(shifted, box, trees, np.array([0, 1, 2]), np.array([3, 4]), np.array([0, 1, 2]), reference)
        self.assertTrue(np.allclose(restored, reference))

    def test_csv_nonfinite_is_not_silently_valid(self):
        with tempfile.TemporaryDirectory() as temp:
            fixture = Path(temp) / "state.csv"
            fixture.write_text('#"Step","Time (ps)","Potential Energy (kJ/mole)","Temperature (K)","Density (g/mL)","Speed (ns/day)"\n500,1,-1,310,1,--\n1000,2,nan,311,1.01,50\n', encoding="utf-8")
            parsed = csv_summary(fixture)
        self.assertEqual(len(parsed["invalid_or_nonfinite"]), 1)
        self.assertEqual(parsed["invalid_or_nonfinite"][0]["problem"], "nonfinite")
        self.assertEqual(len(parsed["expected_speed_placeholders"]), 1)
        self.assertEqual(parsed["columns"]["Temperature (K)"]["mean"], 310.5)
        self.assertEqual(json_safe({"bad": float("nan")}), {"bad": "nan"})

    def test_segment_prefix_and_duplicate_rejection(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "systems" / "toy"
            source.mkdir(parents=True)
            for name in ("system.xml", "complex.cif", "positions_nm.npy", "metadata.json"):
                # Byte fixtures only: never parsed as molecular systems.
                (source / name).write_text("synthetic hash fixture " + name, encoding="utf-8")
            run = root / "run"
            run.mkdir()
            identity = {"case_id": "toy", "input_sha256": {n: analysis.sha(source / n) for n in ("system.xml", "complex.cif", "positions_nm.npy", "metadata.json")},
                        "mode": "production", "seed": 11, "settings": {"timestep_ps": .002}, "equilibration_steps": 100, "target_absolute_step": 500}
            digest = analysis.digest(identity)
            (run / "run_identity.json").write_text(json.dumps(identity), encoding="utf-8")
            (run / "run_summary.json").write_text(json.dumps({"status": "completed", "run_identity_sha256": digest, "completed_absolute_step": 500}), encoding="utf-8")
            segments = []
            for index, steps, start, cap in ((1, [200, 300, 400, 500], 100, 300), (2, [400, 500], 300, 500)):
                folder = run / f"attempt_{index:04d}"
                folder.mkdir()
                lines = ['#"Step","Time (ps)"'] + [f"{s},{s * .002}" for s in steps]
                (folder / "state.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
                segments.append({"attempt": folder.name, "run_identity_sha256": digest,
                                 "state_csv": folder.name + "/state.csv", "trajectory": folder.name + "/unused.dcd",
                                 "include_frame_count": 2, "paired_frames_committed": len(steps), "include_steps_gt": start,
                                 "include_steps_le": cap, "report_interval_steps": 100})
            manifest = {"schema": 2, "run_identity_sha256": digest, "segments": segments}
            path = run / "segment_manifest.json"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            args = SimpleNamespace(case="toy", allow_incomplete=False, legacy_engineering=False)
            records, details = analysis.load_segments(root, run, source, args)
            self.assertEqual(details["retained_frames"], 4)
            self.assertEqual(np.concatenate([r["steps"] for r in records]).tolist(), [200, 300, 400, 500])
            segments[1].update(state_csv="attempt_0001/state.csv", include_steps_gt=100, include_steps_le=300)
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate or backward"):
                analysis.load_segments(root, run, source, args)


if __name__ == "__main__":
    unittest.main()
