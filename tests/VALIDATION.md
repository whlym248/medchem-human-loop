# Public-release validation — 2026-10-03

Executed during release preparation, using only synthetic fixtures:

- Standard-library environment: **12 passed; 5 optional analysis tests skipped**
  because the analysis dependencies were not installed.
- Existing scientific environment: **17 passed, 0 skipped**. Versions: NumPy
  2.4.6, OpenMM 8.6.1, MDTraj 1.11.1, Matplotlib 3.11.1, RDKit 2026.03.6.
- Evidence CLI: **4 candidate-target rows**, **2 structure-file SHA256 checks**,
  output CSV/JSON digests recorded and checked. Missing MD and docking values
  remained missing. The fixture explicitly marks all numerical data synthetic.
- **No MD, OpenMM Context/Simulation, docking, ADMET inference or model training
  was executed.** No real trajectory or external resource was accessed by tests.

Commands from the repository root:

```console
python -m unittest discover -s tests -v
python scripts/assemble_evidence.py --input examples/synthetic_evidence.json --out outputs/demo
```

The tests cover duplicate/missing evidence pairs, nonfinite values, structure
file tampering and path escape, truthful missing/provisional status retention,
output preservation, superseded restart tails, corrupt-checkpoint fallback,
identity mismatch, missing pilot review, invalid durations, PBC reconstruction,
minimum-image distances against brute force, and nonfinite CSV diagnostics.
An additional zero-N/O-candidate test checks empty geometry arrays, zero occupied
pair count and a header-only occupancy table without invoking MDTraj geometry
functions or creating a simulation.

On Windows, starting a Conda interpreter by its full path without activating
the environment initially caused a delayed numerical-library DLL failure. The
tests passed after its existing `Library/bin` directory was placed on the test
process search path. No dependencies were installed or changed. Activate the
scientific environment normally before using the optional scripts.

These checks validate software primitives and bookkeeping. The public release
does not include a complete prepared molecular system or a new end-to-end MD
execution test, and these results are not evidence of chemical improvement.
