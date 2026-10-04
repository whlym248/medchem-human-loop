# Public-release validation — v0.2.0, 2026-10-04

Executed locally using synthetic files, small chemical graphs and synthetic
geometry only:

- Standard-library environment: **65 passed; 20 explicitly skipped** (15 RDKit
  graph tests and 5 optional analysis tests). The missing-RDKit error test passed;
  skips are not chemistry or geometry validation.
- Existing scientific environment: **85 passed, 0 skipped**, including all graph
  and analysis tests. Dependencies were already installed; no MD was run.
- Actual CLI chain: aromatic toy parent → **4 graph proposals** → parent plus
  proposals on 2 toy targets → **10 evidence rows** → frozen review packet →
  explicitly synthetic decision → identity-checked second round.
- Every docking, ADMET/property and MD measurement in that generated-proposal
  chain remained **missing**. Synthetic case matches did not enable real case
  support. Decisions were demonstrations, not real candidate approval.
- Original evidence assembly remains covered with its four-row fixture, whose
  numerical values are manually invented and explicitly synthetic.
- **No MD, OpenMM Context/Simulation, docking, ADMET inference, model training,
  remote school access or paid compute was executed for these checks.**

Scientific dependency versions are unchanged from the v0.1 record below. The new
HPC probe was tested with mocked GPU responses and synthetic state files,
**not on a real school node**.

## Reproduce the lightweight checks

```console
python -m unittest discover -s tests -v
python scripts/assemble_evidence.py --input examples/synthetic_evidence.json --out outputs/demo
python scripts/workflow_round.py prepare --input examples/synthetic_evidence.json --parent toy_parent --out outputs/round1
```

With RDKit, follow the graph-to-evidence chain in [the script guide](../scripts/README.md).
Use fresh output paths. Decision and next-round examples are in
[review rounds](../docs/review-rounds.md).

CI covers standard-library contracts on Python 3.11/3.12 and optional graph
editing on Python 3.12. The graph job installs only `requirements-chem.txt`;
analysis tests are explicitly skipped there. CI runs no simulation or prediction
service. Local results above are distinct from a remote CI run's own result.

## New coverage

- Case schema, sources, duplicate/nonfinite JSON, lexical retrieval and synthetic isolation.
- Restricted edits, atom mapping, fixed neighborhoods, deduplication, unsupported chemistry and scope propagation.
- Proposal IDs, parent linkage, safe filenames, missing metrics and refusal to import fake scores.
- Rule suggestions, missing support and truthful diffusion unavailability.
- Frozen snapshots, path confinement, escaped HTML, endpoint mismatch, explicit decision acknowledgement, history and approved-parent identity/scope.
- Probe's bounded GPU query, partial/changed-state handling, selected state fields and exclusive new receipt creation.

Windows testing revealed that mocking `subprocess.run` also intercepted an
internal `platform` query in one test. The test now stubs host/platform values;
the full suite then passed. This was a test-isolation failure, not a GPU result.
No dependencies were installed or modified. As before, normal Conda activation
must supply the existing scientific environment's DLL search path.

These checks validate software behavior and bookkeeping, not live HPC
deployment, synthesis, affinity, toxicity, convergence or optimization success.

## Previous release — v0.1.0, 2026-10-03

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
