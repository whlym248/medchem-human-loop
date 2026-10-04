# Code and input contracts

These are research scripts, not an installed end-to-end agent. No prepared
molecular systems, trajectories, proprietary data or model weights are bundled.
Self-contained standard-library examples cover **synthetic evidence assembly,
case retrieval, route suggestions and human review rounds**. Graph enumeration
additionally requires RDKit; these examples do not run predictions or MD.

## Evidence demo (Python 3.10+, standard library)

From the repository root:

```console
python scripts/assemble_evidence.py --input examples/synthetic_evidence.json --out outputs/demo
python -m unittest discover -s tests -v
```

The output directory must be new or empty. The demo writes four candidate-target
rows to `candidates.csv`, `results.json`, and a SHA256 `execution_receipt.json`.
Every numerical value in the example is invented. It runs no docking, inference,
simulation, generation or ranking. It preserves missing values and provisional
MD labels and checks file bytes, schema references and finite numbers. It does
not verify that supplied scientific claims are true or that a SMILES matches a
structure file. Use the review-round tool below to record explicit decisions.

## Local workflow modules

| Script | Behavior | Details |
| --- | --- | --- |
| `case_library.py` | Validate a local JSONL case library; retrieve reviewed records by literal tokens, tags and target | [Case library](../docs/case-library.md); no bundled real corpus |
| `enumerate_local_edits.py` | Enumerate four restricted single-site aromatic edits with RDKit | [Local edits](../docs/local-edits.md); no conformers, synthesis or scoring |
| `proposals_to_evidence.py` | Convert supplied proposals into structure files and evidence JSON with every metric missing | [Local edits](../docs/local-edits.md); hashes are not chemical authentication |
| `route_proposal.py` | Explain a rule-based suggestion or why evidence/backend is missing | [Route rules](../docs/route-selection.md); no learned policy |
| `workflow_round.py` | Freeze evidence, display source/uncertainty, record human decisions and link a next round | [Review rounds](../docs/review-rounds.md); no job scheduling |
| `hpc_probe.py` | Read explicit disk, GPU and optional state information; optionally create a new receipt | [HPC agent entry](../docs/hpc/agent-start-here.md); does not certify a running task |

Standard-library round example:

```console
python scripts/workflow_round.py prepare --input examples/synthetic_evidence.json --parent toy_parent --out outputs/round1
python scripts/workflow_round.py status --round outputs/round1
```

With the optional graph environment, a complete **software demonstration** is:

```console
python scripts/enumerate_local_edits.py --input examples/local_edit_request.json --out outputs/local_edits
python scripts/proposals_to_evidence.py --proposals outputs/local_edits/local_edits.json --target toy_target_A --target toy_target_B --scope synthetic_example --out outputs/proposals
python scripts/workflow_round.py prepare --input outputs/proposals/input.json --parent parent --cases examples/cases/synthetic_cases.jsonl --query aromatic --route-request examples/route_request.json --out outputs/proposal_round
```

Open `outputs/proposal_round/report.html`. The four proposals and their parent
have **no docking, ADMET or MD measurements**. The report preserves missing
values; producing valid graphs and a review packet does not establish useful
chemistry. Choose fresh directories when repeating the demo. Real decisions
require an actual reviewer; a demonstration must remain labeled synthetic.

## Optional dependencies

```console
# Only for optional local graph editing, without MD:
python -m pip install -r requirements-chem.txt
# Only on a compute host authorized for MD:
python -m pip install -r requirements-md.txt
# For read-only pilot/trajectory analysis and optional geometry tests:
python -m pip install -r requirements-analysis.txt
```

Use a separate scientific environment and select packages/drivers supported by
the compute host. The files specify compatibility ranges, not a reproducibility
lock. The archived research environment used OpenMM 8.6.1; binary checkpoints are
not portable between arbitrary OpenMM versions, hardware or platforms. Use the
same verified environment for an existing run, never edit its recorded identity
to force acceptance. This release does not install or configure CUDA.

## Prepared-system contract

The runner accepts **prebuilt OpenMM systems**, not an arbitrary PDB or SMILES.
Preparation/force-field assignment is outside this repository. A user-supplied
bundle must have:

```text
bundle/
  systems/<case>/
    system.xml
    complex.cif
    positions_nm.npy
    metadata.json
  inputs/charged_ligands/<case>.json  # needed by trajectory analysis
  package_manifest.json             # optional frozen input hash manifest
  tasks.json                        # only for the array-task wrapper
```

`<case>` is a simple directory identifier, without path separators.

- `system.xml`: OpenMM-serialized system, one particle per atom in the topology,
  exactly one `MonteCarloBarostat`, and already-defined restraint force global
  parameters `k_protein` and `k_ligand`. The schedule sets these parameters to
  force constants in kJ mol⁻¹ nm⁻², ending at zero before production. Their force
  definitions must actually implement the intended position restraints; the
  runner cannot establish that from the names alone. The system's barostat must
  already have the intended pressure/temperature (the supplied protocol intends
  1 bar / 310 K); the runner changes its frequency and random seed, not those
  thermodynamic targets. Review these during preparation.
- `complex.cif`: OpenMM-readable topology with explicit hydrogens, matching atom
  order and valid periodic box vectors.
- `positions_nm.npy`: finite `(n_particles, 3)` float array, in nanometers, in the
  same order; loaded with pickle disabled.
- `metadata.json`: at least `case_id`, `system_sha256`, and, for quality/analysis,
  `ligand_atom_indices` (global zero-based indices), `temperature_K` (310 for the
  supplied protocol). Record preparation provenance and parameter sources too.
- `inputs/charged_ligands/<case>.json`: explicit chemical graph, with `atoms`
  entries containing unique `name`, `atomic_number`, `formal_charge`,
  `is_aromatic`; `bonds` entries contain zero-based `atom1`, `atom2`, `bond_order`
  (1/2/3 for nonaromatic bonds) and `is_aromatic`. Atom names/elements must map to
  the ligand atoms in the topology. All ligand atoms, including H, are required.
  This graph is used to assign RDKit feature roles; a filename containing
  `charged_ligands` does not imply this analysis calculates partial charges.
- Optional `package_manifest.json`: list (or `{"files": [...]}`) of
  `{"path": "systems/<case>/system.xml", "sha256": "..."}` records covering all
  four system input files. Without it, per-input SHA and the metadata system SHA
  still bind runs, but no independent frozen-package comparison is possible.
- Optional `tasks.json`: `{"pilot": [{"case_id": "...", "seed": 11}],
  "production": [{"case_id": "...", "seed": 11, "ns": 20}]}`. The wrapper runs
  exactly one zero-based index; it is not a scheduler and allocates no resources.

Parameters are recorded in `SETTINGS`. Temperature 310 K, a 2 fs timestep,
restraint schedule and minimization limits are inherited protocol choices, not
universal defaults for every biomolecular system. Verify constraints, masses,
force fields, protonation, initial pose and thermodynamic settings beforehand.

## MD execution and review (on an authorized compute host)

The following are templates for users who have independently prepared and
reviewed a compatible bundle. **They launch MD; they are not the demo or tests.**

```console
python scripts/run_md_v2.py pilot --root /path/to/bundle --case CASE --seed 11 --platform CUDA
python scripts/check_pilot_quality.py --bundle-root /path/to/bundle --case CASE --attempt /path/to/pilot/attempt_0001 --out /path/to/pilot_quality.json
```

Read the full quality report and underlying logs/geometry. A successful exit or
fast throughput is not scientific approval. Only after an actual technical
review should the reviewer record the decision:

```console
python scripts/run_md_v2.py review --root /path/to/bundle --case CASE --seed 11 --reviewer REVIEWER --note "Specific reviewed evidence and limitations" --finite-energy-checked --temperature-density-checked --geometry-checked --throughput-checked
python scripts/run_md_v2.py production --root /path/to/bundle --case CASE --seed 29 --platform CUDA --ns 20
```

Review flags are an attestation, not a scoring algorithm. The quality script
never writes approval. Production verifies reviewed pilot completion, input
hashes, settings, platform and OpenMM version. Outputs reside under
`runs_v2/<mode>/<case>/seed_<seed>/`. Each restart creates a new immutable attempt;
use `--resume` only after inspecting the failure. Without a valid checkpoint,
`--resume --restart-from-inputs` is a deliberate new attempt from input, not
continuation. The second flag is rejected when a valid checkpoint exists.

Two checkpoint slots bind bytes and metadata by SHA256. OS advisory locks guard
concurrent writers. A final summary is committed only after its completed attempt
record. Never concatenate full DCD files blindly: `segment_manifest.json` lists
retained committed prefixes and excludes superseded restart tails.

## Read-only analysis

```console
python scripts/analyze_production_md.py --bundle-root /path/to/bundle --case CASE --run-dir /path/to/bundle/runs_v2/production/CASE/seed_29 --out-dir /path/to/new-analysis --highlight-residue 100
```

The output directory must be new/empty. Analysis uses the valid segment manifest,
PBC molecular reconstruction and a protein N/CA/C fit. Defaults: heavy contacts
≤4 Å; direct N/O hydrogen bonds require D–A≤3.5 Å, H–A≤2.5 Å and D–H–A≥135° in
the same frame. `--highlight-residue` is repeatable and refers to topology residue
IDs, not canonical sequence positions; identical residue numbers in different
chains are all selected. All residue contact occupancies remain in the full
table. Hydrogen-bond occupancy is **per donor–acceptor pair**, not residue-union
occupancy; do not sum pair values to infer a residue probability.

If no supported N/O donor–acceptor groups can be enumerated, hydrogen-bond
tables contain headers only, the number of occupied pairs is zero, and the
hydrogen-bond plot states that no qualified direct N/O hydrogen bonds were found.
Contact and RMSD analysis still proceed; absence of supported groups does not
imply absence of all possible chemical interactions.

`--allow-incomplete` explicitly labels provisional snapshots and requires final
reanalyzing after the run finishes. `--legacy-engineering` is restricted to
completed matching legacy pilots. RMSD uses fixed atom identities without
symmetry correction. Minimum-image geometry assumes an orthogonal or suitably
reduced triclinic box. Only direct N/O hydrogen bonds are modeled; water bridges,
sulfur bonding, desolvation and free energies are outside scope. A single short
replicate cannot establish convergence, affinity, activity or optimization.

## Validation scope

Core tests use fake checkpoint bytes, synthetic ledger entries and synthetic
evidence. Optional tests use tiny arrays and CSVs for PBC reconstruction,
minimum-image distance, NaN detection, and restart-prefix rejection. They create
no OpenMM `Context`/`Simulation` and generate no MD trajectories. The old private
end-to-end test depended on archived molecular data and is not included. This
public release therefore does not claim a fresh full MD integration validation.

See `SOURCE_PROVENANCE.json` for original source digests and public-release
changes. Third-party libraries retain their own licenses; their implementation,
parameters and pretrained weights are not included here.
