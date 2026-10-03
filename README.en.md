# MedChem Human Loop

**Traceable computational evidence for local molecular optimization, with human decisions on risks and benefits.**

[中文](README.md) · [Workflow design](docs/workflow-design.md) · [Implementation status](docs/implementation-status.md)

Author: **whlym**. Workflow conception, computational work, analysis, and code preparation used AI assistance. Original code is released under the [MIT License](LICENSE); citation metadata is in [CITATION.cff](CITATION.cff).

The proposed workflow combines literature, docking, molecular dynamics (MD), and ADMET predictions to formulate testable binding hypotheses, suggest local molecular edits, and compare candidates with explicit uncertainties. A human decides whether to proceed, gather evidence, backtrack, or stop.

This release includes the **workflow design, MD execution and review utilities, trajectory analysis, an optional fixed-pose scoring and pretrained ADMET-AI reference interface, and a runnable synthetic evidence example**. It is not a trained reinforcement-learning system or a complete autonomous drug-discovery platform. A medicinal-chemistry case library, constrained local diffusion, and automatic selection between edit methods remain planned.

## Intended workflow

1. Establish molecular identity, target structure, research goals, and constraints.
2. Retrieve literature or inspect docking poses, distances, and chemical context.
3. Form a testable binding hypothesis and identify a local edit region.
4. Use medicinal-chemistry precedents to propose functional-group edits or locally constrained diffusion.
5. Audit chemistry and docking, review a technical MD pilot, and run tracked production MD.
6. Combine trajectory observations with ADMET and other property predictions.
7. Present several options with provenance, expected benefits, risks, and missing evidence.
8. Record a human decision before another iteration.

This describes the design, not a claim that every stage is implemented. See the [status table](docs/implementation-status.md).

## Lightweight quick start

Python 3.10+ is sufficient. No external packages are required for this example.

```bash
python scripts/assemble_evidence.py --input examples/synthetic_evidence.json --out outputs/demo
python -m unittest discover -s tests -v
```

The example validates and exports **synthetic, manually constructed evidence**. It does not run docking, MD, ADMET inference, or molecular generation. Its values are not scientific results. Real project datasets, trajectories, and prepared MD systems are not included.

Validation covers supplied fields and completeness. Evidence statuses such as `completed_reviewed` are declarations from the input provider; the assembler does not authenticate them against actual trajectories or review reports.

The output directory must be new or empty; use another `--out` path when repeating the example. File hashes verify unchanged bytes, not chemical identity or scientific validity.

Release checks passed 17 lightweight tests and verified the four-row example using synthetic inputs, without starting MD. See the [validation record](tests/VALIDATION.md) for environments and coverage.

## Tools

| Script | Purpose |
|---|---|
| `assemble_evidence.py` | Validate supplied evidence and produce a traceable candidate table. |
| `run_md_v2.py` | Run a prebuilt OpenMM system with identity checks, reviewed pilots, attempts, checkpoints, and segment manifests. |
| `run_array_task_v2.py` | Dispatch a declared MD task using the same execution contract. |
| `check_pilot_quality.py` | Inspect pilot logs and final geometry; produce a technical report without automatically approving it. |
| `analyze_production_md.py` | Analyze retained trajectory segments, RMSD, heavy-atom contacts, and direct N/O hydrogen bonds. |
| `extras/fixed_pose_admet/run_reference.py` | Run Vina/Vinardo fixed-pose scoring and pretrained ADMET-AI inference on two user-prepared inputs. |

Core scripts are in `scripts/`. MD dependencies are separated into `requirements-md.txt`; analysis dependencies into `requirements-analysis.txt`. The runner requires compatible user-supplied systems; it does not prepare arbitrary proteins or parameterize ligands. Review the [MD input contract](docs/md-execution.md) and [script documentation](scripts/README.md) before use.

The optional [fixed-pose/ADMET tool](extras/fixed_pose_admet/README.md) was adapted from a previously executed research check. It requires exactly two prepared inputs, installed software, trusted model weights, and input manifests. It performs `--score_only`, not docking search; it does not train or generate molecules. This public adaptation was checked for syntax and CLI behavior only, without rerunning scoring or inference. Data and weights are excluded. Its outputs need human review and conversion before entering the evidence assembler; it is not a complete workflow connector.

## Interpretation limits

- Docking scores and geometric interactions are evidence for hypotheses, not measured affinity or proof of mechanism.
- A short, single-seed trajectory does not establish convergence, efficacy, or successful optimization.
- ADMET classifier scores are not clinical toxicity incidence rates.
- Only committed trajectory prefixes declared by `segment_manifest.json` may be combined. Binary checkpoints must not be assumed portable across platforms.
- No wet-lab validation or superior candidate claim is supplied with this release.

AI assistance does not remove the need to review applicability, evidence conflicts, provenance, and scientific uncertainty. An iterative loop is not automatically reinforcement learning: this repository contains no learned policy, reward training, or demonstrated “chemist intuition” model.

## Authorship and licensing

**whlym** is the author of the workflow concept and contributed computational execution, analysis, and code preparation with AI assistance. Citation is appreciated, using [CITATION.cff](CITATION.cff). This request does not add restrictions to the MIT License. Third-party tools, models, literature, and datasets retain their own licenses. See [AUTHORS.md](AUTHORS.md) and [release provenance](docs/provenance-and-release.md).
