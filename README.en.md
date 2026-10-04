# MedChem Human Loop

**Traceable computational evidence for local molecular optimization, with human decisions on risks and benefits.**

[中文](README.md) · [Workflow design](docs/workflow-design.md) · [Implementation status](docs/implementation-status.md) · [Development plan](docs/development-plan.md) · [HPC agent handoff guide](docs/hpc/agent-start-here.md)

Author: **whlym**. Workflow conception, computational work, analysis, and code preparation used AI assistance. Original code is released under the [MIT License](LICENSE); citation metadata is in [CITATION.cff](CITATION.cff).

The proposed workflow combines literature, docking, molecular dynamics (MD), and ADMET predictions to formulate testable binding hypotheses, suggest local molecular edits, and compare candidates with explicit uncertainties. A human decides whether to proceed, gather evidence, backtrack, or stop.

**v0.2 implements the first workflow modules**: a local case library and lexical retrieval, bounded edits at a manually selected aromatic site, transparent route rules, proposal-to-evidence conversion, and frozen review packets with static HTML and human decision records. Existing MD execution, quality checks, trajectory analysis, and optional fixed-pose scoring/ADMET-AI reference tools remain available.

The bundled cases are synthetic examples; the real literature corpus is **empty**. Editing is limited to explicit rules rather than general chemistry, route suggestions are not learned, and HTML/CLI review is not a full interactive application. Constrained diffusion, model training, complete computation orchestration, and validated automatic benefit/risk advice remain unimplemented.

## Intended workflow

1. Establish molecular identity, target structure, research goals, and constraints.
2. Retrieve literature or inspect docking poses, distances, and chemical context.
3. Form a testable binding hypothesis and identify a local edit region.
4. Use medicinal-chemistry precedents to propose edits: bounded rule-based editing now exists; locally constrained diffusion remains a design goal.
5. Audit chemistry and docking, review a technical MD pilot, and run tracked production MD.
6. Combine trajectory observations with ADMET and other property predictions.
7. Present several options with provenance, expected benefits, risks, and missing evidence.
8. Record a human decision before another iteration.

This describes the design, not a claim that every stage or connection is automated. See the [status table](docs/implementation-status.md) and [development plan](docs/development-plan.md).

## Lightweight quick start

Python 3.10+ is sufficient. No external packages are required for this example.

```bash
python scripts/workflow_round.py prepare --input examples/synthetic_evidence.json --parent toy_parent --out outputs/round1
python scripts/workflow_round.py status --round outputs/round1
python -m unittest discover -s tests -v
```

Open `outputs/round1/report.html` to inspect **synthetic, manually constructed evidence**, missing values, and comparison limits. Preparation freezes evidence/structure files; the status command checks integrity and decision records. It does not run docking, MD, ADMET inference, or molecular generation. Its values are not scientific results. Real project datasets, trajectories, and prepared MD systems are not included.

Validation covers supplied fields and completeness. Evidence statuses such as `completed_reviewed` are declarations from the input provider; the assembler does not authenticate them against actual trajectories or review reports.

The round output directory must not exist; use another `--out` path when repeating the example. File hashes verify unchanged bytes, not chemical identity or scientific validity. See [review rounds](docs/review-rounds.md) for decision recording and another iteration. Optional RDKit generation is documented in [bounded edits](docs/local-edits.md) and the [script guide](scripts/README.md); it is separate from this standard-library example.

See the [validation record](tests/VALIDATION.md) for checks actually run, environments, coverage, and optional-dependency skips. Lightweight tests and examples start no MD. Skipped chemistry tests are not chemical validation in that environment.

## Tools

| Script | Purpose |
|---|---|
| `case_library.py` | Validate local JSONL cases and retrieve reviewed lexical/tag matches, preserving source and assay context; only synthetic records are bundled. |
| `enumerate_local_edits.py` | Apply four bounded transformations at a selected aromatic C–H site with protected-neighborhood checks; requires RDKit. |
| `route_proposal.py` | Produce transparent route/evidence-gathering suggestions; constrained diffusion remains blocked. |
| `proposals_to_evidence.py` | Export local proposals to evidence inputs with calculation values explicitly missing. |
| `workflow_round.py` | Freeze evidence, produce static HTML, record human decisions, and check the approved parent's identity between rounds. |
| `assemble_evidence.py` | Validate supplied evidence and produce a traceable candidate table. |
| `run_md_v2.py` | Run a prebuilt OpenMM system with identity checks, reviewed pilots, attempts, checkpoints, and segment manifests. |
| `run_array_task_v2.py` | Dispatch a declared MD task using the same execution contract. |
| `check_pilot_quality.py` | Inspect pilot logs and final geometry; produce a technical report without automatically approving it. |
| `analyze_production_md.py` | Analyze retained trajectory segments, RMSD, heavy-atom contacts, and direct N/O hydrogen bonds. |
| `extras/fixed_pose_admet/run_reference.py` | Run Vina/Vinardo fixed-pose scoring and pretrained ADMET-AI inference on two user-prepared inputs. |
| `hpc_probe.py` | Collect read-only operational evidence on an authorized compute host; does not log in, configure resources, or submit jobs. |

Core scripts are in `scripts/`. Optional dependencies are separated into `requirements-chem.txt` for editing, `requirements-md.txt` for MD, and `requirements-analysis.txt` for quality/trajectory analysis. The runner requires compatible user-supplied systems; it does not prepare arbitrary proteins or parameterize ligands. Review the [MD input contract](docs/md-execution.md) and [script documentation](scripts/README.md) before use. A new operator should start with the [HPC agent guide](docs/hpc/agent-start-here.md), which distinguishes historical incidents from facts that require a fresh check and contains no live account credentials.

The optional [fixed-pose/ADMET tool](extras/fixed_pose_admet/README.md) was adapted from a previously executed research check. It requires exactly two prepared inputs, installed software, trusted model weights, and input manifests. It performs `--score_only`, not docking search; it does not train or generate molecules. This public adaptation was checked for syntax and CLI behavior only, without rerunning scoring or inference. Data and weights are excluded. Its outputs need human review and conversion before entering the evidence assembler; it is not a complete workflow connector.

## Interpretation limits

- Docking scores and geometric interactions are evidence for hypotheses, not measured affinity or proof of mechanism.
- A short, single-seed trajectory does not establish convergence, efficacy, or successful optimization.
- ADMET classifier scores are not clinical toxicity incidence rates.
- Only committed trajectory prefixes declared by `segment_manifest.json` may be combined. Binary checkpoints must not be assumed portable across platforms.
- No wet-lab validation or superior candidate claim is supplied with this release.

AI assistance does not remove the need to review applicability, evidence conflicts, provenance, and scientific uncertainty. New case retrieval, editing, and routing use lexical matching, fixed chemical-graph transformations, and an explicit decision table; they neither train nor call a generative model. An iterative loop is not automatically reinforcement learning: this repository contains no learned policy, reward training, or demonstrated “chemist intuition” model.

Next work includes licensed and verified real cases, chemical retrieval, broader edit validation, complete remote computation integration, constrained diffusion, and evaluated risk/benefit recommendations. See the [module-by-module plan](docs/development-plan.md) for acceptance boundaries.

## Authorship and licensing

**whlym** is the author of the workflow concept and contributed computational execution, analysis, and code preparation with AI assistance. Citation is appreciated, using [CITATION.cff](CITATION.cff). This request does not add restrictions to the MIT License. Third-party tools, models, literature, and datasets retain their own licenses. See [AUTHORS.md](AUTHORS.md) and [release provenance](docs/provenance-and-release.md).
