# Third-party software and data

Copyright and the MIT license in this repository apply to its original code and documentation, not to third-party dependencies, model weights, molecular structures, force fields, publications, or external datasets. No third-party source library or pretrained weight file is vendored here.

The scientific tools call external libraries. Consult the upstream project and the exact installed release for license terms, version requirements and academic citations:

| Software/resource | Role | Upstream |
| --- | --- | --- |
| OpenMM | Molecular simulation engine and serialized system/state interfaces | https://github.com/openmm/openmm |
| AutoDock Vina | Optional reference tool for fixed-pose Vina/Vinardo scoring | https://github.com/ccsb-scripps/AutoDock-Vina |
| NumPy | Numerical arrays and geometry | https://github.com/numpy/numpy |
| MDTraj | Trajectory reading and analysis | https://github.com/mdtraj/mdtraj |
| Matplotlib | Analysis plots | https://github.com/matplotlib/matplotlib |
| RDKit | Optional molecular identity/chemistry tooling in the wider workflow | https://github.com/rdkit/rdkit |
| Open Force Field Toolkit and NAGL | External ligand preparation and estimated charges; not supplied by this release | https://github.com/openforcefield/openff-toolkit and https://github.com/openforcefield/openff-nagl |
| ADMET-AI | External property prediction referenced in the workflow; no weights or predictions distributed here | https://github.com/swansonk14/admet_ai |
| Chemprop, PyTorch, Lightning and pandas | Optional ADMET reference runtime dependencies | https://github.com/chemprop/chemprop, https://github.com/pytorch/pytorch, https://github.com/Lightning-AI/pytorch-lightning and https://github.com/pandas-dev/pandas |

This list records software relationships, not endorsement or affiliation. The current scripts and installation instructions determine which dependencies are actually needed. The standard-library demonstration does not import the scientific dependencies.

Prepared protein/ligand systems, raw research data, trajectories, model files and papers are deliberately excluded. Obtain your own inputs from sources that permit your intended use and redistribution. Record the input source, license, preparation, exact model/force-field versions and checksums. A parameter file generated using external software is not automatically covered by this repository's license.

For research use, cite both this project (see `CITATION.cff`) and the relevant upstream methods. Citation is a scholarly request, not an extra restriction on the MIT license.
