# NeuroFEM SpiNNaker2
This repo contains code for the NeuroFEM algorithm implemented on SpiNNaker2.
## repo structure
```
.
└── src
    ├── neurofem                                # python code for neurofem
    │   ├── simulation                          # main NeuroFEM source code
    │   ├── config                              # configuration classes for SpiNNaker2 and NeuroFEM
    │   ├── problems
    │   │   ├── poisson                         # generation of linear systems for Poisson problems
    │   │   └── utils                           # utility functions for problem generation
    │   └── calc                                # mse, error calculation, etc.
    ├── requirements                            # directory for dependencies
    └── requirements.txt                        # dependencies file
```
## install
To set up the environment, first download the dependencies:
```bash
cd src
git clone --branch neurofem_v6 git@gitlab.com:spinnaker2/py-spinnaker2.git requirements/spinnaker2
git clone --branch neurofem_v6 git@gitlab.com:spinnaker2/s2-sim2lab-app.git requirements/spinnaker2/src/spinnaker2/libs
```
Then (create a virtual env and) install the dependencies:
<!-- ToDo: check if libs_install.sh is needed (only on main branch) -->
```bash
python3 -m venv ./venv && source ./venv/bin/activate
pip install -U -r requirements.txt
cp -r requirements/spinnaker2/src/spinnaker2/libs/chip/app-pe/s2app/neurofem_2048/ venv/lib/python3.12/site-packages/spinnaker2/libs/chip/app-pe/s2app/
```
## usage
To run NeuroFEM on SpiNNaker2 use the NeurofemSimulation class:
```python
from neurofem import NeurofemSimulation, NeurofemConfig, SpinnakerConfig
from neurofem.problems import generate_poisson_unitdisk_variable_f

matrix, rhs, _ = generate_poisson_unitdisk_variable_f(nrefs=5)
spinn_config = SpinnakerConfig(s2_ip="192.168.0.X")
neurofem_config = NeurofemConfig(matrix=matrix)
simulation = NeurofemSimulation(spinn_config, neurofem_config, matrix)

solution = simulation.run()
```
