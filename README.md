# IFE_AmSC

This project provides tools for organizing experimental and simulation data into structured, AI-ready formats using a unified, well-defined and interoperable schemas.

The goal of this project is to standardize heterogeneous scientific datasets to enable reproducible data processing and machine learning workflows.

---

## Environment Setup

### Create the environment

```bash
conda env create -f environment.yml
```

### Activate the environment
```bash
conda activate IFE_AmSC
```
### On NERSC
  
  
## Installation on NERSC

This project uses a conda environment to manage dependencies. To set up the environment:

1. Activate python
```bash
module load python
```

2. Clone this repo
```bash
git clone git@github.com:cbhowmic/IFE_AmSC.git
cd IFE_AmSC
```

3. Create and activate the environment:
```bash
mkdir -p /global/cfs/cdirs/m3239/$(whoami)/sw/perlmutter/IFE_AmSC
conda env create --prefix /global/cfs/cdirs/m3239/$(whoami)/sw/perlmutter/IFE_AmSC -f environment.yml
conda activate /global/cfs/cdirs/m3239/$(whoami)/sw/perlmutter/IFE_AmSC
```

4. Enable the environment in notebooks
```bash
python -m ipykernel install --user --name env --display-name 'IFE AmSC'
```

4. Go to `https://jupyter.nersc.gov/`, select "Login node", open the notebook `... .ipynb` and select the kernel `IFE AmSC`.

## Project Structure
```text
IFE_AmSC/
├── environment.yml
├── README.md
├── Laser/
├── RHINO/
├── IPM/
├── ...

Each module corresponds to a different application area within the project.

Each subdirectory contains application-specific scripts that write simulation or experimental data into AI-ready ADIOS datasets. These datasets are then placed into campaigns that can be shared across users and queried for specific metadata or scientific information.
