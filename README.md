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
