# RHINO

Reduced Hydrogen INventory Optimization (RHINO) is a discrete-time, processing-time-based model of the fusion fuel cycle, used to study tritium inventory and startup-inventory requirements for fusion power plants.

This AI-ready workflow for RHINO provides tools for transforming RHINO simulation outputs (recorded in .pkl format) into
AI-ready data, suitable for downstream tasks.

The current workflow standardizes raw outputs, organizes the resulting datasets
into searchable campaigns, extracts ML-ready features, and uses those features
for surrogate modeling. A local demo connects surrogate inference with archived
simulation data and an interactive visualization.

Surrogate training and evaluation optionally integrate with MLflow for
experiment tracking, artifact logging, and model registration. See the
[surrogate training guide](ML/surrogate_training/README.md) for setup and usage.

## Component Guides

- [Shim and data standardization](AI_ready_workflow/1_shim/README.md)
- [Campaign management](AI_ready_workflow/2_campaign/README.md)
- [Feature extraction](AI_ready_workflow/3_feature_extraction/README.md)
- [Surrogate training](ML/surrogate_training/README.md)
- [Surrogate inference and visualization demo](demo/surrogate/README.md)
- [Repository-wide unified schema outline](../schema.md)
- [Schema design for RHINO](AI_ready_workflow/1_shim/schema.md)

## Packaging Status

The Python package under `src/` and its command-line interfaces are under active
development. Use the repository-level [`environment.yml`](../environment.yml)
for the current shared development environment.
