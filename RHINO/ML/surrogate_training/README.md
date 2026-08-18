# RHINO Surrogate Training

This layer trains and evaluates a multi-output neural-network surrogate using
the CSV produced by `AI_ready_workflow/3_feature_extraction`. It preserves the
feature-table provenance needed to reproduce the test split and writes model,
metric, prediction, and plotting artifacts.

## Relationship to the Feature Layer

The current feature table begins with campaign identity and simulation-time
metadata, followed by declared model inputs and outputs:

```text
archive, datasetid, run_id,
simulation_date, simulation_time, simulation_datetime,
<model inputs>, <model outputs>
```

Identity and datetime columns are not passed to the neural network or included
in normalization. By default, `trainSurrogate.py` reads the `inputs` and
`outputs` keys from the feature layer's `feature_spec.json`; `metadata` entries
are intentionally excluded. The available metadata columns are retained in
`split_indices.json` and copied into `test_predictions.csv`, allowing every
prediction to be traced back to its campaign run.

The trainer rejects missing columns, nonnumeric model columns, NaN or infinite
model values, and datasets too small to create non-empty training, validation,
and test splits. Missing feature values should therefore be resolved in the
feature/campaign layers rather than silently dropped during training.

## Requirements

- Python 3.10 or newer
- NumPy and pandas
- PyTorch
- Matplotlib for evaluation plots
- MLflow when remote tracking or model registration is enabled
- A completed `rhino_features.csv`

The RHINO package and its I/O dependencies are required to rebuild the
upstream BP5 and feature products, but training itself reads the resulting CSV.
Install the training dependencies from the RHINO package root with:

```bash
python -m pip install -e '.[ml]'
```

Use `'.[io,ml]'` when the same environment must also rebuild BP5 products.

## Training Workflow

```text
feature_spec.json + rhino_features.csv
                  |
                  v
Validate selected model columns and finite values
                  |
                  v
Seeded 80% / 15% / 5% train-validation-test split
                  |
                  v
Normalize from training statistics only
                  |
                  v
Train MLP and retain the best validation state
                  |
                  v
Model, normalization, dataset fingerprint, and split metadata
```

The default network has three hidden layers, 100 units per hidden layer, ReLU
activations, and one output node per declared target. It uses Adam and
mean-squared error on normalized targets.

## Files

### `trainSurrogate.py`

Loads and validates the feature CSV, selects input/output columns, creates
seeded data splits, normalizes values, trains the model, and saves:

```text
artifacts/
├── rhino_surrogate.pt
├── training_history.csv
└── split_indices.json
```

`rhino_surrogate.pt` contains model weights, architecture, column names,
normalization statistics, and best validation loss. `split_indices.json`
contains source row indexes, the resolved CSV and feature-spec paths, CSV
SHA-256 fingerprint, row count, metadata columns, random seed, and—when
enabled—the MLflow run, model URI, registered-model name, and version.

### `testSurrogate.py`

Reconstructs the held-out split, verifies the CSV fingerprint and model-column
contract, generates predictions, and saves:

```text
artifacts/
├── metrics.json
├── test_predictions.csv
└── parity_plots.png
```

Metrics include aggregate and per-output MSE, RMSE, MAE, and maximum absolute
error. Per-output metrics are the meaningful comparison when targets have
different units. The prediction CSV includes campaign identity, simulation
datetime, model inputs, original source-row number, true values, predictions,
and absolute errors.

### `mlflow_model.py`

Defines the serving interface registered with MLflow. It accepts a pandas
table with raw, named RHINO input columns, applies the saved input
normalization, evaluates the PyTorch network, reverses output normalization,
and returns a table with named outputs in physical units. This prevents API
clients from having to reproduce training normalization.

## Usage

The default paths are resolved from the script locations, so the commands can
be run from any working directory. The expected default feature table is:

```text
RHINO/AI_ready_workflow/3_feature_extraction/outputs/rhino_features.csv
```

Train using the feature JSON's declared inputs and outputs:

```bash
python RHINO/ML/surrogate_training/trainSurrogate.py
```

Use a different feature CSV or explicitly select model columns:

```bash
python RHINO/ML/surrogate_training/trainSurrogate.py \
  --features /path/to/rhino_features.csv \
  --inputs tritium_burning_rate burn_fraction \
  --outputs plant_doubling_time_days \
            minimum_startup_inventory_g \
            tritium_in_isotope_separation
```

Other useful training controls include `--epochs`, `--batch-size`,
`--hidden-dim`, `--hidden-layers`, `--lr`, `--seed`, `--feature-spec`, and
`--outdir`. Run the training command with `--help` for the complete interface.

Evaluate the saved model on its held-out rows:

```bash
python RHINO/ML/surrogate_training/testSurrogate.py
```

If the exact training CSV was moved without being changed, provide its new
location. The SHA-256 fingerprint must still match:

```bash
python RHINO/ML/surrogate_training/testSurrogate.py \
  --features /new/location/rhino_features.csv
```

Existing artifacts created by the earlier scripts do not contain the new CSV
fingerprint or metadata-column declarations. Retrain after rebuilding the
feature CSV to obtain the complete provenance and prediction schema.

## MLflow Tracking and Model Registry

MLflow is optional: commands without `--mlflow` retain the local-only behavior
described above. When enabled, one training execution creates one MLflow run
and logs:

- Network, optimizer, split, column, and seed parameters
- Training and validation loss at every epoch
- Best validation loss
- Dataset and feature-spec SHA-256 provenance
- Native checkpoint, training history, split metadata, and feature spec
- A signature-bearing serving model with normalization included

By default, the serving model is registered as `rhino-surrogate`. Pass
`--skip-model-registration` to log it only as a run model, or use
`--registered-model-name` to choose another registry name. Model Registry
requires a database-backed MLflow store.

### Start a small tracking server

For personal testing, use a persistent SQLite database and artifact directory:

```bash
mlflow server \
  --backend-store-uri sqlite:////path/to/mlflow/mlflow.db \
  --artifacts-destination file:///path/to/mlflow/artifacts \
  --host 127.0.0.1 \
  --port 5000
```

Use PostgreSQL and managed/shared artifact storage for a multi-user production
service. Put remotely accessible servers behind HTTPS and authentication; do
not expose an unauthenticated tracking server publicly.

### Train and log

Point the client at the server and enable logging explicitly:

```bash
export MLFLOW_TRACKING_URI=http://127.0.0.1:5000

python RHINO/ML/surrogate_training/trainSurrogate.py \
  --mlflow \
  --mlflow-experiment rhino-surrogate \
  --mlflow-run-name baseline-2026-04-29
```

`--mlflow-tracking-uri` may be used instead of the environment variable. If
MLflow logging fails, the MLflow run is marked failed by its run context.

### Evaluate and attach results

Evaluation reads the originating run ID from `split_indices.json`, reopens that
run, and logs aggregate/per-output test metrics and evaluation artifacts:

```bash
python RHINO/ML/surrogate_training/testSurrogate.py --mlflow
```

After reviewing the metrics, an alias can be assigned explicitly:

```bash
python RHINO/ML/surrogate_training/testSurrogate.py \
  --mlflow \
  --model-alias candidate
```

Use an alias such as `champion` only after the version satisfies the project's
validation policy. The alias and registered version are not created when
training uses `--skip-model-registration`.

### Serve a registered version

Registration stores and versions the model but does not create a continuously
running inference endpoint. For local validation, serve a selected alias with:

```bash
mlflow models serve \
  --model-uri 'models:/rhino-surrogate@candidate' \
  --port 8080
```

The endpoint accepts raw physical inputs by name:

```bash
curl http://127.0.0.1:8080/invocations \
  -H 'Content-Type: application/json' \
  -d '{
    "dataframe_split": {
      "columns": ["tritium_burning_rate", "burn_fraction"],
      "data": [[82.4, 0.035]]
    }
  }'
```
