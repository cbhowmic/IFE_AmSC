#!/usr/bin/env python

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import matplotlib.pyplot as plt


class SurrogateMLP(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim=100, num_hidden_layers=3):
        super().__init__()

        layers = []
        in_dim = input_dim

        for _ in range(num_hidden_layers):
            layers.append(nn.Linear(in_dim, hidden_dim))
            layers.append(nn.ReLU())
            in_dim = hidden_dim

        layers.append(nn.Linear(in_dim, output_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test trained RHINO surrogate model."
    )

    parser.add_argument(
        "--model",
        default="artifacts/rhino_surrogate.pt",
        help="Path to trained surrogate artifact.",
    )

    parser.add_argument(
        "--splits",
        default="artifacts/split_indices.json",
        help="Path to saved train/val/test split indices.",
    )

    parser.add_argument(
        "--outdir",
        default="artifacts",
        help="Directory where test results will be saved.",
    )

    return parser.parse_args()


def compute_metrics(y_true, y_pred):
    error = y_pred - y_true

    metrics = {
        "mse": float(np.mean(error ** 2)),
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "mae": float(np.mean(np.abs(error))),
        "max_abs_error": float(np.max(np.abs(error))),
    }

    return metrics


def save_parity_plots(y_true, y_pred, output_columns, outdir):
    n_outputs = len(output_columns)

    fig, axes = plt.subplots(
        nrows=1,
        ncols=n_outputs,
        figsize=(5 * n_outputs, 4),
        dpi=300,
    )

    if n_outputs == 1:
        axes = [axes]

    for i, col in enumerate(output_columns):
        true_i = y_true[:, i]
        pred_i = y_pred[:, i]

        axes[i].scatter(true_i, pred_i, s=15)

        min_val = min(true_i.min(), pred_i.min())
        max_val = max(true_i.max(), pred_i.max())
        axes[i].plot([min_val, max_val], [min_val, max_val])

        axes[i].set_xlabel("True")
        axes[i].set_ylabel("Predicted")
        axes[i].set_title(col)

    fig.tight_layout()
    fig.savefig(outdir / "parity_plots.png")
    plt.close(fig)


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    artifact = torch.load(args.model, map_location="cpu")

    with open(args.splits, "r") as f:
        split_indices = json.load(f)

    features_file = split_indices["features_file"]
    input_columns = split_indices["input_columns"]
    output_columns = split_indices["output_columns"]
    test_idx = split_indices["test_idx"]

    df = pd.read_csv(features_file)

    X = df[input_columns].to_numpy(dtype=np.float32)
    Y = df[output_columns].to_numpy(dtype=np.float32)

    X_test = X[test_idx]
    Y_test = Y[test_idx]

    x_mean = artifact["x_mean"].numpy()
    x_std = artifact["x_std"].numpy()
    y_mean = artifact["y_mean"].numpy()
    y_std = artifact["y_std"].numpy()

    X_test_n = (X_test - x_mean) / x_std

    config = artifact["model_config"]

    model = SurrogateMLP(
        input_dim=config["input_dim"],
        output_dim=config["output_dim"],
        hidden_dim=config["hidden_dim"],
        num_hidden_layers=config["num_hidden_layers"],
    )

    model.load_state_dict(artifact["model_state_dict"])
    model.eval()

    with torch.no_grad():
        X_test_t = torch.tensor(X_test_n, dtype=torch.float32)
        Y_pred_test_n = model(X_test_t).numpy()

    Y_pred_test = Y_pred_test_n * y_std + y_mean

    metrics = compute_metrics(Y_test, Y_pred_test)

    with open(outdir / "metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    pred_df = pd.DataFrame()

    for i, col in enumerate(output_columns):
        pred_df[f"true_{col}"] = Y_test[:, i]
        pred_df[f"pred_{col}"] = Y_pred_test[:, i]
        pred_df[f"abs_error_{col}"] = np.abs(Y_pred_test[:, i] - Y_test[:, i])

    pred_df.to_csv(outdir / "test_predictions.csv", index=False)

    save_parity_plots(Y_test, Y_pred_test, output_columns, outdir)
    
    print("\nTest metrics:")
    print(json.dumps(metrics, indent=2))

    print("\nSaved:")
    print(outdir / "metrics.json")
    print(outdir / "test_predictions.csv")
    print(outdir / "parity_plots.png")


if __name__ == "__main__":
    main()