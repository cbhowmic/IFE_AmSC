#!/usr/bin/env python

import argparse
import json
import copy
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader


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


def split_dataset(X, Y, train_frac=0.80, val_frac=0.15, seed=42):
    rng = np.random.default_rng(seed)

    n = X.shape[0]
    idx = np.arange(n)
    rng.shuffle(idx)

    n_train = int(train_frac * n)
    n_val = int(val_frac * n)

    train_idx = idx[:n_train]
    val_idx = idx[n_train:n_train + n_val]
    test_idx = idx[n_train + n_val:]

    return (
        X[train_idx], Y[train_idx],
        X[val_idx], Y[val_idx],
        X[test_idx], Y[test_idx],
        train_idx, val_idx, test_idx,
    )


def normalize(train, val, test):
    mean = train.mean(axis=0, keepdims=True)
    std = train.std(axis=0, keepdims=True)
    std = np.where(std < 1e-12, 1.0, std)

    return (
        (train - mean) / std,
        (val - mean) / std,
        (test - mean) / std,
        mean,
        std,
    )


def make_loader(X, Y, batch_size=64, shuffle=False):
    X_t = torch.tensor(X, dtype=torch.float32)
    Y_t = torch.tensor(Y, dtype=torch.float32)
    ds = TensorDataset(X_t, Y_t)
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle)


def run_epoch(model, loader, criterion, optimizer=None, device="cpu"):
    model.train() if optimizer else model.eval()

    total_loss = 0.0
    total_samples = 0

    for xb, yb in loader:
        xb = xb.to(device)
        yb = yb.to(device)

        if optimizer:
            optimizer.zero_grad()

        pred = model(xb)
        loss = criterion(pred, yb)

        if optimizer:
            loss.backward()
            optimizer.step()

        total_loss += loss.item() * xb.size(0)
        total_samples += xb.size(0)

    return total_loss / total_samples


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train RHINO surrogate model from feature-layer CSV."
    )

    parser.add_argument(
        "--features",
        default="../feature_extraction/rhino_features.csv",
        help="Path to ML-ready feature CSV from feature extraction layer.",
    )

    parser.add_argument(
        "--inputs",
        nargs="+",
        default=[
            "tritium_burning_rate",
            "burn_fraction",
        ],
        help="Input feature columns.",
    )

    parser.add_argument(
        "--outputs",
        nargs="+",
        default=[
            "tritium_in_isotope_separation",
            "plant_doubling_time_days",
            "minimum_startup_inventory_g",
        ],
        help="Output target columns.",
    )

    parser.add_argument(
        "--outdir",
        default="artifacts",
        help="Directory where model and results will be saved.",
    )

    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--hidden-dim", type=int, default=100)
    parser.add_argument("--hidden-layers", type=int, default=3)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--seed", type=int, default=42)

    return parser.parse_args()

    
def main():
    args = parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.features)

    X = df[args.inputs].to_numpy(dtype=np.float32)
    Y = df[args.outputs].to_numpy(dtype=np.float32)

    X_train, Y_train, X_val, Y_val, X_test, Y_test, train_idx, val_idx, test_idx = split_dataset(
        X, Y, seed=args.seed
    )

    X_train_n, X_val_n, X_test_n, x_mean, x_std = normalize(
        X_train, X_val, X_test
    )
    Y_train_n, Y_val_n, Y_test_n, y_mean, y_std = normalize(
        Y_train, Y_val, Y_test
    )

    train_loader = make_loader(X_train_n, Y_train_n, args.batch_size, shuffle=True)
    val_loader = make_loader(X_val_n, Y_val_n, args.batch_size)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = SurrogateMLP(
        input_dim=X.shape[1],
        output_dim=Y.shape[1],
        hidden_dim=args.hidden_dim,
        num_hidden_layers=args.hidden_layers,
    ).to(device)

    split_indices = {
        "train_idx": train_idx.tolist(),
        "val_idx": val_idx.tolist(),
        "test_idx": test_idx.tolist(),
        "features_file": args.features,
        "input_columns": args.inputs,
        "output_columns": args.outputs,
        "seed": args.seed,
    }

    with open(outdir / "split_indices.json", "w") as f:
        json.dump(split_indices, f, indent=2)

    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    best_val_loss = float("inf")
    best_state = copy.deepcopy(model.state_dict())

    history = {"train_loss": [], "val_loss": []}

    for epoch in range(1, args.epochs + 1):
        train_loss = run_epoch(model, train_loader, criterion, optimizer, device)
        val_loss = run_epoch(model, val_loader, criterion, None, device)

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = copy.deepcopy(model.state_dict())

        if epoch % 20 == 0:
            print(
                f"Epoch {epoch:4d} | "
                f"train_loss={train_loss:.6e} | "
                f"val_loss={val_loss:.6e}"
            )

    model.load_state_dict(best_state)

    artifact = {
        "model_state_dict": model.state_dict(),
        "model_config": {
            "input_dim": X.shape[1],
            "output_dim": Y.shape[1],
            "hidden_dim": args.hidden_dim,
            "num_hidden_layers": args.hidden_layers,
        },
        "input_columns": args.inputs,
        "output_columns": args.outputs,
        "x_mean": torch.tensor(x_mean, dtype=torch.float32),
        "x_std": torch.tensor(x_std, dtype=torch.float32),
        "y_mean": torch.tensor(y_mean, dtype=torch.float32),
        "y_std": torch.tensor(y_std, dtype=torch.float32),
        "best_val_loss": best_val_loss,
    }

    torch.save(artifact, outdir / "rhino_surrogate.pt")

    pd.DataFrame(history).to_csv(outdir / "training_history.csv", index=False)

    print("\nSaved:")
    print(outdir / "rhino_surrogate.pt")
    print(outdir / "training_history.csv")
    print(outdir / "split_indices.json")
    

if __name__ == "__main__":
    main()