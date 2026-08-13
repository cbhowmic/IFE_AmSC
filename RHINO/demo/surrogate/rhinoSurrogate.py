#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpmd_api as io
import torch 
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader
from pathlib import Path
import glob
import copy


# In[2]:


def _get_first_snapshot(series):
    """Return the first snapshot from an openPMD RHINO series.

    Parameters
    ----------
    series : openPMD.Series
        OpenPMD RHINO series object.

    Returns
    -------
    object
        First snapshot in the series.

    Raises
    ------
    ValueError
        If the series does not contain any snapshots.
    """
    snapshots = series.snapshots()
    if not snapshots:
        raise ValueError("The series does not contain any snapshots.")
    return snapshots[0]

def _validate_species(it, species):
    """Validate that a species exists in a snapshot and is not ``'Times'``.

    Parameters
    ----------
    it : object
        Snapshot object from an openPMD RHINO series.
    species : str
        Species name to validate.

    Returns
    -------
    list[str]
        List of available species names.

    Raises
    ------
    ValueError
        If the species is invalid.
    """
    avail_species = [sp for sp in it.particles if sp != "Times"]
    if species == "Times" or species not in avail_species:
        raise ValueError(
            f"Invalid species {species!r}. Available species are: {avail_species}"
        )
    return avail_species

def _validate_subsystem(it, species, subsystem):
    """Validate that a subsystem exists for a given species.

    Parameters
    ----------
    it : object
        Snapshot object from an openPMD RHINO series.
    species : str
        Species name whose subsystems are being queried.
    subsystem : str
        Subsystem name to validate.

    Returns
    -------
    list[str]
        List of available subsystem names for the species.

    Raises
    ------
    ValueError
        If the subsystem is invalid.
    """
    avail_subsystems = [sub for sub in it.particles[species]["subsystems"]]
    if subsystem not in avail_subsystems:
        raise ValueError(
            f"Invalid subsystem {subsystem!r}. "
            f"Available subsystems for species {species!r} are: {avail_subsystems}"
        )
    return avail_subsystems

def print_summary(series):
    """Print a human-readable summary of an openPMD RHINO series.

    The summary includes:
    - series-level metadata
    - metadata for the first snapshot
    - available particle records and record components in that snapshot

    Parameters
    ----------
    series : openPMD.Series
        OpenPMD RHINO series object to inspect.

    Returns
    -------
    None
        This function prints information to standard output.

    Raises
    ------
    ValueError
        If the series does not contain any snapshots.
    """
    it = _get_first_snapshot(series)

    print("===============")
    print("Series metadata")
    print("===============")

    for a in series.attributes:
        print(f"{a}={series.get_attribute(a)}")

    print("\t==================")
    print("\tIteration metadata")
    print("\t==================")
    for a in it.attributes:
        print(f"\t{a}={it.get_attribute(a):<20}")

    print("\t\t==============")
    print("\t\tAvailable data")
    print("\t\t==============")
    for p in it.particles:
        print(f"\t\t{p}")
        for r in it.particles[p]:
            print(f"\t\t\t{r}")
            for c in it.particles[p][r]:
                print(f"\t\t\t\t{c}")

def get_times(series):
    """Return the time array stored in an openPMD RHINO series.

    This function reads the scalar ``Times/data`` record from the first snapshot
    in the series and returns it as a NumPy array. Time values are expected to
    be in days.

    Parameters
    ----------
    series : openPMD.Series
        OpenPMD RHINO series object containing time data.

    Returns
    -------
    numpy.ndarray
        One-dimensional array of time values in days.

    Raises
    ------
    ValueError
        If the series does not contain any snapshots.
    """
    it = _get_first_snapshot(series)
    times = it.particles["Times"]["data"][io.Record_Component.SCALAR].load_chunk()
    return times

def get_steady_state(series, species="Tritium", subsystem=None):
    """Return the steady-state mass inventory for a species.

    This function loads the ``mass_steady`` data for the requested species from
    the first snapshot in the series. If a subsystem name is provided, only the
    value for that subsystem is returned.

    Parameters
    ----------
    series : openPMD.Series
        OpenPMD RHINO series object containing steady-state mass data.
    species : str, optional
        Name of the species to query. Default is ``"Tritium"``.
    subsystem : str or None, optional
        Name of a subsystem for which to extract the steady-state value.
        If ``None``, the full array is returned.

    Returns
    -------
    numpy.ndarray or scalar
        If ``subsystem is None``, returns the full steady-state mass array.
        Otherwise returns the steady-state mass value for the selected subsystem.

    Raises
    ------
    ValueError
        If the series does not contain any snapshots, if the species is invalid,
        or if the subsystem is not found.
    """
    it = _get_first_snapshot(series)
    _validate_species(it, species)

    mass_steady = it.particles[species]["mass_steady"][io.Record_Component.SCALAR]
    mass_steady_data = mass_steady.load_chunk()
    series.flush()

    if subsystem is None:
        return mass_steady_data

    _validate_subsystem(it, species, subsystem)
    subsystem_id = it.particles[species]["subsystems"][subsystem].get_attribute("id")
    return mass_steady_data[subsystem_id]


def get_time_inventory(series, species="Tritium", subsystem=None):
    """Return the time-dependent mass inventory for a species.

    This function loads the ``mass`` data for the requested species from the
    first snapshot in the series. If a subsystem name is provided, only the
    time history for that subsystem is returned.

    Parameters
    ----------
    series : openPMD.Series
        OpenPMD RHINO series object containing mass inventory data.
    species : str, optional
        Name of the species to query. Default is ``"Tritium"``.
    subsystem : str or None, optional
        Name of a subsystem for which to extract the time-dependent mass.
        If ``None``, the full array is returned.

    Returns
    -------
    numpy.ndarray
        If ``subsystem is None``, returns the full mass array, typically indexed
        by subsystem and time. Otherwise returns the one-dimensional time series
        for the selected subsystem.

    Raises
    ------
    ValueError
        If the series does not contain any snapshots, if the species is invalid,
        or if the subsystem is not found.
    """
    it = _get_first_snapshot(series)
    _validate_species(it, species)

    mass = it.particles[species]["mass"][io.Record_Component.SCALAR]
    mass_data = mass.load_chunk()
    series.flush()

    if subsystem is None:
        return mass_data

    _validate_subsystem(it, species, subsystem)
    subsystem_id = it.particles[species]["subsystems"][subsystem].get_attribute("id")
    return mass_data[subsystem_id, :]

def get_subsystem_id(series, subsystem, species="Tritium"):
    """Return the integer ID associated with a subsystem name.

    Parameters
    ----------
    series : openPMD.Series
        OpenPMD RHINO series object containing subsystem metadata.
    subsystem : str
        Name of the subsystem to look up.
    species : str, optional
        Name of the species whose subsystem mapping should be used.
        Default is ``"Tritium"``.

    Returns
    -------
    int
        Integer ID associated with the requested subsystem.

    Raises
    ------
    ValueError
        If the series does not contain any snapshots, if the species is invalid,
        or if the subsystem is not available for that species.
    """
    it = _get_first_snapshot(series)
    _validate_species(it, species)
    _validate_subsystem(it, species, subsystem)
    return it.particles[species]["subsystems"][subsystem].get_attribute("id")


def get_subsystem_name(series, subsystem_id, species="Tritium"):
    """Return the subsystem name associated with a subsystem ID.

    Parameters
    ----------
    series : openPMD.Series
        OpenPMD RHINO series object containing subsystem metadata.
    subsystem_id : int
        Integer subsystem ID to look up.
    species : str, optional
        Name of the species whose subsystem mapping should be used.
        Default is ``"Tritium"``.

    Returns
    -------
    str
        Name of the subsystem associated with the given ID.

    Raises
    ------
    ValueError
        If the series does not contain any snapshots, if the species is invalid,
        if the subsystem ID is out of range, or if no subsystem matches the ID.
    """
    it = _get_first_snapshot(series)
    _validate_species(it, species)

    n_subsystems = len(it.particles[species]["subsystems"])
    if subsystem_id not in range(n_subsystems):
        raise ValueError(
            f"Invalid subsystem_id {subsystem_id!r}. "
            f"Valid IDs are: {list(range(n_subsystems))}"
        )

    for k, v in it.particles[species]["subsystems"].items():
        if v.get_attribute("id") == subsystem_id:
            return k

    raise ValueError(
        f"No subsystem found with id {subsystem_id!r} for species {species!r}."
    )

def get_subsystems_dict(series, species="Tritium"):
    """Return a dictionary mapping subsystem names to subsystem IDs.

    Parameters
    ----------
    series : openPMD.Series
        OpenPMD RHINO series object containing subsystem metadata.
    species : str, optional
        Name of the species whose subsystem mapping should be returned.
        Default is ``"Tritium"``.

    Returns
    -------
    dict[str, int]
        Dictionary whose keys are subsystem names and whose values are
        subsystem IDs.

    Raises
    ------
    ValueError
        If the series does not contain any snapshots or if the requested
        species is not available.
    """
    it = _get_first_snapshot(series)
    _validate_species(it, species)

    subsystems = {}
    for k, v in it.particles[species]["subsystems"].items():
        subsystems[k] = v.get_attribute("id")

    return subsystems

def get_series_input_attributes(series):
    inputs = {}
    for a in series.attributes: 
        if a.startswith("input:"):
            inputs[a] = series.get_attribute(a)
    return inputs 

def get_series_output_attributes(series):
    outputs = {}
    for a in series.attributes: 
        if a.startswith("output:"):
            outputs[a] = series.get_attribute(a)
    return outputs 

def get_steady_state_time(series):
    return series.get_attribute("output:Steady state time (days)")

def get_inputs_for_surrogate(series, input_specs):
    inputs = get_series_input_attributes(series)
    values = []

    for spec in input_specs:
        if spec["source"] != "series_input_attribute":
            raise ValueError(f"Unsupported input source: {spec['source']}")
        values.append(inputs[spec["attribute_key"]])
        x = np.array(values, dtype=np.float32)

    return x

def get_outputs_for_surrogate(series, output_specs):
    outputs = get_series_output_attributes(series)
    values = []
    for spec in output_specs:
        if spec["source"] == "series_output_attribute": 
            values.append(outputs[spec["attribute_key"]])

        elif spec["source"] == "particle_record_steady_state":
            values.append( get_steady_state(series, 
                                            species = spec["species"], 
                                            subsystem = spec["subsystem"]) 
                         )

        else:
            raise ValueError(f"Unsupported output source: {source}")         

    y = np.array(values, dtype=np.float32)
    if y.ndim == 0:
        y = y[None]

    return y


# In[ ]:


INPUT_SPECS = [
    {
        "key": "Ndotminus",
        "display_name": "Tritium burning rate [g/d]",
        "source": "series_input_attribute",
        "attribute_key": "input:Ndotminus:Tritium burned per day"
    },
    {
        "key": "beta",
        "display_name": "Burn fraction [-]",
        "source": "series_input_attribute",
        "attribute_key": "input:beta:Burn fraction"
    },
]

OUTPUT_SPECS = [
    {
        "key": "tritium_in_isotope_separation",
        "display_name": "Tritium in isotope separation [g]",
        "source": "particle_record_steady_state",
        "species": "Tritium",
        "subsystem": "Isotope_Seperation",
    },
    {
        "key": "plant_doubling_time_days",
        "display_name": "Plant doubling time [d]",
        "source": "series_output_attribute",
        "attribute_key": "output:plant_doubling_time (days)",
    },
    {
        "key": "minimum_startup_inventory_g",
        "display_name": "Minimum startup inventory [g]",
        "source": "series_output_attribute",
        "attribute_key": "output:I_startup (g)",
    },
]


# In[4]:


torch.manual_seed(42)
np.random.seed(42)


# In[5]:


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)


# In[6]:


def load_one_sample(series):

    x = get_inputs_for_surrogate(series, INPUT_SPECS)
    y  = get_outputs_for_surrogate(series, OUTPUT_SPECS)

    return x, y


# In[8]:


ROOT_PATH=Path("/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino/surrogate_bp_output")
SCENARIOS=["2026-04-29", "2026-04-30", "2026-05-01"]

X_list = []
Y_list = []

for scenario in SCENARIOS:
    DATA_PATH = ROOT_PATH / scenario 
    for bp5_file in sorted(DATA_PATH.glob("*.bp5")):

        series = io.Series(
            bp5_file,
            io.Access.read_only, 
            '{"verify_homogeneous_extents": false}'
        )
        x, y = load_one_sample(series)
        X_list.append(x)
        Y_list.append(y.reshape(-1))  
        series.close()


# In[9]:


X = np.stack(X_list, axis=0)       # shape: (N, d_in)
Y = np.stack(Y_list, axis=0)       # shape: (N, d_out)

_, INPUT_DIM = X.shape
N_SAMPLES, OUTPUT_DIM = Y.shape
print(f"Number of: inputs = {INPUT_DIM}, outputs {OUTPUT_DIM}, samples {N_SAMPLES}")
HIDDEN_DIM=100
N_HIDDEN=3


# In[10]:


print("X shape:", X.shape)
print("Y shape:", Y.shape)
print("First X row:", X[0])
print("First Y row shape:", Y[0].shape)


# In[11]:


print(X[:,0].min(), X[:,0].max())
print(X[:,1].min(), X[:,1].max())

print(Y[:,0].min(), Y[:,0].max())
print(Y[:,1].min(), Y[:,1].max())
print(Y[:,2].min(), Y[:,2].max())


# In[12]:


# Reproducible shuffle
rng = np.random.default_rng(42)

N = X.shape[0]
indices = np.arange(N)
rng.shuffle(indices)

# 70 / 15 / 15 split
n_train = int(0.80 * N)
n_val = int(0.15 * N)
n_test = N - n_train - n_val

train_idx = indices[:n_train]
val_idx = indices[n_train:n_train + n_val]
test_idx = indices[n_train + n_val:]

X_train = X[train_idx]
Y_train = Y[train_idx]

X_val = X[val_idx]
Y_val = Y[val_idx]

X_test = X[test_idx]
Y_test = Y[test_idx]

print("X_train shape:", X_train.shape)
print("Y_train shape:", Y_train.shape)
print("X_val shape:  ", X_val.shape)
print("Y_val shape:  ", Y_val.shape)
print("X_test shape: ", X_test.shape)
print("Y_test shape: ", Y_test.shape)


# In[13]:


# training-set stats only
x_mean = X_train.mean(axis=0, keepdims=True)
x_std  = X_train.std(axis=0, keepdims=True)

y_mean = Y_train.mean(axis=0, keepdims=True)
y_std  = Y_train.std(axis=0, keepdims=True)

# protect against divide-by-zero
x_std = np.where(x_std < 1e-12, 1.0, x_std)
y_std = np.where(y_std < 1e-12, 1.0, y_std)

# normalize
X_train_n = (X_train - x_mean) / x_std
X_val_n   = (X_val   - x_mean) / x_std
X_test_n  = (X_test  - x_mean) / x_std

Y_train_n = (Y_train - y_mean) / y_std
Y_val_n   = (Y_val   - y_mean) / y_std
Y_test_n  = (Y_test  - y_mean) / y_std

print("X_train_n mean:", X_train_n.mean(axis=0))
print("X_train_n std: ", X_train_n.std(axis=0))

print("Y_train_n mean:", Y_train_n.mean(axis=0))
print("Y_train_n std:", Y_train_n.std(axis=0))


# In[15]:


fig, ax = plt.subplots(1,1)
plt.scatter(X_train_n[:,0], Y_train_n[:,0])
plt.scatter(X_train_n[:,0], Y_train_n[:,1])
plt.scatter(X_train_n[:,0], Y_train_n[:,2])

#plt.scatter(X_train_n[:,1], Y_train_n[:,0])
#plt.scatter(X_train_n[:,1], Y_train_n[:,1])
#plt.scatter(X_train_n[:,1], Y_train_n[:,2])

plt.show()


# In[16]:


# Convert numpy arrays to torch tensors
X_train_t = torch.tensor(X_train_n, dtype=torch.float32)
Y_train_t = torch.tensor(Y_train_n, dtype=torch.float32)

X_val_t = torch.tensor(X_val_n, dtype=torch.float32)
Y_val_t = torch.tensor(Y_val_n, dtype=torch.float32)

X_test_t = torch.tensor(X_test_n, dtype=torch.float32)
Y_test_t = torch.tensor(Y_test_n, dtype=torch.float32)

print("X_train_t shape:", X_train_t.shape)
print("Y_train_t shape:", Y_train_t.shape)

# Wrap tensors into datasets
train_ds = TensorDataset(X_train_t, Y_train_t)
val_ds   = TensorDataset(X_val_t, Y_val_t)
test_ds  = TensorDataset(X_test_t, Y_test_t)

# Create dataloaders
batch_size = 64

train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
val_loader   = DataLoader(val_ds, batch_size=batch_size, shuffle=False)
test_loader  = DataLoader(test_ds, batch_size=batch_size, shuffle=False)

print("Number of training batches:", len(train_loader))
print("Number of validation batches:", len(val_loader))
print("Number of test batches:", len(test_loader))

# Peek at one batch
xb, yb = next(iter(train_loader))
print("One batch X shape:", xb.shape)
print("One batch Y shape:", yb.shape)


# In[17]:


class SurrogateMLP(nn.Module):
    def __init__(self, input_dim=INPUT_DIM, output_dim=OUTPUT_DIM, hidden_dim=100, num_hidden_layers=5):
        super().__init__()

        layers = []
        in_dim = input_dim

        for _ in range(num_hidden_layers):
            layers.append(nn.Linear(in_dim, hidden_dim))
            layers.append(nn.ReLU())
            in_dim = hidden_dim

        layers.append(nn.Linear(hidden_dim, output_dim))

        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)

model = SurrogateMLP(
    input_dim=INPUT_DIM,
    output_dim=OUTPUT_DIM,
    hidden_dim=100,
    num_hidden_layers=3
)

print(model)


# In[19]:


# Use GPU if available, otherwise CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

# Move model to device
model = model.to(device)

# Loss function for regression
criterion = nn.MSELoss()

# Optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)

print("Loss:", criterion)
print("Optimizer:", optimizer)


# In[20]:


def run_epoch(model, loader, criterion, optimizer=None, device="cpu"):
    if optimizer is None:
        model.eval()
    else:
        model.train()

    total_loss = 0.0
    total_samples = 0

    for xb, yb in loader:
        xb = xb.to(device)
        yb = yb.to(device)

        if optimizer is not None:
            optimizer.zero_grad()

        preds = model(xb)
        loss = criterion(preds, yb)

        if optimizer is not None:
            loss.backward()
            optimizer.step()

        batch_size = xb.size(0)
        total_loss += loss.item() * batch_size
        total_samples += batch_size

    return total_loss / total_samples


num_epochs = 200

train_losses = []
val_losses = []

best_val_loss = float("inf")
best_state = copy.deepcopy(model.state_dict())

for epoch in range(1, num_epochs + 1):
    train_loss = run_epoch(model, train_loader, criterion, optimizer=optimizer, device=device)
    val_loss   = run_epoch(model, val_loader, criterion, optimizer=None, device=device)

    train_losses.append(train_loss)
    val_losses.append(val_loss)

    if val_loss < best_val_loss:
        best_val_loss = val_loss
        best_state = copy.deepcopy(model.state_dict())

    if epoch % 20 == 0:
        print(f"Epoch {epoch:3d} | train_loss = {train_loss:.6f} | val_loss = {val_loss:.6f}")

# restore best model
model.load_state_dict(best_state)

print("\nBest validation loss:", best_val_loss)


# In[21]:


fig, ax = plt.subplots(1,1)
ax.plot(train_losses)
plt.show()


# In[22]:


model.eval()

with torch.no_grad():
    X_test_device = X_test_t.to(device)
    Y_test_device = Y_test_t.to(device)

    Y_pred_test_n = model(X_test_device).cpu().numpy()
    Y_test_n_np = Y_test_device.cpu().numpy()

# MSE and RMSE in normalized space
test_mse_n = np.mean((Y_pred_test_n - Y_test_n_np) ** 2)
test_rmse_n = np.sqrt(test_mse_n)
test_mae_n = np.mean(np.abs(Y_pred_test_n - Y_test_n_np))

print("Normalized test MSE: ", test_mse_n)
print("Normalized test RMSE:", test_rmse_n)
print("Normalized test MAE: ", test_mae_n)

# Convert predictions back to physical units
Y_pred_test = Y_pred_test_n * y_std + y_mean

# Metrics in original units
test_mse = np.mean((Y_pred_test - Y_test) ** 2)
test_rmse = np.sqrt(test_mse)
test_mae = np.mean(np.abs(Y_pred_test - Y_test))

print("\nPhysical-units test MSE: ", test_mse)
print("Physical-units test RMSE:", test_rmse)
print("Physical-units test MAE: ", test_mae)

print("\nPrediction array shape:", Y_pred_test.shape)


# In[23]:


plt.clf()
fig,ax = plt.subplots(ncols=3,nrows=1,figsize=(9,3), dpi=300)

for idx in range(OUTPUT_DIM):

    ax[idx].scatter(X_train[:,1], Y_train[:,idx], s=30, label= "train data")    
    ax[idx].scatter(X_test[:,1],  Y_test[:,idx], s=20, label= "test data")    
    ax[idx].scatter(X_test[:,1],  Y_pred_test[:,idx], s=10, label= "prediction")
    ax[idx].legend()
    ax[idx].set_xlabel(INPUT_SPECS[1]["display_name"])
    ax[idx].set_ylabel(OUTPUT_SPECS[idx]["display_name"])

plt.tight_layout()
plt.show()


# In[24]:


for idx in range(OUTPUT_DIM):
    plt.figure()

    y_true = Y_test[:, idx]
    y_pred = Y_pred_test[:, idx]

    plt.scatter(y_true, y_pred, s=10)

    # perfect prediction line
    min_val = min(y_true.min(), y_pred.min())
    max_val = max(y_true.max(), y_pred.max())
    plt.plot([min_val, max_val], [min_val, max_val])

    plt.xlabel("True")
    plt.ylabel("Predicted")
    plt.title(f"Output dimension {idx}")

    plt.show()


# In[25]:


# worst-case error across all outputs
max_abs_error = np.max(np.abs(Y_pred_test - Y_test))
print("Max absolute error:", max_abs_error)


# In[26]:


print("X shape:", X.shape)
print("Y shape:", Y.shape)
print("Y dtype:", Y.dtype)
print("Y nbytes (MB):", Y.nbytes / 1024**2)

print("Any NaNs in Y?", np.isnan(Y).any())
print("Any infs in Y?", np.isinf(Y).any())

for j in range(Y.shape[1]):
    col = Y[:, j]
    print(
        f"Output {j:2d} | min={col.min():.6e} max={col.max():.6e} "
        f"mean={col.mean():.6e} std={col.std():.6e}"
    )


# In[27]:


for j in [1, Y.shape[1]-2, Y.shape[1]-1]:
    col = Y[:, j]
    print(f"\nOutput {j}")
    print("first 20 values:", col[:20])
    print("unique count (rounded):", len(np.unique(np.round(col, 12))))
    print("std:", col.std())


# In[28]:


print("Y.shape:", Y.shape)
print("First row shape:", Y[0].shape)
print("Last row shape:", Y[-1].shape)
print("Y.shape:", Y.shape, "Y.nbytes MB:", Y.nbytes / 1024**2)


# In[29]:


y_std_raw = Y_train.std(axis=0)
print(y_std_raw)


# In[30]:


tol = 1e-6
varying_mask = y_std_raw > tol
constant_mask = ~varying_mask

print("Varying outputs:", np.where(varying_mask)[0])
print("Constant outputs:", np.where(constant_mask)[0])


# In[31]:


artifact = {
    "model_state_dict": model.state_dict(),
    "model_config": {
        "input_dim": INPUT_DIM,
        "output_dim": OUTPUT_DIM,
        "hidden_dim": HIDDEN_DIM,
        "num_hidden_layers": N_HIDDEN,
    },
    "input_specs": INPUT_SPECS,
    "output_specs": OUTPUT_SPECS,
    "x_mean": torch.tensor(x_mean, dtype=torch.float32),
    "x_std": torch.tensor(x_std, dtype=torch.float32),
    "y_mean": torch.tensor(y_mean, dtype=torch.float32),
    "y_std": torch.tensor(y_std, dtype=torch.float32),
}
torch.save(artifact, "rhino_surrogate.pt")

