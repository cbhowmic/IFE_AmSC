# LASER HDF5 to BP5 Shim

The LASER shim converts each source HDF5 acquisition into one openPMD/ADIOS2
BP5 dataset. A source acquisition can contain many valid laser shots.

## NERSC batch conversion

Submit the supplied Slurm job from the repository root:

```bash
sbatch LASER/convert_laser.slurm
```

The job defaults to the `m3239` LASER data location and the `IFE_AmSC` Conda
environment. Override locations or the environment at submission time when
needed:

```bash
sbatch \
  --export=ALL,LASER_DATA_ROOT=/path/to/laser,LASER_BP_OUT_DIR=/path/to/bp_output,LASER_CONDA_ENV=IFE_AmSC \
  LASER/convert_laser.slurm
```

Each HDF5 file is converted in a separate Python process so memory is released
between acquisitions. Existing BP5 outputs are skipped, making the job safe to
restart after a timeout or failure. Failures are reported and cause the final
job status to be nonzero.

## Interactive conversion

For small tests on an appropriate node:

```bash
python LASER/laserWrite_multi.py /path/to/laser \
  --recursive \
  --out-dir /path/to/bp_output \
  --skip-existing
```

Use `--stop-on-error` to stop immediately when an input cannot be converted.

## Count outputs and shots

```bash
python LASER/count_outputs.py /path/to/bp_output
```

This reports the number of BP5 acquisitions and sums the
`input:num_valid_samples` attribute to obtain the number of valid shots.
