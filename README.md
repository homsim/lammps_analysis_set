# About

This repository is a set of functions to analyse molecular dynamics simulation data. It works best with files computed by *LAMMPS* (especially its log- and trajectory-files). 
The repo is heavily dependent on the API of the Open Visualization Tool *Ovito* (*A. Stukowski, Modelling Simul. Mater. Sci. Eng. 18, 015012 (2010)*). The `ovito_env.yml` file provides the needed libraries to import in a conda virtual environment:
```bash
conda env create -f ovito_env.yml
```

# Methods and Usage

## Spectra

A set of methods to calculate the mode spectrum by fourier transforming some given autocorrelation function.

```python
from spectra import Spectra

spec = Spectra(vacf_file = 'vacf.out', dipole_file = 'dipole.out',
                plot = True, dt = 0.25)
```

The `*.out` files are expected to have the following format:

```
"time [fs]"\t"ACF"
```


### Spectra.velocity()

Calculate the mode spectrum from the velocity autocorrelation function provided as `vacf_file`. 

```python
spec.velocity()
```

### Spectra.velocity()

Calculate the mode spectrum from the autocorrelation function of the change in dipole moment provided as `dipole_file`. 

```python
spec.dipole()
```

## Temperatures

A set of methods to calculate the mode temperatures of diatomic molecules. It calculates center-of-mass, vibrational and rotational temperatures by either creating bonds dependend on the vdW-radii of the molecules' atoms or by reading in a ReaxFF-bond files (can be computed by `fix reaxff/bonds` in *LAMMPS*).

```python
from temperatures import Temperatures

temps = Temperatures(traj_file = 'test.lammpstrj', bond_file = 'bonds.reaxff', 
                    plot = True, dt = 0.25, types = [2, 2])
```

