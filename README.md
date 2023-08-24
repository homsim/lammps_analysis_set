# About

This repository is a set of functions to analyse molecular dynamics simulation data. It works best with files computed by *LAMMPS* (especially its log- and trajectory-files). 
The repo is heavily dependent on the API of the Open Visualization Tool *Ovito* (*A. Stukowski, Modelling Simul. Mater. Sci. Eng. 18, 015012 (2010)*). The `ovito_env.yml` file provides the needed libraries to import in a conda virtual environment (will later be replaced by some sort of dependency handler):
```bash
conda env create -f ovito_env.yml
```
So far no unit conversion tool is implemented. The expected units in all provided files follow the [*real* unit style of *LAMMPS*](https://docs.lammps.org/units.html#description). Above all, this means that positions are in units of Angstrom and velocites in Angstrom/femtosecond.

# Methods and Usage

## Spectra

A set of methods to calculate the mode spectrum by Fourier transforming some given autocorrelation function (ACF).

```python
from spectra import Spectra

spec = Spectra(vacf_file = 'vacf.out', dipole_file = 'dipole.out', dt = 0.25)
```

The `*.out` files are expected to have the following format:

```
# Some header
[time in fs]\t[ACF in arb. units]
```


### Spectra.velocity()

Calculate the mode spectrum from the velocity autocorrelation function provided as `vacf_file`. 

```python
spec.velocity(plot = True)
```

### Spectra.dipole()

Calculate the mode spectrum from the autocorrelation function of the change in dipole moment provided as `dipole_file`. The method also applies the harmonic approximation ([https://doi.org/10.1063/1.1774986](https://doi.org/10.1063/1.1774986)) with a default temperature of 300 K (can be changed by `T` variable).

```python
spec.dipole(plot = True)
```

## Temperatures

A set of methods to calculate the mode temperatures of diatomic molecules. It calculates center-of-mass, vibrational and rotational temperatures by either creating bonds dependend on the vdW-radii of the molecules' atoms or by reading in a ReaxFF-bond files (can be computed by `fix reaxff/bonds` in *LAMMPS*).

```python
from temperatures import Temperatures

temps = Temperatures(traj_file = 'test.lammpstrj', bond_file = 'bonds.reaxff', 
                    dt = 0.25, types = [2, 2])
```

The temperatures are computed via the respective mode velocities. $$ x = {-b \pm \sqrt{b^2-4ac} \over 2a} $$
