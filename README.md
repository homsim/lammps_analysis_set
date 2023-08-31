# About

This repository is a set of functions to analyse molecular dynamics simulation data. It works best with files computed by *LAMMPS* (especially its log- and trajectory-files). 
The repo is heavily dependent on the API of the Open Visualization Tool *Ovito* (*A. Stukowski, Modelling Simul. Mater. Sci. Eng. 18, 015012 (2010)*). The `ovito_env.yml` file provides the needed libraries to import in a conda virtual environment (to be replaced by a proper dependency handler):
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

The temperatures are computed via the respective mode velocities. 
### center-of-mass
```math
\boldsymbol{v}_{com}  = \frac{m_A \boldsymbol{v}_A + m_B \boldsymbol{v}_B}{m_A + m_B}
```
```math
E_{com} = \frac{1}{2} (m_A + m_B) v_{com}^2
```
```math
T_{com} = \frac{2 E_{com}}{3 k_B}
```
### vibrational
```math
\boldsymbol{v}_{vib,A}  = ((\boldsymbol{v}_A - \boldsymbol{v}_{com}) \cdot \hat{\boldsymbol{d}}_{AB}) \cdot \hat{\boldsymbol{d}}_{AB} 
```
with 
```math
\hat{\boldsymbol{d}}_{AB} = \frac{\boldsymbol{p}_B - \boldsymbol{p}_A}{ \left| \left | \boldsymbol{p}_B - \boldsymbol{p}_A \right| \right|} 
```
```math
E_{vib} = \frac{1}{2} (m_A v_{vib,A}^2 + m_B v_{vib,B}^2) 
```
```math
T_{vib} = \frac{2 E_{vib}}{ k_B}
```

### rotational
```math
\boldsymbol{v}_{rot,A} = \boldsymbol{v}_A - \boldsymbol{v}_{com} - \boldsymbol{v}_{str,A} 
```
```math
I = m_A R_A^2 + m_B R_B^2
```
```math
\boldsymbol{\omega}_A = \frac{\boldsymbol{R}_A \times \boldsymbol{v}_{rot,A}}{R_A} 
```
with
```math
\boldsymbol{R}_A = \boldsymbol{p}_A - \boldsymbol{COM}
```
```math
E_{rot} = \frac{1}{2} I \omega_A^2
```
```math
T_{rot} = \frac{E_{rot}}{k_B}
```
