# About

This repository is heavily dependent on the API of the Open Visualization Tool *ovito* (*A. Stukowski, Modelling Simul. Mater. Sci. Eng. 18, 015012 (2010)*). The `ovito_env.yml` file provides the needed libraries to import in a conda virtual environment:
`conda env create -f ovito_env.yml`

# Methods and Usage

## Spectra

A set of methods to calculate the mode spectrum given some autocorrelation function.

```
from spectra import Spectra

spec = Spectra(vacf_file = 'vacf.out', dipole_file = 'dipole_test')
```

### Spectra.velocity()

Calculate 
