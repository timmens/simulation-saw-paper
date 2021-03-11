# TODO

- [ ] I noticed that the authors 


## Simulation Code

#### A Wavelet Method for Panel Models with Jump Discontinuities in the Parameters

### Introduction

Here we provide computer code which can be used to reproduce the simulation results in
the paper **A Wavelet Method for Panel Models with Jump Discontinuities in the Parameters**
_(Bada O., Kneip A., Liebl D., Gualtieri J. and Sickles R. C.)_.


The repository is structured as follows. In the folder ``src`` you will find code
that can be executed to produce simulation results and tables that are used in the
paper. The data which is used in the paper is also stored in the folder ``bld``. Instead
of the common approach to leave the ``bld`` folder empty such that everyone who wishes
to compare results has to run the codes themselves, we fill this folder since some codes
are computationally expensive to run and rely on proprietary software (MATLAB).


The implementation of the method presented in the paper is hosted here: [sawr](https://github.com/timmens/sawr)
and can be installed into R using

```R
library("devtools")
install_github("timmens/sawr")
```


### Simulation
