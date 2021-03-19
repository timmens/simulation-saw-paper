# TODO

- [ ] 5. R1: Setup DGP with time fixed effects.
- [ ] 5. R3: Simulation results n = 30. Might need to focus on SAW-estimator and compare quantity defined in Thm 1.
- [ ] 6. R3: DGP with endogeneous regressors, maybe delete DGP2 or DGP3 (Dominik favors deleting DGP3). Need to check if MATLAB codes can handle this case and if yes udpate code there to.


## Introduction

Here we provide computer code which can be used to reproduce the simulation results in
the paper **A Wavelet Method for Panel Models with Jump Discontinuities in the Parameters**
_(Bada O., Kneip A., Liebl D., Gualtieri J. and Sickles R. C.)_.


## How to Reproduce

Reproduction of all results is a bit tricky. This is due to the fact that the comparison
method is written in MATLAB and hence for a complete reproduction one needs to run the
MATLAB files stored in ``src/matlab``; more on that below. Given that the MATLAB
simulation results have been produced the rest is easy. You simply need to source the
file ``src/main.R``. It will run the simulations of our method and produce the overall
tex file which we used in the paper. The results can then be found in ``bld/tex`` and
``bld/R``.

**Note:**
Notice that the ``bld`` we provide in this repository is not empty, which is usually
the case in reproduction repositiories. We do this since the runtime of the MATLAB codes
is long and as MATLAB is proprietary not everyone has access to it. For a fresh run of
the simulation of our method remove the folder ``bld/R`` and ``bld/tex`` then source
``src/main.R``.

## Our Method 

The implementation of out method is hosted here: [sawr](https://github.com/timmens/sawr)
and can be installed into R using

```R
library("devtools")
install_github("timmens/sawr")
```

## Contact & Contribution

If you find bugs in our codes feel free to open an issue / pull-request or simply
contact me. See [here](https://github.com/timmens).
