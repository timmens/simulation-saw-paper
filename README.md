## Introduction

Here we provide computer code which can be used to reproduce the simulation results in
the paper

***A Wavelet Method for Panel Models with Jump Discontinuities in the Parameters***

by Bada O., Kneip A., Liebl D., Gualtieri J. and Sickles R. C.


## How to Reproduce

Reproduction of all results is a bit tricky. This is due to the fact that the comparison
method is written in MATLAB and hence for a complete reproduction one needs to run the
MATLAB files stored in ``src/matlab``; more on that below. Given that the MATLAB
simulation results have been produced the rest is easy: You simply need open the Rproject
``simulation.Rproj`` and then source the file ``src/main.R``.  It will run the
simulations of our method and produce the overall tex files which are used in
the paper. The results can then be found in ``bld/tex`` and ``bld/R``. Note
that you do not have to use the ``Rproj`` file as long as you adjust the
relevant paths.

**Run the MATLAB files:**

To run the matlab code you first need to download the software corresponding to the
method presented in [Qian and Su, 2016b](https://www.sciencedirect.com/science/article/abs/pii/S0304407615002377?via%3Dihub).
The software is available on the authors [website](http://jhqian.org/software/index.htm)
(download structb_panel.rar and unpack). The files need to be stored in ``src/matlab/matlab_src``.
Then in your MATLAB session you need to add the paths ``src/matlab/matlab_src`` and
``src/matlab/simulation_src``. At last simply run the files ``monte_carlo_study_dgp*``.
This runs the simulation and writes the results to ``bld/matlab``.

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
