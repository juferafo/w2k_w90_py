# WIEN2k and wannier90 python package

This repository contains an object-oriented python package which can be used to extract
and process programmatically the data obtained from a [WIEN2k](http://www.wien2k.at)
and/or [wannier90](http://www.wannier.org) calculations.
It also contains a density-matrix class and a quantum-mechanics subpackage which supports, 
among other features, occupation matrix manipulation such as basis transformation or 
angular/spin momentum projections.

## Requisites

### WIEN2K

This package is being implemented and tested for the versions WIEN2k_17.1 and WIEN2k_16.1

### Python

Although this package is tested for Python2.7, the [`__future__`](https://docs.python.org/3/howto/pyporting.html)
module was used to prevent compatibility problems with Python3.

### Python packages

* [NumPy](https://pypi.python.org/pypi/numpy)
* [Matplotlib](https://matplotlib.org)
* [SymPy](http://www.sympy.org/en/index.html)
