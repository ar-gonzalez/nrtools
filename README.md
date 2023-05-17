# Numerical Relativity Tools

*[AG, 17/05/2023]*

Python package to handle and produce initial data and evolutions with Elliptica and BAM for BHNS mergers (with the intention to generalize to more).

## Installation
`python setup.py install`

## Requirements
Works with Python 3 in addition to `numpy` and `matplotlib`. To produce 2D movies one needs `pyvista`, `itertools`, and `imageio`.
For waveform extraction, [watpy](https://git.tpi.uni-jena.de/core/watpy) is also necessary.

## Features
- `Initial_Data()` class that sets up a generic par file (`Parameter_File()` class) and batch script to produce initial data with Elliptica. 
	- Checks status of run
	- Produces metadata
	- Checks contraints and other variables
	- Output class to check out resulting binary's parameters
- `Evolution()` class to evolve the initial data and handle output plots and 2D movies
	- Waveform output to interface with WATPy and CoReDB waveforms
