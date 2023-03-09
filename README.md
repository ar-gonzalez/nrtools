# Numerical Relativity Tools

~*Under Construction*~

Python package to handle and produce initial data and evolutions with Elliptica and BAM for BHNS mergers (with the intention to generalize to more).

## Installation
`python setup.py install`

## Requirements
Works with Python 3 in addition to `numpy` and `matplotlib`.

## Features
- `Initial_Data()` class that sets up a generic par file (`Parameter_File()` class) and batch script to produce initial data with Elliptica. 
	- Checks status of run
	- Produces metadata
- Coming soon:
	- Diagnostics class to check contraints and other variables
	- Evolution class to evolve the initial data and handle outputs (plots and movies)
	- Waveform class to interface with WATPy and CoReDB waveforms