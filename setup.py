#!/usr/bin/env python
from setuptools import setup, find_packages
import re
import sys

__ascii_art__ = """
.__   __. .______     .___________. ______     ______    __          _______.
|  \ |  | |   _  \    |           |/  __  \   /  __  \  |  |        /       |
|   \|  | |  |_)  |   `---|  |----|  |  |  | |  |  |  | |  |       |   (----`
|  . `  | |      /        |  |    |  |  |  | |  |  |  | |  |        \   \    
|  |\   | |  |\  \----.   |  |    |  `--'  | |  `--'  | |  `----.----)   |   
|__| \__| | _| `._____|   |__|     \______/   \______/  |_______|_______/  
"""

# Read version from pyproject.toml
ini = open('pyproject.toml').read()
vrs = r"^version = ['\"]([^'\"]*)['\"]"
mo  = re.search(vrs, ini, re.M)
version = mo.group(1)

setup(
    name='nrtools',
    version=version,
    description='Numerical Relativity Tools',
    author='AG',
    author_email='alejandra.gonzalez@uni-jena.de',
    url = 'https://git.tpi.uni-jena.de/agonzalez/nrtools',
    packages = find_packages(),
    requires = ['numpy', 'matplotlib'],
)

if 'install' in sys.argv:
    print(__ascii_art__)
    print("Check the tutorials folder to get started.")
    print("~May the force be with you~")
    