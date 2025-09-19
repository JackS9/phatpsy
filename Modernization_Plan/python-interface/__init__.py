# python/phatpsy/__init__.py
"""
PHATPSY: Projected Hamiltonian Approach to Polyatomic Systems
Modern Python interface to the PHATPSY quantum chemistry package
"""

from .core import Atom, Molecule, PhatpsyParameters
from .calculation import PhatpsyCalculation, SCFResult
from .analysis import PopulationAnalysis, PropertyCalculator
from .io import read_xyz, write_output, load_calculation
from .visualization import plot_molecular_orbitals, plot_density

__version__ = "2.0.0"
__author__ = "Jack A. Smith and contributors"
