# tests/test_phatpsy_validation.py
"""
Comprehensive validation suite for modernized PHATPSY
Ensures accuracy against original Fortran implementation and literature values
"""

import numpy as np
import pytest
import tempfile
import os
from pathlib import Path
import subprocess

import phatpsy
from phatpsy import Atom, Molecule, PhatpsyCalculation
from phatpsy.analysis import PopulationAnalysis
from phatpsy.io import read_xyz

class TestOriginalValidation:
    """Validate against original PHATPSY Fortran results"""
    
    def setup_method(self):
        """Setup test environment"""
        self.test_data_dir = Path(__file__).parent / "data"
        self.reference_results = self.load_reference_data()
        
    def load_reference_data(self):
        """Load reference results from original PHATPSY calculations"""
        # Data from Jack Smith's 1978 dissertation
        return {
            'NO': {
                'total_energy': -129.2437,  # hartree
                'ionization_potentials': [9.9, 17.2, 17.5, 18.1, 18.4, 21.9, 21.4],
                'bond_length': 1.10,  # angstrom
                'vibrational_frequency': 1993,  # cm^-1
                'atomic_charges': {'N': -0.23, 'O': 0.23}
            },
            'H2': {
                'total_energy': -1.1315,
                'bond_length': 0.74,
                'dissociation_energy': 4.7
            },
            'CO': {
                'total_energy': -112.7854,
                'bond_length': 1.13,
                'atomic_charges': {'C': 0.15, 'O': -0.15}
            }
        }
