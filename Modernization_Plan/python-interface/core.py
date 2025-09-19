# python/phatpsy/core.py
import numpy as np
from typing import List, Tuple, Optional, Dict, Any
from . import _phatpsy_core  # C++ extension module

class Atom:
    """Python wrapper for PHATPSY Atom class"""
    
    def __init__(self, atomic_number: int, position: np.ndarray):
        self._cpp_atom = _phatpsy_core.Atom(atomic_number, position)
        
    @property
    def atomic_number(self) -> int:
        return self._cpp_atom.atomic_number()
    
    @property
    def position(self) -> np.ndarray:
        return self._cpp_atom.position()
    
    def add_basis_function(self, n: int, l: int, m: int, 
                          exponent: float, coefficient: float = 1.0):
        """Add a Slater-type orbital basis function"""
        self._cpp_atom.add_basis_function(n, l, m, exponent, coefficient)
    
    def set_occupation(self, orbital_idx: int, alpha_occ: float, 
                      beta_occ: float = 0.0):
        """Set orbital occupation numbers (supports fractional occupations)"""
        self._cpp_atom.set_occupation(orbital_idx, alpha_occ, beta_occ)
    
    def get_atomic_charge(self) -> float:
        """Get net atomic charge using PHATPSY population analysis"""
        return self._cpp_atom.net_atomic_charge()

class Molecule:
    """Python wrapper for PHATPSY Molecule class"""
    
    def __init__(self, title: str = "PHATPSY Calculation"):
        self._cpp_molecule = _phatpsy_core.Molecule(title)
        self.atoms = []
    
    def add_atom(self, atomic_number: int, position: np.ndarray) -> Atom:
        """Add an atom to the molecule"""
        atom = Atom(atomic_number, position)
        self.atoms.append(atom)
        self._cpp_molecule.add_atom(atom._cpp_atom)
        return atom
    
    def from_xyz(self, xyz_string: str):
        """Create molecule from XYZ format string"""
        lines = xyz_string.strip().split('\n')
        natoms = int(lines[0])
        self.title = lines[1] if len(lines) > 1 else "PHATPSY Calculation"
        
        for i in range(2, 2 + natoms):
            parts = lines[i].split()
            symbol = parts[0]
            coords = np.array([float(x) for x in parts[1:4]])
            
            # Convert atomic symbol to number
            atomic_number = self._symbol_to_number(symbol)
            self.add_atom(atomic_number, coords)
