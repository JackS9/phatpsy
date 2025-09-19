# python/phatpsy/calculation.py
class PhatpsyCalculation:
    """Main calculation class for PHATPSY computations"""
    
    def __init__(self, molecule: Molecule, parameters: Optional[Dict] = None):
        self.molecule = molecule
        self.parameters = PhatpsyParameters()
        
        if parameters:
            for key, value in parameters.items():
                self.parameters.set(key, value)
        
        self._cpp_calc = _phatpsy_core.PhatpsyCalculation(
            molecule._cpp_molecule, self.parameters._cpp_params
        )
    
    def run_scf(self) -> 'SCFResult':
        """Run self-consistent field calculation"""
        print("Starting PHATPSY SCF calculation...")
        
        # Initialize atomic calculations
        self._initialize_atomic_calculations()
        
        # Main SCF loop with PHATPSY methodology
        result = self._cpp_calc.run_scf()
        
        return SCFResult(result)
    
    def run_ewmo(self) -> np.ndarray:
        """Run Energy-Weighted Maximum Overlap calculation"""
        if not self.parameters.use_ewmo():
            raise ValueError("EWMO calculation not enabled")
        
        mo_coefficients = self._cpp_calc.compute_ewmo_orbitals()
        return mo_coefficients
    
    def _initialize_atomic_calculations(self):
        """Initialize separated-atom calculations (PHATPSY methodology)"""
        print("Initializing atomic calculations...")
        for i, atom in enumerate(self.molecule.atoms):
            print(f"  Atom {i+1}: {self._number_to_symbol(atom.atomic_number)}")
            # Setup initial atomic configuration
            self._setup_atomic_configuration(atom)
