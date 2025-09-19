<?php
// webmo/plugins/phatpsy/phatpsy.php
/**
 * PHATPSY Plugin for WebMO
 * Projected Hamiltonian Approach to Polyatomic Systems
 * 
 * This plugin integrates the modernized PHATPSY quantum chemistry
 * package with the WebMO web-based interface.
 * 
 * @author Jack A. Smith and contributors
 * @version 2.0.0
 */

class PhatpsyPlugin extends QMEngine {
    
    public $name = "PHATPSY";
    public $version = "2.0.0";
    public $description = "Projected Hamiltonian Approach to Polyatomic Systems";
    public $author = "Jack A. Smith";
    public $url = "https://github.com/JackS9/phatpsy";
    
    // Supported calculation types
    public $calculation_types = array(
        'sp' => 'Single Point Energy',
        'opt' => 'Geometry Optimization',
        'freq' => 'Vibrational Frequencies',
        'pop' => 'Population Analysis',
        'esca' => 'ESCA Core Binding Energies',
        'surf' => 'Surface Chemistry'
    );
    
    // Basis sets supported by PHATPSY
    public $basis_sets = array(
        'STO-DZ' => 'Slater Double-Zeta',
        'STO-DZP' => 'Slater Double-Zeta + Polarization',
        'STO-TZ' => 'Slater Triple-Zeta',
        'STO-TZP' => 'Slater Triple-Zeta + Polarization'
    );
    
    // PHATPSY-specific methods
    public $methods = array(
        'PHATPSY-HF' => 'Standard PHATPSY Hartree-Fock',
        'PHATPSY-GCHF' => 'Grand-Canonical Hartree-Fock',
        'PHATPSY-EWMO' => 'Energy-Weighted Maximum Overlap',
        'PHATPSY-LSP' => 'Localized Spin-Polarized (Future)'
    );
