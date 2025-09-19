// include/phatpsy/core/molecule.hpp
#pragma once
#include "atom.hpp"
#include <memory>

class Molecule {
private:
    std::vector<std::unique_ptr<Atom>> atoms_;
    std::string title_;
    int total_charge_;
    int multiplicity_;
    bool unrestricted_;
    
    // PHATPSY-specific data
    Eigen::MatrixXd overlap_matrix_;
    Eigen::MatrixXd pseudopotential_matrix_;
    std::vector<Eigen::MatrixXd> atomic_fock_matrices_;
    
public:
    Molecule(const std::string& title = "PHATPSY Calculation");
    
    // Molecular structure
    void add_atom(int atomic_number, const Eigen::Vector3d& position);
    void add_atom(std::unique_ptr<Atom> atom);
    
    size_t num_atoms() const { return atoms_.size(); }
    const Atom& atom(size_t i) const { return *atoms_[i]; }
    Atom& atom(size_t i) { return *atoms_[i]; }
    
    // Properties
    void set_charge_and_multiplicity(int charge, int multiplicity);
    void set_unrestricted(bool unrestricted) { unrestricted_ = unrestricted; }
    
    int total_charge() const { return total_charge_; }
    int multiplicity() const { return multiplicity_; }
    bool is_unrestricted() const { return unrestricted_; }
    
    // Molecular properties
    int total_electrons() const;
    size_t total_basis_functions() const;
    double nuclear_repulsion_energy() const;
    
    // PHATPSY-specific methods
    void build_overlap_matrix();
    void build_pseudopotential_matrix();
    void compute_atomic_fock_matrices();
    
    const Eigen::MatrixXd& overlap_matrix() const { return overlap_matrix_; }
    const Eigen::MatrixXd& pseudopotential_matrix() const { return pseudopotential_matrix_; }
};
