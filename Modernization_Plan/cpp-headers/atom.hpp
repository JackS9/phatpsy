// include/phatpsy/core/atom.hpp
#pragma once
#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>

namespace phatpsy {

struct AtomicBasisFunction {
    int n, l, m;              // Quantum numbers
    double exponent;          // STO exponent
    double coefficient;       // Contraction coefficient
    Eigen::Vector3d center;   // Atomic center position
};

class Atom {
private:
    int atomic_number_;
    double nuclear_charge_;
    Eigen::Vector3d position_;
    std::vector<AtomicBasisFunction> basis_functions_;
    std::vector<double> occupation_alpha_, occupation_beta_;
    
public:
    Atom(int z, const Eigen::Vector3d& pos);
    
    // Accessors
    int atomic_number() const { return atomic_number_; }
    double nuclear_charge() const { return nuclear_charge_; }
    const Eigen::Vector3d& position() const { return position_; }
    
    // Basis function management
    void add_basis_function(int n, int l, int m, double exp, double coeff = 1.0);
    const std::vector<AtomicBasisFunction>& basis_functions() const { return basis_functions_; }
    size_t num_basis_functions() const { return basis_functions_.size(); }
    
    // Occupation management
    void set_occupation(size_t orbital_idx, double alpha_occ, double beta_occ = 0.0);
    double get_occupation_alpha(size_t orbital_idx) const;
    double get_occupation_beta(size_t orbital_idx) const;
    
    // PHATPSY-specific methods
    double total_electronic_charge() const;
    double net_atomic_charge() const;
    std::vector<double> compute_atomic_density_matrix() const;
};
