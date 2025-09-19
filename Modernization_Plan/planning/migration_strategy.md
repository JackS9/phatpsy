# PHATPSY Component Migration Strategy

## Migration Decision Matrix

| Component | Fortran File(s) | Migration Strategy | Priority | Complexity | Risk |
|-----------|----------------|-------------------|----------|------------|------|
| **Mathematical Utilities** | `erase.f`, `zero.f`, `sum.f`, `asum.f` | **Rewrite** - C++ Templates | HIGH | LOW | LOW |
| **Angular Momentum** | `cgcoef.f`, `gencgc.f`, `ndxcgc.f` | **Rewrite** - Optimized C++ | HIGH | MEDIUM | MEDIUM |
| **Integral Evaluation** | `oneint.f`, `twoint.f` | **Hybrid** - Core algorithms + Modern optimization | HIGH | HIGH | HIGH |
| **Pseudopotentials** | `epkpot.f`, `pkepot.f` | **Direct Port** - Preserve exact algorithms | HIGH | MEDIUM | MEDIUM |
| **SCF Engine** | `scf.f`, `atomic.f` | **Rewrite** - Modern SCF with PHATPSY features | HIGH | HIGH | MEDIUM |
| **EWMO Calculator** | `ewmo.f` | **Direct Port** - Unique PHATPSY algorithm | CRITICAL | HIGH | HIGH |
| **Analysis Tools** | `analys.f`, `popout.f` | **Rewrite** - Enhanced with Python integration | MEDIUM | MEDIUM | LOW |
| **I/O Systems** | `output.f`, `contrl.f` | **Rewrite** - Modern file formats | MEDIUM | LOW | LOW |
| **Matrix Operations** | `invert.f`, `jacobi.f`, `eig.f` | **Replace** - Use Eigen3/LAPACK | MEDIUM | LOW | LOW |

## Detailed Migration Strategies

### 1. Mathematical Utilities - Template-Based Rewrite
**Rationale:** Simple functions that benefit greatly from C++ templates and type safety

```cpp
// Before (Fortran): erase.f
//       SUBROUTINE ERASE(A, N)
//       DO 10 I = 1, N
//         A(I) = 0.0D0
//  10   CONTINUE

// After (C++): math/utilities.hpp
template<typename T>
void zero_array(T* array, size_t n) {
    std::fill(array, array + n, T{0});
}

template<typename Container>
void zero_container(Container& container) {
    std::fill(container.begin(), container.end(), typename Container::value_type{0});
}
```

### 2. EWMO Calculator - Exact Algorithmic Port
**Rationale:** This is the unique PHATPSY contribution - must preserve exactly

```cpp
class EWMOCalculator {
public:
    struct EWMOResult {
        Eigen::MatrixXd molecular_orbitals;
        Eigen::VectorXd orbital_energies;
        Eigen::VectorXd occupation_numbers;
        bool converged;
    };
    
    EWMOResult compute_ewmo_orbitals(
        const std::vector<Atom>& atoms,
        const std::vector<Eigen::MatrixXd>& atomic_fock_matrices,
        const Eigen::MatrixXd& overlap_matrix
    ) const;
};
```

This comprehensive plan provides a clear roadmap for modernizing PHATPSY while preserving its unique theoretical advantages and integrating it seamlessly with WebMO's educational and research platform.
