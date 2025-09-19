# PHATPSY Routine Dictionary

## Overview
This document provides a comprehensive dictionary of all routines in the PHATPSY quantum chemistry package, organized by functional categories. PHATPSY is a Fortran-based package for molecular orbital calculations using a Projected Hamiltonian Approach.

## Main Program Structure

### **Main Program**
- **`PHATPSY`** (`main.f`) - Main program entry point that reads control data and initiates calculations

### **Primary Control Routines**
- **`CONTRL`** (`contrl.f`) - Main driving routine that orchestrates the entire calculation
- **`ALOOP`** (`aloop.f`) - Main loop for atomic orbital calculations
- **`BLOOP`** (`bloop.f`) - Secondary loop for basis function generation
- **`ATOMIC`** (`atomic.f`) - Handles atomic calculations and integral generation
- **`SCF`** (`scf.f`) - Self-Consistent Field calculation driver
- **`EWMO`** (`ewmo.f`) - Energy-Weighted Maximum Overlap calculations

## Core Calculation Routines

### **Integral Calculations**
- **`ONEINT`** (`oneint.f`) - One-electron integral calculations (overlap, kinetic, potential)
- **`TWOINT`** (`twoint.f`) - Two-electron integral calculations (Coulomb, exchange)
- **`GENSAB`** (`gensab.f`) - S-orbital basis generation for integrals
- **`GENVAB`** (`genvab.f`) - V-orbital basis generation for integrals
- **`WFCT`** (`wfct.f`) - Wave function calculations for integrals

### **Matrix Operations**
- **`JACOBI`** (`jacobi.f`) - Jacobi diagonalization for symmetric matrices
- **`LOWDIN`** (`lowdin.f`) - Lowdin orthogonalization and matrix transformations
- **`EIG`** (`eig.f`) - General eigenvalue problem solver
- **`EISPAK`** (`eispak.f`) - EISPACK eigenvalue/eigenvector routines
- **`TRED3`** (`tred3.f`) - Tridiagonal matrix reduction
- **`IMTQLV`** (`imtqlv.f`) - Implicit QL algorithm for eigenvalues
- **`TINVIT`** (`tinvit.f`) - Inverse iteration for eigenvectors
- **`TRBAK3`** (`trbak3.f`) - Back transformation of eigenvectors

### **SCF Components**
- **`FOKMAT`** (`fokmat.f`) - Fock matrix construction
- **`DENMAT`** (`denmat.f`) - Density matrix construction
- **`ANALYS`** (`analys.f`) - Analysis of molecular properties
- **`MAXOVL`** (`maxovl.f`) - Maximum overlap calculations
- **`NEWORD`** (`neword.f`) - New orbital ordering calculations
- **`UTHU`** (`uthu.f`) - UTHU transformation calculations
- **`TRISQ`** (`trisq.f`) - Triangular matrix operations

### **Special Functions and Scaling**
- **`ASCALE`** (`ascale.f`) - Scaled A-function generation
- **`BSCALE`** (`bscale.f`) - Scaled B-function generation
- **`ANMBNM`** (`anmbnm.f`) - Modified A and B function generation
- **`SFCT`** (`sfct.f`) - S-function calculations
- **`GENPLM`** (`genplm.f`) - Associated Legendre polynomial generation
- **`GENNRM`** (`gennrm.f`) - Normalization factor generation

## Symmetry and Transformation Routines

### **Symmetry Operations**
- **`SYMOP`** (`symop.f`) - General symmetry operations
- **`SYMXYZ`** (`symxyz.f`) - XYZ coordinate symmetry operations
- **`SYMPRO`** (`sympro.f`) - Symmetry projection operations
- **`SYMPR1`** (`sympr1.f`) - Primary symmetry operations
- **`SYMPR2`** (`sympr2.f`) - Secondary symmetry operations
- **`SYMPAK`** (`sympak.f`) - Symmetry package utilities

### **Coordinate and Rotation**
- **`XYZMAP`** (`xyzmap.f`) - XYZ coordinate mapping
- **`GENROT`** (`genrot.f`) - Rotation matrix generation
- **`ROTSAB`** (`rotsab.f`) - S-orbital rotation operations
- **`DEPTH`** (`depth.f`) - Depth calculations for coordinates

### **Clebsch-Gordan Coefficients**
- **`GENCGC`** (`gencgc.f`) - Clebsch-Gordan coefficient generation
- **`PUTCGC`** (`putcgc.f`) - Clebsch-Gordan coefficient output
- **`CGCOEF`** (`cgcoef.f`) - Clebsch-Gordan coefficient retrieval
- **`CHKCGC`** (`chkcgc.f`) - Clebsch-Gordan coefficient validation
- **`NDXCGC`** (`ndxcgc.f`) - Clebsch-Gordan coefficient indexing

## Basis and Coefficient Generation

### **Basis Functions**
- **`GENBC`** (`genbc.f`) - Binomial coefficient generation
- **`GENFAC`** (`genfac.f`) - Factorial generation
- **`PUTDC`** (`putdc.f`) - D-coefficient output
- **`DCOEF`** (`dcoef.f`) - D-coefficient calculations

### **Indexing and Mapping**
- **`NDXD`** (`ndxd.f`) - D-function indexing
- **`MAPNDX`** (`mapndx.f`) - General index mapping
- **`ARRMAP`** (`arrmap.f`) - Generates index of arrays dynamically allocated ASSARD (ASSign ARrays Dynamically)

## Utility and Support Routines

### **Matrix Utilities**
- **`ADDMAT`** (`addmat.f`) - Matrix addition
- **`SUBMAT`** (`submat.f`) - Matrix subtraction
- **`INSERT`** (`insert.f`) - Matrix element insertion
- **`INSRTD`** (`insrtd.f`) - Matrix insertion with dimensions
- **`INVERT`** (`invert.f`) - Matrix inversion
- **`MATINV`** (`matinv.f`) - Alternative matrix inversion
- **`SCMULT`** (`scmult.f`) - Scalar matrix multiplication
- **`ZERO`** (`zero.f`) - Zero array checking

### **Array Operations**
- **`DCOPY`** (`dcopy.f`) - Array copying
- **`ERASE`** (`erase.f`) - Array erasure
- **`QERASE`** - Logical array erasure
- **`ASUM`** (`asum.f`) - Array sum calculations
- **`ABSUM`** (`absum.f`) - A and B function sums
- **`DOTPRD`** (`dotprd.f`) - Dot product calculations
- **`SUM`** (`sum.f`) - Sum of evenly spaced elements

### **Input/Output and Control**
- **`OUTPUT`** (`output.f`) - Matrix output formatting
- **`OUTPAK`** (`outpak.f`) - Matrix output packaging
- **`OUTVEC`** (`outvec.f`) - Vector output
- **`PUTMAT`** (`putmat.f`) - Matrix output
- **`PUTONE`** (`putone.f`) - One-electron matrix output
- **`PUTSYM`** (`putsym.f`) - Symmetry matrix output
- **`TWOOUT`** (`twoout.f`) - Two-electron integral output
- **`TIMOUT`** (`timout.f`) - Time output and management
- **`GETSET`** (`timout.f`) - Time get/set operations

## Analysis and Output Routines

### **Molecular Analysis**
- **`CHGAN`** (`analys.f`) - Charge analysis
- **`POPOUT`** (`popout.f`) - Population output
- **`PUNCHV`** (`punchv.f`) - Vector punching for restart
- **`PLOTV`** (`plotv.f`) - Vector plotting
- **`EXCITE`** (`excite.f`) - Excitation calculations
- **`OCCMET`** (`occmet.f`) - Occupation matrix operations

### **Energy and Properties**
- **`AVERAG`** (`averag.f`) - Expectation value calculations
- **`PKEPOT`** (`pkepot.f`) - Phillips-Kleinman potential calculations
- **`EPKPOT`** (`epkpot.f`) - Extended Phillips-Kleinman potential calculations
- **`RBETA`** (`rbeta.f`) - Beta parameter calculations
- **`RBETA2`** (`rbeta.f`) - Alternative beta calculations

## Advanced Mathematical Routines

### **Eigenvalue Methods**
- **`HESEIG`** (`heseig.f`) - Hessenberg eigenvalue calculations
- **`HESVEC`** (`hesvec.f`) - Hessenberg eigenvector calculations
- **`HESBRG`** (`hesbrg.f`) - Hessenberg matrix reduction
- **`JACOBI`** (`jacobi.f`) - Jacobi diagonalization

### **Specialized Calculations**
- **`OGEN`** (`ogen.f`) - Omega coefficient generation
- **`CTSC`** (`ctsc.f`) - Coefficient array to Overlap matrix transformation
- **`CSCT`** (`ctsc.f`) - Overlap matrix to Coefficient array transformation  
- **`PCCT`** (`pcct.f`) - Projection matrix to Coefficient array transformation
- **`TRITRI`** (`tritri.f`) - Triangular matrix operations
- **`NORMLZ`** (`normlz.f`) - Vector normalization
- **`SCHMID`** (`schmid.f`) - Schmidt orthogonalization

## Error Handling and Debugging

### **Error Management**
- **`BOMB`** (`bomb.f`) - Error termination routine
- **`ERR243`** (`err243.f`) - Error code 243 handler
- **`ERRTRA`** - Error trace routine
- **`SETBOM`** - Set bomb routine

### **Debugging and Tracing**
- **`DUMP`** - Debug dump routine
- **`QTRACE`** - Trace flag for debugging
- **`QDEBUG`** - Debug mode flag

## Time and Performance

### **Timing Routines**
- **`ETIME`** - Execution time measurement
- **`SECNDS`** (`timout.f`) - Seconds since reference time
- **`TIMEND`** - Time limit checking

## Specialized Functions

### **Physical Constants**
- **`NOBLE`** (`noble.f`) - Noble gas orbital counting
- **`ROUND`** (`round.f`) - Number rounding function

### **Array Management**
- **`DMPAB`** (`dmpab.f`) - Double precision Matrix Product of A and B functions
- **`DMPATB`** (`dmpatb.f`) - Double precision Matrix Product of A, T, and B functions
- **`INTNRM`** (`intnrm.f`) - Integral normalization

## Common Naming Conventions

Throughout the PHATPSY codebase, several consistent naming conventions are used:

### **Matrix and Array Naming:**
- **S** - Overlap matrix (S-matrix)
- **C** - Coefficient array (molecular orbital coefficients)
- **P** - Projection or density matrix
- **T** - Transpose operations (often in routine names)
- **H** - Hamiltonian matrix
- **F** - Fock matrix
- **V** - Potential matrix (nuclear attraction)
- **K** - Kinetic energy matrix

### **Function and Variable Naming:**
- **A, B** - Special functions (A-functions and B-functions)
- **W** - Wave functions or W-functions
- **V** - Potential functions or V-functions
- **D** - D-coefficients or D-functions
- **L, M** - Angular momentum quantum numbers
- **N** - Principal quantum number
- **ETA** - Exponent parameters

### **Routine Naming Patterns:**
- **GEN*** - Generate routines (e.g., GENCGC, GENSAB, GENVAB)
- **PUT*** - Output routines (e.g., PUTCGC, PUTDC, PUTMAT)
- **NDX*** - Index routines (e.g., NDXCGC, NDXD)
- **SYM*** - Symmetry routines (e.g., SYMOP, SYMXYZ, SYMPRO)

## Routine Categories Summary

### **Core Calculation Routines (High Priority)**
- Integral calculations (ONEINT, TWOINT, GENSAB, GENVAB)
- Matrix operations (JACOBI, LOWDIN, EIG, EISPAK)
- SCF components (FOKMAT, DENMAT, ANALYS)
- Special functions (ASCALE, BSCALE, ANMBNM)

### **Symmetry and Transformation (Medium Priority)**
- Symmetry operations (SYMOP, SYMXYZ, SYMPRO)
- Coordinate transformations (GENROT, ROTSAB)
- Clebsch-Gordan coefficients (GENCGC, PUTCGC, CGCOEF)

### **Utility and Support (Lower Priority)**
- Matrix utilities (ADDMAT, SUBMAT, INSERT, INVERT)
- Array operations (DCOPY, ERASE, ASUM, DOTPRD)
- Input/Output (OUTPUT, OUTPAK, PUTMAT, TIMOUT)

### **Analysis and Output (Medium Priority)**
- Molecular analysis (CHGAN, POPOUT, PLOTV)
- Energy calculations (AVERAG, PKEPOT, RBETA)
- Advanced math (HESEIG, HESVEC, SCHMID)

### **Error Handling and Debugging (Low Priority)**
- Error management (BOMB, ERR243, ERRTRA)
- Debugging (DUMP, QTRACE, QDEBUG)

## Integration Points

The routines are integrated through several key calling patterns:
- **CONTRL** calls major calculation drivers (ATOMIC, ALOOP, SCF, EWMO)
- **SCF** orchestrates self-consistent field calculations
- **ATOMIC** handles atomic-level calculations and integral generation
- **BLOOP** manages basis function generation and integration
- **LOWDIN** provides matrix transformation services
- **EISPAK** provides eigenvalue/eigenvector services

This dictionary provides a comprehensive overview of all PHATPSY routines, their purposes, and their relationships within the quantum chemistry calculation framework. 