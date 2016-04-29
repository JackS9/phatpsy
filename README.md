# phatpsy
PHATPSY - Projected Hamiltonian Approach to Polyatomic Systems

 * To get a copy:

  git clone https://github.com/JackS9/phatpsy.git
  
 * To compile (without Makefile):

  gfortran -c *.f

 * To build (without Makefile):

  gfortran *.o -o phatpsy

 * To test (without Makefile):

  ./phatpsy \<n2.stdin \>n2.stdout

 * To use Makefile and access MKL Libraries

  module load mkl  # once only after logging in
  
  make clean  # optional
  
  make

  make test

