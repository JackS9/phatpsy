# phatpsy
PHATPSY - Projected Hamiltonian Approach to Polyatomic Systems

 * To get a copy:

  git clone https://github.com/JackS9/phatpsy.git
  
 * To compile:

  gfortran -c *.f

 * To build:

  gfortran *.o -o phatpsy

 * To test:

  ./phatpsy \<n2.stdin \>n2.stdout

 * To use MKL Libraries and Makefile first time

  module load mkl  # once only

  make

  make test

 * To use MKL Libraries and Makefile subsequent times 

  make clean  # optional

  make

  make test
