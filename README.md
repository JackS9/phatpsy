# phatpsy
PHATPSY - Projected Hamiltonian Approach to Polyatomic Systems

 * To get a copy:

  git clone https://github.com/JackS9/phatpsy.git
  
 * To compile (without Makefile):

  gfortran -c src/*.f

 * To build (without Makefile):

  gfortran src/*.o -o bin/phatpsy

 * To test (without Makefile):

  ./bin/phatpsy \<examples/n2.stdin \>data/n2.stdout

 * To use Makefile and access MKL Libraries

   + module load mkl  # once only after logging in
   + cd src
   + make clean  # optional
   + make
   + make test

