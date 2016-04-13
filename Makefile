# 
# PHATPSY
# - Jack A Smith
#
FC 	= gfortran
FFLAGS 	= -O2 -pg -static
#
GNUDIR 	= /opt/lapack/gnu/lib
GNU_BLAS_LIBS 	= -lblas
GNU_LAPACK_LIBS = -llapack
GNU_SCALAPACK_LIBS = -lscalapack
GNU_FLIBS 	= -L${GNUDIR} ${GNU_BLAS_LIBS} ${GNU_LAPACK_LIBS} ${GNU_SCALAPACK_LIBS}
#
MKLDIR	= ${MKLROOT}/lib/intel64
# Intel ifort compiler
#MKL_BLAS_LIBS 	= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
# GNU gfortran compiler
# Sequential BLAS
#MKL_BLAS_LIBS 	= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
# Parallel BLAS
MKL_BLAS_LIBS 	= -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -lm -ldl
MKL_LAPACK_LIBS = -lmkl_blacs_openmpi_lp64
MKL_SCALAPACK_LIBS = -lmkl_scalapack_lp64 
MKL_FLIBS 	= -L${MKLDIR} ${MKL_BLAS_LIBS} ${MKL_LAPACK_LIBS} ${MKL_SCALAPACK_LIBS}
#
FLIBS	= ${MKL_FLIBS}
OBJS 	= absum.o addmat.o aloop.o analys.o anmbnm.o arrmap.o ascale.o asum.o atomic.o averag.o bloop.o \
	  bomb.o bscale.o cdct.o cgcoef.o chkcgc.o chrgbo.o contrl.o ctsc.o dampv.o dcoef.o dcopy.o denmat.o \
	  depth.o dmpab.o dmpatb.o dotprd.o eig.o eispak.o epkpot.o erase.o err243.o esca.o ewmo.o excite.o \
	  flip.o fokmat.o genbc.o gencgc.o gendc.o genfac.o gennrm.o genplm.o genrot.o gensab.o genvab.o \
	  hesbrg.o heseig.o hesvec.o imtqlv.o insert.o insrtd.o intnrm.o invert.o jacobi.o lowdin.o \
	  main.o mapndx.o matinv.o maxovl.o ndxcgc.o ndxd.o neword.o noble.o norhes.o normlz.o occmet.o \
	  ogen.o oneint.o order.o outpak.o output.o outvec.o pcct.o pkepot.o plotv.o popout.o punchv.o \
	  putcgc.o putdc.o putone.o putsym.o rbeta.o rotsab.o round.o scf.o schmid.o scmult.o sfct.o sum.o \
	  symop.o sympak.o sympr1.o sympr2.o sympro.o symxyz.o teindx.o timout.o \
	  tinvit.o trbak3.o tred3.o trisq.o tritri.o \
	  twoint.o twoout.o uthu.o wfct.o xyzmap.o zero.o 
#
phatpsy: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(FLIBS) -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $<
test:
	./phatpsy <n2.stdin >n2.stdout
clean:
	rm *.o phatpsy
