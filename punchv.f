      SUBROUTINE PUNCHV(STEVEC,EVAL,OCC,L,NSTO,MSTO,NORB,ISPIN,IP)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C
C     PUNCHV...
C
C        THIS ROUTINE PUNCHES OUT THE FINAL ATOMIC EIGENVECTORS IN A
C     FORM READABLE BY ATOMIC (INPUT ROUTINE).  THE OUTPUT INCLUDES THE
C     NAMELIST THAT PRECEDES THE VECTORS.
C
C     VARIABLE DEFINITIONS:
C
C        STEVEC(*,*,*)... ATOMIC EIGENVECTORS OVER STO BASIS.
C        EVAL(*,*)....... CORRESPONDING ATOMIC ORBITAL EIGENVALUES.
C        OCC(*,*)........ ATOMIC ORBITAL OCCUPATION NUMBERS.
C        L(*)............ L-QUANTUM NUMBERS OF STO'S.
C        NSTO............ NUMBER OF STO'S IN BASIS (NOT INCLUDING ML'S).
C        MSTO............ NUMBER OF STO'S IN BASIS (INCLUDING ML'S).
C        NORB............ NUMBER OF ATOMIC ORBITALS (EIGENVECTORS).
C        ISPIN........... =1 --> SPIN-RESTRICTED, =2 --> UNRESTRICTED.
C        IP.............. FORTRAN I/O UNIT FOR PUNCHING.
C
C-----------------------------------------------------------------------
      DIMENSION STEVEC(MSTO,NORB,ISPIN), EVAL(NORB,ISPIN),
     X          OCC(NORB,ISPIN),         L(NSTO),
     X          EIGVAL(2),               OCCNO(2),
     X          COEF(2)
      NAMELIST /ATORB/ NCOEF,EIGVAL,OCCNO
      NCOEF=MSTO
      DO 30 IORB=1,NORB
      EIGVAL(1)=EVAL(IORB,1)
      EIGVAL(2)=EVAL(IORB,ISPIN)
      OCCNO(1)=OCC(IORB,1)
      OCCNO(2)=OCC(IORB,ISPIN)
      WRITE (IP,ATORB)
      JSTO=0
      DO 20 ISTO=1,NSTO
      LI=L(ISTO)
      LLIP1=2*LI+1
      DO 10 MLI=1,LLIP1
      ML=MLI-LI-1
      JSTO=JSTO+1
      COEF(1)=STEVEC(JSTO,IORB,1)
      COEF(2)=STEVEC(JSTO,IORB,ISPIN)
      WRITE (IP,1000) ISTO,ML,COEF
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
      RETURN
 1000 FORMAT(2I3,2F10.6)
      END
