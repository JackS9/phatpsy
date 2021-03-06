      SUBROUTINE DENMAT(DENS,EVEC,OCC,NBAS,N2BAS,NVEC)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C
C     DENMAT...
C
C        THIS ROUTINE CONSTRUCTS A ONE-ELECTRON DENSITY MATRIX FROM AN
C     EIGENVECTOR ARRAY OVER SOME BASIS AND A SET OF OCCUPATION NUMBERS.
C
C     VARIABLE DEFINITIONS:
C
C        EVEC(*,*)... EIGENVECTOR ARRAY (COLUMN VECTORS).
C        OCC(*)...... OCCUPATION NUMBERS FOR EACH VECTOR.
C        DENS(*)..... ONE-ELECTRON DENSITY MATRIX OVER SAME BASIS.
C        NBAS........ ORDER OF BASIS (NUMBER OF ROWS).
C        N2BAS....... =NBAS*(NBAS+1)/2, DIMENSION OF DENS(*).
C        NVEC........ NUMBER OF VECTORS (COLUMNS).
C
C     ROUTINES CALLED:  DABS
C
C-----------------------------------------------------------------------
      DIMENSION DENS(N2BAS),EVEC(NBAS,NVEC),OCC(NVEC)
      DATA ZERO/0.0D0/,SIGNIF/1.D-8/
      IJ=0
      DO 30 I=1,NBAS
      DO 20 J=1,I
      IJ=IJ+1
      SUMOCC=ZERO
      DO 10 K=1,NVEC
      SUMOCC=SUMOCC+EVEC(I,K)*EVEC(J,K)*OCC(K)
   10 CONTINUE
C...  ELIMINATE SMALL ELEMENTS TO HELP PREVENT SYMMETRY MIXING.
      IF (DABS(SUMOCC).LT.SIGNIF) SUMOCC=ZERO
      DENS(IJ)=SUMOCC
   20 CONTINUE
   30 CONTINUE
      RETURN
      END
