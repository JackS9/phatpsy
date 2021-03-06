      SUBROUTINE CHRGBO(DENS,OVLP,POPUL,TOTPOP,NBAS,N2BAS)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
C-----------------------------------------------------------------------
C
C     CHRGBO...
C
C        THIS ROUTINE CONSTRUCTS THE NET-CHARGE BOND-ORDER MATRIX AND
C     THE GROSS ORBITAL POPULATIONS FROM A GENERAL DENSITY (FOCK-DIRAC)
C     MATRIX AND A  CORRESPONDING OVERLAP (METRIC) MATRIX.
C
C     DEFINITIONS:
C
C        DENS(*)..... SYMMETRY PACKED DENSITY MATRIX.  EQUAL TO THE
C                     NET-CHARGE BOND-ORDER MATRIX ON RETURN.
C        OVLP(*)..... SYMMETRY PACKED OVERLAP MATRIX FOR THE BASIS.
C        POPUL(*).... GROSS ORBITAL (BASIS) POPULATIONS.
C                     (THE OVERLAP POPULATIONS ARE WEIGHTED BY THE
C                      CORRESPONDING NET-CHARGES).
C        TOTPOP...... TOTAL POPULATION, SUM OF POPUL(*).
C        NBAS........ ORDER OF THE BASIS.
C        N2BAS....... =NBAS*(NBAS+1)/2, DIMENSION OF DENS(*) AND OVLP(*)
C
C     ENTRY:  MULPOP...
C
C        THIS ENTRY POINT RETURNS THE GROSS POPULATIONS AS DEFINED
C     BY MULLIKEN (WEIGHT = 1/2).
C
C-----------------------------------------------------------------------
      DIMENSION DENS(N2BAS),OVLP(N2BAS),POPUL(NBAS)
      DATA ZERO/0.D0/,HALF/0.5D0/,TWO/2.D0/
      QMULL=.FALSE.
      GO TO 10
C
C...
      ENTRY MULPOP(DENS,OVLP,POPUL,TOTPOP,NBAS,N2BAS)
      QMULL=.TRUE.
   10 CONTINUE
      IF (NBAS.LT.2) GO TO 30
      DO 20 I=2,NBAS
      II=(I*(I-1))/2
      IM1=I-1
      DO 20 J=1,IM1
      IJ=II+J
      DENS(IJ)=DENS(IJ)*OVLP(IJ)*TWO
   20 CONTINUE
   30 CONTINUE
      TOTPOP=ZERO
      DO 50 I=1,NBAS
      II=(I*(I-1))/2
      III=II+I
      POPI=DENS(III)
      POP=POPI
      DO 40 J=1,NBAS
      IF (J.EQ.I) GO TO 40
      JJ=(J*(J-1))/2
      JJJ=JJ+J
      POPJ=DENS(JJJ)
      IJ=II+J
      IF (J.GT.I) IJ=JJ+I
      WGT=POPI/(POPI+POPJ)
      IF (QMULL) WGT=HALF
      POP=POP+DENS(IJ)*WGT
   40 CONTINUE
      POPUL(I)=POP
      TOTPOP=TOTPOP+POP
   50 CONTINUE
      RETURN
      END
