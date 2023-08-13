      SUBROUTINE GENPLM(PLM,BETA,NPLM,LMAX)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C
C     GENPLM...
C
C        THIS ROUTINE GENERATES ALL THE ASSOCIATED LEGENDRE POLYNOMIALS
C     UP TO A MAXIMUM L-VALUE AND EVALUATES THEM FOR A SPECIFIC ANGLE.
C
C     VARIABLE DEFINITIONS:
C
C        PLM(LM).... P(L,M), WHERE LM=(L*(L+1))/2+ABS(M)+1.
C        BETA....... AZIMUTHAL ANGLE (TILT FROM Z-AXIS).
C        LMAX....... MAXIMUM L-VALUE.
C        NPLM....... =(LMAX+1)*(LMAX+2)/2, DIMENSION OF PLM(*).
C
C     ROUTINES CALLED:  DCOS, DSQRT, DFLOAT, MOD
C
C-----------------------------------------------------------------------
      DIMENSION PLM(NPLM)
      DATA ZERO/0.0D0/,ONE/1.0D0/,SIGNIF/1.0D-10/
      IINDEX(IARG)=(IARG*(IARG-1))/2
      LMXP1=LMAX+1
      COSB=DCOS(BETA)
      PLM(1)=ONE
      PLM(2)=COSB
      IF ((ONE-COSB*COSB).LT.SIGNIF) GO TO 40
      DO 10 LP1=2,LMAX
      L=LP1-1
      RL=L
      LP1M=IINDEX(L+2)+1
      LM=IINDEX(LP1)+1
      LM1M=IINDEX(L)+1
      PLM(LP1M)=((RL+RL+ONE)*COSB*PLM(LM)-RL*PLM(LM1M))/(RL+ONE)
   10 CONTINUE
      FACTOR=ONE/DSQRT(ONE-COSB*COSB)
      DO 30 LP1=2,LMXP1
      L=LP1-1
      DO 20 MP1=1,L
      M=MP1-1
      LM=IINDEX(LP1)+M+1
      LMP1=LM+1
      LM1M=IINDEX(L)+M+1
      TERM=DFLOAT(L-M)*COSB*PLM(LM)
      IF (MP1.LE.L) TERM=TERM-DFLOAT(L+M)*PLM(LM1M)
   20 PLM(LMP1)=FACTOR*TERM
   30 CONTINUE
      GO TO 70
   40 LM=1
      DO 60 LP1=1,LMXP1
      DO 50 MP1=1,LP1
      PLM(LM)=ZERO
   50 LM=LM+1
      L0=IINDEX(LP1)+1
      PLM(L0)=ONE
      IF ((COSB.LT.ZERO).AND.(MOD(LP1,2).EQ.0)) PLM(L0)=-ONE
   60 CONTINUE
   70 CONTINUE
      RETURN
      END