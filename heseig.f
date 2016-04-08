      SUBROUTINE HESEIG(A,N,NRT,X,XP,EIG,ETA,NZC,*)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N),X(N),XP(N),EIG(NRT)
      REAL*8 TOL
C	  PARAMETER (TOL = $3920000000000000)
	  PARAMETER (TOL = 2.0D0*16.0D0**(-8))
      NROOT=MIN0(NRT,NZC)
      IF (ETA.LT.0.) GO TO 100
      R=A(1,1)-1.D0
      DO 10 I=2,N
      Y=A(I,I)
      IF (I.LT.N) Y=Y-1.D0
      DO 5 J=2,I
5     Y=Y-DABS(A(J-1,I))
      IF (Y.LT.R) R=Y
10    CONTINUE
11    N1=N-1
      DO 50 M=1,NROOT
      XP(NZC)=1.D0
      NIT=0
12    X(NZC)=R-A(NZC,NZC)
      DO 32 I=2,NZC
      NI=NZC-I+1
      X(NI)=(R-A(NI,NI))*X(NI+1)-A(NI,NZC)
      XP(NI)=X(NI+1)+(R-A(NI,NI))*XP(NI+1)
      IF (I.EQ.2) GO TO 32
      NRI=I-2
      DO 30 J=1,NRI
      X(NI)=X(NI)-A(NI,NI+J)*X(NI+J+1)
30    XP(NI)=XP(NI)-A(NI,NI+J)*XP(NI+J+1)
32    CONTINUE
      IF (M.EQ.1) GO TO 40
      Y=XP(NI)/X(NI)
      DO 35 I=2,M
35    Y=Y-1.D0/(R-EIG(I-1))
      Y=1.D0/Y
      GO TO 42
40    Y=X(NI)/XP(NI)
42    R=R-Y
44    NIT=NIT+1
      IF (DABS(Y/R).LE.TOL) GO TO 45
      IF (NIT.GT.25) RETURN1
      GO TO 12
45    EIG(M)=R
50    R=R+ETA
      RETURN
100   R=A(I,I)+1.D0
      DO 110 I=2,N
      Y=A(I,I)
      IF (I.LT.N) Y=Y+1.D0
      DO 105 J=2,I
105   Y=Y+DABS(A(J-1,I))
110   IF (R.LT.Y) R=Y
      GO TO 11
      END
