      SUBROUTINE HESBRG(A,N,MPR)
C CONVERTS MATRIX IN A TO UPPER HESSENBERG FORM
C             BY A SIMILARITY TRANSFORMATION.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),MPR(N)
      N1=N-1
      DO 50 M=1,N
      PIV=0.
      MPIV=M+1
      DO 30 I=1,N
      KMAX=MIN0(I-1,M)
      IF (KMAX.LT.2) GO TO 12
      DO 10 K=2,KMAX
10    A(I,M)=A(I,M)-A(I,K-1)*A(K,M)
12    IF (M.EQ.N) GO TO 30
      IF (M.EQ.1) GO TO 25
      DO 20 K=M,N1
20    A(I,M)=A(I,M)+A(I,K+1)*A(K+1,M-1)
25    IF (I.LE.M) GO TO 30
      IF (DABS(A(I,M)).LE.DABS(PIV)) GO TO 30
      MPIV=I
      PIV=A(I,M)
30    CONTINUE
      IF (M.EQ.N) RETURN
      MPR(M+1)=MPIV
      IF(MPIV.EQ.M+1) GO TO 40
      DO 33 I=1,N
      TEM=A(MPIV,I)
      A(MPIV,I)=A(M+1,I)
33    A(M+1,I)=TEM
      DO 36 I=1,N
      TEM=A(I,MPIV)
      A(I,MPIV)=A(I,M+1)
36    A(I,M+1)=TEM
40    M2=M+2
      IF (M2.GT.N) GO TO 50
      IF (PIV.NE.0.D0) GO TO 45
      DO 44 I=M2,N
44    A(I,M)=0.
      GO TO 50
45    DO 48 I=M2,N
48    A(I,M)=A(I,M)/PIV
50    CONTINUE
      RETURN
      END
