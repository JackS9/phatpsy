      SUBROUTINE NORHES(H,N,FAC,NN)
      REAL*8 A,H(N,N),FAC(N)
      N1=N-1
      NN=N
      DO 20 J=2,N
      A=H(J,J-1)
      IF (A.EQ.0.D0) GO TO 30
      FAC(J)=A
      H(J,J-1)=1.D0
      DO 10 K=2,J
10    H(K-1,J)=H(K-1,J)*A
      IF (J.EQ.N) RETURN
      H(J+1,J)=H(J+1,J)*A
      A=1.D0/A
      DO 20 K=J,N1
20    H(J,K+1)=H(J,K+1)*A
30    NN=J-1
      RETURN
      END
