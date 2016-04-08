C
C     ------------------------------------------------------------------
C
      SUBROUTINE TRED3(N,NV,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,IZ,JK,NV
      REAL*8 A(NV),D(N),E(N),E2(N)
      REAL*8 F,G,H,HH,SCALE
      REAL*8 DSQRT,DABS,DSIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX;
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT;
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C                FOR I=N STEP -1 UNTIL 1 DO --
      DO  90 II = 1, N
      I = N + 1 - II
      L = I - 1
      IZ = (I * L) / 2
      H = 0.0D0
      SCALE = 0.0D0
      IF (L .LT. 1) GO TO 20
C                SCALE ROW (ALGOL TOL THEN NOT NEEDED)
      DO 10 K = 1, L
      IZ = IZ + 1
      D(K) = A(IZ)
      SCALE = SCALE + DABS(D(K))
   10 CONTINUE
C
      IF (SCALE .NE. 0.0D0) GO TO 30
   20 E(I) = 0.0D0
      E2(I) = 0.0D0
      GO TO 80
C
   30 DO 40 K = 1, L
      D(K) = D(K) / SCALE
      H = H + D(K) * D(K)
   40 CONTINUE
C
      E2(I) = SCALE * SCALE * H
      F = D(L)
      G = -DSIGN(DSQRT(H),F)
      E(I) = SCALE * G
      H = H - F * G
      D(L) = F - G
      A(IZ) = SCALE * D(L)
      IF (L .EQ. 1) GO TO 80
      F = 0.0D0
C
      DO 60 J = 1, L
      G = 0.0D0
      JK = (J * (J-1)) / 2
C                FORM ELEMENT OF A*U
      DO 50 K = 1, L
      JK = JK + 1
      IF (K .GT. J) JK = JK + K - 2
      G = G + A(JK) * D(K)
   50 CONTINUE
C                FORM ELEMENT OF P
      E(J) = G / H
      F = F + E(J) * D(J)
   60 CONTINUE
C
      HH = F / (H + H)
      JK = 0
C                FORM REDUCED A
      DO 70 J = 1, L
      F = D(J)
      G = E(J) - HH * F
      E(J) = G
C
      DO 70 K = 1, J
      JK = JK + 1
      A(JK) = A(JK) - F * E(K) - G * D(K)
   70 CONTINUE
C
   80 D(I) = A(IZ+1)
      A(IZ+1) = SCALE * DSQRT(H)
   90 CONTINUE
C
      RETURN
C                LAST CARD OF TRED3
      END