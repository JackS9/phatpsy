      SUBROUTINE SCHMID (S,NROWS)
C...........VERSION = 03/11/73/03
C......................................................................
C     THIS SUBROUTINE CARRIES OUT A SCHMIDT ORTHONORMALIZATION OF A
C POSITIVE DEFINITE OVERLAP MATRIX AND RETURNS THE TRANSFORMATION
C OVERWRITTEN ON THE OVERLAP MATRIX.  NO ADDITIONAL ARRAY STORAGE
C BEYOND THAT FOR THE OVERLAP MATRIX IS REQUIRED.  SHOULD THE OVERLAP
C MATRIX PROVE TO BE SINGULAR, THE LINEARLY-DEPENDENT FUNCTIONS WILL BE
C DELETED FROM THE BASIS BY INSERTING A ZERO COLUMN IN THE
C TRANSFORMATION MATRIX AND A WARNING MESSAGE WILL BE PRINTED.  THE
C ARGUMENTS ARE:
C
C S(*).........OVERLAP MATRIX STORED IN PACKED LOWER-TRIANGULAR FORM.
C NROWS........NUMBER OF ROWS IN THE OVERLAP MATRIX.
C
C IN THE COMMENTS BELOW DESCRIBING THE ALGORITHM, THE OVERLAP MATRIX
C S(*,*) IS ASSUMED TO CORRESPOND TO A BASIS PHI(*):
C     S(*,*) = <PHI(*)!PHI(*)>
C AND THE SCHMIDT TRANSFORMATION WILL PRODUCE AN UPPER-TRIANGULAR
C MATRIX T(*,*) SUCH THAT
C     CHI(*) = PHI(*) T(*,*)
C AND
C     <CHI(*)!CHI(*)> = T#(*,*) S(*,*) T(*,*)
C                     = 1(*,*)
C WHERE CHI IS THE ORTHONORMALIZED BASIS, T#(*,*) IS THE HERMITIAN
C ADJOINT OF T(*,*), AND 1(*,*) IS THE UNIT MATRIX.
C THIS VERSION IMPLEMENTS ACCUMULATION IN A HIGHER PRECISION THAN
C THAT OF S(*).  ON THE IBM 360 AND 370 MODELS WITH EXTENDED PRECISON
C FACILITIES, S(*) IS DOUBLE PRECISION.  FOR CDC 6000 AND 7000
C SERIES, S(*) IS SINGLE PRECISION AND THE ACCUMULATION IS DOUBLE
C PRECISION.  THE ONLY CHANGES THAT NEED BE MADE ARE TO THE TYPE
C STATEMENTS, THE DATA STATEMENT, AND THE ARITHMETIC STATEMENT FUNCTION
C DEFINITIONS.
C NOTE: FUNCTION EXTEND EXTENDS THE PRECISION OF ITS ARGUMENT.  FOR
C CDC USE, SPECIFY EXTEND(ARG) = DBLE(ARG).
C.......................................................................
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 S
      DIMENSION S(1)
      DATA ONE/1.D0/, ZERO/0.D0/, TOLER/1.D-12/
C.......................................................................
C THE FOLLOWING STATEMENT FUNCTION(S) DEFINE THE BUILT-IN FUNCTION(S)
C REQUIRED BY THIS ROUTINE.
C.......................................................................
      ZSQRT(ARG) = DSQRT(ARG)
      EXTEND(ARG) = ARG
C.......................................................................
      S(1) = ONE/ZSQRT(EXTEND(S(1)))
      IF (NROWS .LE. 1) RETURN
      KPN = 1
      DO 700 N = 2,NROWS
      KPNP1 = KPN + N
      OVRLAP = ZERO
      NM1 = N - 1
      NK = KPN + NM1
C......................................................................
C GENERATE Q(N,N-1), Q(N,N-2), ..., Q(N,1), OVERWRITING ON S(N,N-1),
C S(N,N-2), ..., S(N,1).  Q(N,K) IS DEFINED BY
C     Q(N,K) = <PHI(N)!CHI(K)>
C            = SUM(J=1,K) S(N,J)*T(J,K)
C IN TERMS OF Q(*,*) AND T(*,*), THE NEW FUNCTION CHI(N)' IS DEFINED BY
C     CHI(N)' = PHI(N) - SUM(J=1,N-1) PHI(J) SUM(K=J,N-1) T(J,K) Q(N,K)
C AND THE OVERLAP OF THIS FUNCTION WITH ITSELF IS GIVEN BY
C     <CHI(N)'!CHI(N)'> = S(N,N) - SUM(K=1,N-1) Q(N,K) Q(N,K)
C SO THAT
C     CHI(N) = CHI(N)'/<CHI(N)'!CHI(N)'>**(1/2)
C WE COMPUTE THIS NORMALIZATION FACTOR AS WE GENERATE Q(N,*).
C......................................................................
      KPK = (N*(N-1))/2
      DO 200 KBACK = 1,NM1
      K = N - KBACK
      KPK = KPK - K
      JK = KPK
      NJ = KPN
      SUM = ZERO
      DO 100 J = 1,K
      NJ = NJ + 1
      JK = JK + 1
  100 SUM = SUM + EXTEND(S(NJ))*EXTEND(S(JK))
      S(NK) = SUM
      OVRLAP = OVRLAP + SUM*SUM
  200 NK = NK - 1
      NN = KPNP1
      OVRLAP = EXTEND(S(NN)) - OVRLAP
C......................................................................
C IF THE OVERLAP <CHI(N)'!CHI(N)'> IS BELOW OUR PRESET TOLERANCE,
C DELETE THE FUNCTION FROM THE BASIS BY INSERTING A ZERO COLUMN IN THE
C TRANSFORMATION MATRIX AND PRINT A WARNING MESSAGE.
C......................................................................
      IF (OVRLAP .GT. TOLER) GO TO 400
      WRITE (6,10000) N,NROWS
      JN = KPN
      DO 300 J = 1,N
      JN = JN + 1
  300 S(JN) = ZERO
      GO TO 700
  400 FACTOR = ONE/ZSQRT(OVRLAP)
C......................................................................
C NOW COMPUTE
C     T(N,N)' = 1
C     T(J,N)' = SUM(K=J,N-1) T(J,K) Q(N,K)       (J .LE. N)
C AND THEN NORMALIZE THE COLUMN OF THE TRANSFORMATION MATRIX ACCORDING
C TO
C     T(*,N) = T(*,N)'/<CHI(N)'!CHI(N)'>**(1/2)
C THE N-TH COLUMN OF T(*,*) OVERLAYS THE N-TH ROW OF S(*,*) IN CORE.
C......................................................................
      JN = KPN
      KPJ = 0
      DO 600 J = 1,NM1
      JN = JN + 1
      NK = KPN + J
      JK = KPJ + J
      SUM = ZERO
      DO 500 K = J,NM1
      SUM = SUM - EXTEND(S(JK))*EXTEND(S(NK))
      NK = NK + 1
  500 JK = JK + K
      S(JN) = SUM*FACTOR
  600 KPJ = KPJ + J
      S(JN+1) = FACTOR
  700 KPN = KPNP1
      RETURN
10000 FORMAT (65H-*** WARNING ***  OVERLAP MATRIX PASSED TO SCHMID SINGU
     &LAR AT ROW,I6,21H.  OVERLAP MATRIX HAS,I6,6H ROWS.)
      END
