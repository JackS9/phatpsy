      SUBROUTINE LOWDIN(OVRLAP,ONTRAN,SCRAT1,NBAS,NDIM,ITYPE,QINDEF)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
C-----------------------------------------------------------------------
C
C     LOWDIN...
C
C        THIS ROUTINE RETURNS ONE OF THE FOLLOWING ORTHONORMAL TRANS-
C     FORMATION MATRICES
C
C        ITYPE             TYPE
C
C                           -1   T    -1
C          0 ) . . . . . U D    U  = S    (INVERSE)
C
C                           -1/2      -1/2
C          1 ) . . . . . U D       = S    (CANONICAL)
C
C                           -1/2 T    -1/2
C          2 ) . . . . . U D    U  = S    (SYMMETRIC)
C
C                           +1/2 T    +1/2
C          3 ) . . . . .   D    U  = S    (CANONICAL)
C
C                           +1/2 T    +1/2
C          4 ) . . . . . U D    U  = S    (SYMMETRIC)
C
C     WHERE U IS THE UNITARY MATRIX WHICH DIAGONALIZES THE OVERLAP
C     MATRIX S
C
C             T
C            U S U = D
C
C     VARIABLE DEFINITIONS:
C
C        OVRLAP(*)...... OVERLAP (METRIC) MATRIX S FOR THE NONORTHOGONAL
C                        BASIS.  THE DIAGONAL MATRIX D ON RETURN.
C        ONTRAN(*,*).... ORTHONORMAL TRANSFORMATION MATRIX.
C        SCRAT1(*,*).... SCRATCH MATRIX USED FOR THE UNITARY MATRIX U.
C        NBAS........... ORDER OF BASIS, COLUMN DIMENSION OF ONTRAN(*).
C        NDIM........... ROW DIMENSION OF ONTRAN(*).
C
C        QINDEF......... =T,IF LINEAR DEPENDENCY (NEGATIVE EIGENVALUE).
C
C     NOTE:  IF ITYPE=1 THEN ONTRAN(*,*) AND SCRAT1(*,*) MAY BE THE SAME.
C
C     NOTE:  NOT RECOMMENDED FOR VERY LARGE MATRICES SINCE
C            JACOBI IS NOT ADEQUATE AND MUCH MORE SCRATCH SPACE
C            WOULD BE REQUIRED TO IMPROVE UPON THIS.
C
C     ROUTINES CALLED:  JACOBI; DSQRT, DABS, MOD
C
C-----------------------------------------------------------------------
      DIMENSION OVRLAP(1),ONTRAN(NDIM,NBAS),SCRAT1(NBAS,NBAS)
C     DATA ZERO/0.0D0/,ONE/1.0D0/,TOLER/1.0D-25/
      DATA ZERO/0.0D0/,ONE/1.0D0/,TOLER/1.0D-50/
      QINDEF=.FALSE.
      CALL JACOBI(NBAS,OVRLAP,SCRAT1,NBAS)
      DO 40 I=1,NBAS
      DO 30 J=1,NBAS
      IF (DABS(SCRAT1(I,J)*TOLER).GT.ONE) WRITE (6,1000)
      IF (MOD(ITYPE,2).NE.0) GO TO 20
      SUM=ZERO
      DO 10 K=1,NBAS
      KK=(K*(K+1))/2
      DIAG=OVRLAP(KK)
C     IF (DIAG.LT.TOLER.AND.ITYPE.GT.0) GO TO 50
      IF (DIAG.LT.TOLER) GO TO 50
      IF (ITYPE.EQ.0) DIAG=ONE/DIAG
      IF (ITYPE.EQ.2) DIAG=ONE/DSQRT(DIAG)
      IF (ITYPE.EQ.4) DIAG=DSQRT(DIAG)
      SUM=SUM+SCRAT1(I,K)*DIAG*SCRAT1(J,K)
   10 CONTINUE
      ONTRAN(I,J)=SUM
      GO TO 30
   20 JJ=(J*(J+1))/2
      DIAG=OVRLAP(JJ)
      IF (DIAG.LT.TOLER) GO TO 50
      DIAG=DSQRT(DIAG)
      IF (ITYPE.EQ.1) ONTRAN(I,J)=SCRAT1(I,J)/DIAG
      IF (ITYPE.EQ.3) ONTRAN(J,I)=DIAG*SCRAT1(I,J)
   30 CONTINUE
   40 CONTINUE
      RETURN
   50 WRITE (6,2000) DIAG
      QINDEF=.TRUE.
      RETURN
 1000 FORMAT('-*** WARNING - ORTHONORMAL TRANSFORMATION UNSTABLE ***')
 2000 FORMAT('-*** LINEAR DEPENDENCY ENCOUNTERED, METRIC EIGENVALUE =',
     X       1PD10.1,' ***')
      END
