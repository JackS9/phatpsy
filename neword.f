      SUBROUTINE NEWORD(VECTOR,SCRAT,IORDER,NBAS)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C
C     NEWORD...
C
C        THIS ROUTINE REORDERS THE ELEMENTS OF A VECTOR ACCORDING TO
C     AN ORDERING ARRAY.
C
C     DEFINITIONS:
C
C        VECTOR(*)... VECTOR WHOSE ELEMENTS ARE TO BE REORDERED.
C        SCRAT(*).... A SCRATCH ARRY FOR STORING THE VECTOR IN ITS
C                     ORIGINAL ORDER (RETAINS IT ON RETURN ALSO).
C        IORDER(*)... THE ORDERING ARRAY.
C        NBAS........ LENGTH OF THE VECTOR.
C
C     ROUTINES CALLED:  DCOPY
C
C-----------------------------------------------------------------------
      DIMENSION VECTOR(NBAS),SCRAT(NBAS),IORDER(NBAS)
      CALL DCOPY(VECTOR,NBAS,SCRAT)
      DO 10 I=1,NBAS
      J=IORDER(I)
      IF ((J.LE.0).OR.(J.GT.NBAS)) GO TO 10
      VECTOR(I)=SCRAT(J)
   10 CONTINUE
      RETURN
      END
