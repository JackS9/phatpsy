      SUBROUTINE GENROT(ALPHA,BETA,GAMMA,ROTATE)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C
C     GENROT...
C
C        THIS ROUTINE GENERATES THE ROTATION MATRIX FOR TRANSFORMING
C     A SET OF CARTESIAN COORDINATES GIVEN THE EULER ANGLES.
C
C     VARIABLE DEFINITIONS:
C
C        ALPHA....... COUTERCLOCKWISE ROTATION ABOUT THE Z-AXIS.
C        BETA........ COUNTERCLOCKWISE ROTATION ABOUT THE NEW
C                     X-AXIS (AFTER ROTATION ALPHA).
C        GAMMA....... COUNTERCLOCKWISE ROTATION ABOUT THE NEW
C                     Z-AXIS (AFTER ROTATIONS ALPHA AND BETA).
C        ROTATE(*)... ROTATION MATRIX.
C
C     ROUTINES CALLED:  DCOS, DSIN
C
C-----------------------------------------------------------------------
      DIMENSION ROTATE(3,3)
      COSA=DCOS(ALPHA)
      SINA=DSIN(ALPHA)
      COSB=DCOS(BETA)
      SINB=DSIN(BETA)
      COSG=DCOS(GAMMA)
      SING=DSIN(GAMMA)
      ROTATE(1,1)=COSG*COSA-COSB*SINA*SING
      ROTATE(1,2)=COSG*SINA+COSB*COSA*SING
      ROTATE(1,3)=SING*SINB
      ROTATE(2,1)=-SING*COSA-COSB*SINA*COSG
      ROTATE(2,2)=-SING*SINA+COSB*COSA*COSG
      ROTATE(2,3)=COSG*SINB
      ROTATE(3,1)=SINB*SINA
      ROTATE(3,2)=-SINB*COSA
      ROTATE(3,3)=COSB
      RETURN
      END