      SUBROUTINE EXCITE(OCCNO,NORB,ISPIN)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
C-----------------------------------------------------------------------
C
C     EXCITE...
C
C        THIS ROUTINE GENERATES THE OCCUPATION NUMBERS FOR AN EXCITED
C     (OR IONIZED) STATE FROM THE OCCUPATION NUMBERS OF A GROUND STATE.
C
C     VARIABLE DEFINITIONS:
C
C        OCCNO(*).... GROUND STATE OCCUPATION NUMBERS ON INPUT AND
C                     EXCITED STATE OCCUPATION NUMBERS ON OUTPUT.
C        NORB........ NUMBER OF ORBITALS.
C        ISPIN....... =1 IF RESTRICTED CASE, =2 IF UNRESTRICTED CASE.
C        NCREAT...... NUMBER OF CREATION OPERATIONS.
C        NANNIH...... NUMBER OF ANNIHILATION OPERATIONS.
C        ICREAT(*)... ORBITALS TO BE CREATED (OCCUPIED).
C                     (NEGATIVE FOR BETA SPIN-ORBITAL)
C        IANNIH(*)... ORBITALS TO BE ANNIHILATED (UNOCCUPIED).
C                     (NEGATIVE FOR BETA SPIN-ORBITAL)
C        QTRANS...... =T --> TRANSITION STATE, HALF OCCUPATIONS.
C
C     ROUTINES CALLED:  IABS
C
C     COMMON USAGE:
C
C        /STATE/  USES - NCREAT,NANNIH,ICREAT(*),IANNIH(*),QTRANS
C
C-----------------------------------------------------------------------
      COMMON /STATE/ NCREAT,NANNIH,ICREAT(7),IANNIH(7),QTRANS
      DIMENSION OCCNO(NORB,ISPIN)
      DATA HALF/0.5D0/
      ORBOCC=HALF*ISPIN
      IF (QTRANS) ORBOCC=HALF*ORBOCC
      IF (NCREAT.LE.0) GO TO 20
      DO 10 NCR=1,NCREAT
      ICR=ICREAT(NCR)
      IORB=IABS(ICR)
      IF ((IORB.EQ.0).OR.(IORB.GT.NORB)) GO TO 10
      ISP=1
      IF (ICR.LT.0) ISP=ISPIN
      OCCNO(IORB,ISP)=OCCNO(IORB,ISP)+ORBOCC
   10 CONTINUE
   20 CONTINUE
      IF (NANNIH.LE.0) GO TO 40
      DO 30 NAN=1,NANNIH
      IAN=IANNIH(NAN)
      IORB=IABS(IAN)
      IF ((IORB.EQ.0).OR.(IORB.GT.NORB)) GO TO 30
      ISP=1
      IF (IAN.LT.0) ISP=ISPIN
      OCCNO(IORB,ISP)=OCCNO(IORB,ISP)-ORBOCC
   30 CONTINUE
   40 CONTINUE
      RETURN
      END
