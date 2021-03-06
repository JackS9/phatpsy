      SUBROUTINE OUTPUT (MATRIX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,
     &                   NCTL)
C.......................................................................
C
C OUTPUT PRINTS A REAL*8 MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C
C AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        MATRIX(*,*).........MATRIX TO BE OUTPUT
C
C        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
C
C        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
C
C        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
C
C        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
C
C        ROWDIM..............ROW DIMENSION OF MATRIX(*,*)
C
C        COLDIM..............COLUMN DIMENSION OF MATRIX(*,*)
C
C        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE
C                                                   2 FOR DOUBLE SPACE
C                                                   3 FOR TRIPLE SPACE
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
C
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P5D24.15 FORMAT FOR
C
C THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
C
C FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
C
C COLUMNS.
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C
C REVISED:  FEBRUARY 26, 1971
C
C.......................................................................
      INTEGER*4 ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      REAL*8 MATRIX(ROWDIM,COLDIM)
      CHARACTER*1 COLUMN(6),ASA(3),BLANK,CTL
      DATA COLUMN/'C','O','L','U','M','N'/,ASA/' ','0','-'/,BLANK/' '/
	  DATA KCOL/5/
      CTL = BLANK
      IF ((NCTL.LE.3).AND.(NCTL.GT.0)) CTL = ASA(NCTL)
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
      LAST = MIN0(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
      WRITE (6,1000) (COLUMN,I,I = BEGIN,LAST)
      DO 1 K = ROWLOW,ROWHI
      DO 4 I=BEGIN,LAST
      IF (MATRIX(K,I).NE.0.D00) GO TO 5
    4 CONTINUE
      GO TO 1
    5 WRITE (6,2000) CTL,K,(MATRIX(K,I), I = BEGIN,LAST)
    1 CONTINUE
    2 LAST = MIN0(LAST+KCOL,COLHI)
    3 RETURN
 1000 FORMAT ('-',12X,5(7X,6A1,I4,7X))
 2000 FORMAT (A1,' ROW',I4,2X,5F24.15)
      END
