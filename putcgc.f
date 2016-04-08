      SUBROUTINE PUTCGC(CGC,NCGC,LMAX,QREALY,IW)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
C-----------------------------------------------------------------------
C
C     PUTCGC...
C
C        THIS ROUTINE WRITES OUT ALL THE VALID CLEBSCH-GORDON COEF-
C     FICIENTS GENERATED BY 'GENCGC'.
C
C     WHERE,
C
C        LMAX...MAXIMUM L-VALUE.
C        QREALY. =T --> COEFFICIENTS CORRESPOND TO REAL Y(L,M)'S.
C                =F --> COEFFICIENTS CORRESPOND TO COMPLEX Y(L,M)'S.
C        IW.....FORTRAN I/O UNIT FOR WRITING.
C
C     ROUTINES CALLED:  CGCOEF; IABS
C
C-----------------------------------------------------------------------
      DIMENSION CGC(NCGC)
      WRITE (IW,1000)
      LMXP1=LMAX+1
      DO 20 L1P1=1,LMXP1
      L1=L1P1-1
      LL1P1=2*L1+1
      DO 20 L2P1=1,LMXP1
      L2=L2P1-1
      LL2P1=2*L2+1
      LLMN=IABS(L1-L2)+1
      LLMX=L1+L2+1
      DO 20 LLP1=LLMN,LLMX
      LL=LLP1-1
      DO 20 M1PLP1=1,LL1P1
      M1=M1PLP1-L1P1
      DO 20 M2PLP1=1,LL2P1
      M2=M2PLP1-L2P1
      MM=M1+M2
      QFIRST=.TRUE.
   10 IF (M1*M2.LT.0) MM=-IABS(MM)
      IF (IABS(MM).GT.LL) GO TO 20
      CG=CGCOEF(CGC,MM,M1,M2,LL,L1,L2,QREALY)
      WRITE (IW,2000) MM,M1,M2,LL,L1,L2,CG
      IF ((.NOT.QFIRST).OR.(.NOT.QREALY).OR.(M1*M2.EQ.0)) GO TO 20
      QFIRST=.FALSE.
      MM=IABS(M1-M2)
      GO TO 10
   20 CONTINUE
      RETURN
 1000 FORMAT('- MM M1 M2 LL L1 L2  CLEBSCH-GORDON COEFFICIENT'/
     X       '+ __ __ __ __ __ __  __________________________'/)
 2000 FORMAT(' ',6I3,F17.8)
      END
