      SUBROUTINE FOKMAT(FC,FX,D,TEIBUF)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
C-----------------------------------------------------------------------
C
C     FOKMAT...
C
C        THIS ROUTINE CONSTRUCTS THE TWO-ELECTRON COULOMB AND EXCHANGE
C     PARTS OF THE FOCK MATRIX, SEPARATELY OR SUMMED, BY READING IN AND
C     LOOPING OVER THE NONZERO TWO-ELECTRON INTEGRALS (ALONG WITH THE
C     POINTERS AND WEIGHTS AS DESCRIBED IN 'TWOINT'), MULTIPLYING IT
C     BY THE PROPER DENSITY MATRIX ELEMENTS AND WEIGHT FACTORS, AND
C     ADDING THE RESULTS TO THE PROPER ELEMENTS OF THE FOCK MATRIX
C     WHICH ARE TO BE ADDED TO THE ONE-ELECTRON HAMILTONIAN MATRIX.
C
C     VARIABLE DEFINITIONS:
C
C        FC(*,*)...... THE COULOMB PART OF THE FOCK MATRIX.
C        FX(*,*)...... THE EXCHANGE PART OF THE FOCK MATRIX.
C        D(*,*)....... THE DENSITY MATRIX (FOCK-DIRAC) OVER STO BASIS.
C        ISPIN = 1 ... CLOSED SHELL CASE.
C              = 2 ... OPEN SHELL CASE.
C        TEIBUF.... THE TWO-ELECTRON INTEGRAL BUFFER.
C        LENBUF.... THE LENGTH (DIMENSION) OF TEIBUF.
C        TEI....... THE VALUE OF THE TWO-ELECTRON INTEGRAL JUST READ.
C        M2STO..... MSTO*(MSTO+1)/2, DIMENSION OF F(*) AND D(*) WHERE
C                   MSTO IS THE NUMBER OF STO'S (INCLUDING ML-VALUES).
C
C        NOTE:  THE COULOMB AND EXCHANGE ARE SUMMED IF FC(*,*)=FX(*,*).
C
C     ROUTINES CALLED:  DERASE, TEINDX(GETTEI)
C
C     COMMON USAGE:
C
C        /PARMS/  USES - IPARM(27)(=M2STO),  IPARM(33)(=ISPIN)
C
C        /NDXTEI/ USES - IAC,KAC,IBD,KBD,IAD,KAD,IBC,KBC,
C                        IAB,KAB,ICD,KCD
C
C        /IODATA/ USES - LENBUF
C
C-----------------------------------------------------------------------
      COMMON /PARMS/ APARM(20),IPARM(50),QPARM(50)
      EQUIVALENCE (IPARM(27),M2STO),  (IPARM(33),ISPIN)
      COMMON /NDXTEI/ IAC,KAC,IBD,KBD,IAD,KAD,IBC,KBC,IAB,KAB,ICD,KCD
      DIMENSION FC(M2STO,ISPIN), FX(M2STO,ISPIN),
     X          D(M2STO,ISPIN),  TEIBUF(LENBUF)
      COMMON /IODATA/ IUNIT(20),LENBUF
      CALL DERASE(FC,M2STO*ISPIN)
      CALL DERASE(FX,M2STO*ISPIN)
      CALL TEINDX(0,TEIBUF)
   10 CALL GETTEI(TEI,TEIBUF)
      IF (IAC.EQ.0) GO TO 30
      FAC=KAC*(D(IBD,1)+D(IBD,ISPIN))*TEI
      FBD=KBD*(D(IAC,1)+D(IAC,ISPIN))*TEI
      DO 20 ISP=1,ISPIN
      FC(IAC,ISP)=FC(IAC,ISP)+FAC
      FC(IBD,ISP)=FC(IBD,ISP)+FBD
      FAD=-KAD*D(IBC,ISP)*TEI
      FX(IAD,ISP)=FX(IAD,ISP)+FAD
      FBC=-KBC*D(IAD,ISP)*TEI
      FX(IBC,ISP)=FX(IBC,ISP)+FBC
      FAB=-KAB*D(ICD,ISP)*TEI
      FX(IAB,ISP)=FX(IAB,ISP)+FAB
      FCD=-KCD*D(IAB,ISP)*TEI
      FX(ICD,ISP)=FX(ICD,ISP)+FCD
   20 CONTINUE
      GO TO 10
   30 RETURN
      END