      SUBROUTINE ROTSAB(STOVLP,OVRLAP,LA,NSTOA,MSTOA,LB,NSTOB,MSTOB,D)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
C-----------------------------------------------------------------------
C
C     ROTSAB...
C
C        THIS ROUTINE ROTATES A SET OF OVERLAPS BETWEEN STO'S ON TWO
C     CENTERS WHICH WERE ORIGINALLY CALCULATED FOR THE SPECIAL DIATOMIC
C     CASE (COMMON Z-AXIS).
C
C     VARIABLE DEFINITIONS:
C
C        STOVLP(*,*)..... ROTATED 2-CENTER OVERLAP INTEGRALS.
C        OVRLAP(*,*,*)... DIATOMIC (COMMON Z-AXIS) OVERLAP INTEGRALS.
C        LA(*)........... L-QUANTUM NUMBERS FOR BASIS ON ATOM A.
C        NSTOA........... NUMBER OF STO'S IN ATOM A'S BASIS (NOT
C                         INCLUDING ML-VALUES).
C        MSTOA........... NUMBER OF STO'S (INCLUDING ML-VALUES).
C        LB(*),NSTOB,MSTOB... SIMILARLY FOR ATOM B.
C        D(*)............ D-COEFFICIENTS.
C
C     ROUTINES CALLED:  DCOEF; MIN0
C
C     COMMON USAGE:
C
C        /PARMS/  USES - IPARM(22)(=LLMXP1), IPARM(28)(=L3MX)
C
C        /TABLES/ USES - ZERO
C
C-----------------------------------------------------------------------
      COMMON /PARMS/ APARM(20),IPARM(50),QPARM(50)
      EQUIVALENCE (IPARM(22),LLMXP1),(IPARM(28),L3MX)
      COMMON /TABLES/ REALS(10),FACT(22),FFAC(19),BINOM(91),HALFPI(8),
     X                ZERO,HALF,ROOT(10),ROOTPI,CONST(10),CONVRT(10)
      DIMENSION STOVLP(MSTOA,MSTOB),OVRLAP(LLMXP1,NSTOA,NSTOB),
     X          LA(NSTOA),LB(NSTOB),D(L3MX,2)
	 
      IA=0
	  
      DO 50 I=1,NSTOA
      LAI=LA(I)
      LLAIP1=2*LAI+1
	  
      DO 40 IM=1,LLAIP1
      MAI=IM-LAI-1
      IA=IA+1
      JB=0
	  
      DO 30 J=1,NSTOB
      LBJ=LB(J)
      LLBJP1=2*LBJ+1
      LMIN=MIN0(LAI,LBJ)
      LLMNP1=2*LMIN+1
	  
      DO 20 JM=1,LLBJP1
      MBJ=JM-LBJ-1
      JB=JB+1
      SUM=ZERO
	  
      DO 10 IML=1,LLMNP1
      ML=IML-LMIN-1
      SUM=SUM+DCOEF(D,LAI,MAI,ML)*DCOEF(D,LBJ,MBJ,ML)*OVRLAP(IML,I,J)
   10 CONTINUE
   
      STOVLP(IA,JB)=SUM
   20 CONTINUE   
   30 CONTINUE
   
   40 CONTINUE
   50 CONTINUE
   
      RETURN
      END