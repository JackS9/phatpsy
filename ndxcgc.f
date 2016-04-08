      FUNCTION NDXCGC(MM,M1,M2,LL,L1,L2)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
C-----------------------------------------------------------------------
C
C     NDXCGC...
C
C        THIS FUNCTION RETURNS THE COMPACT SINGLE INDEX USED FOR
C     STORING AND RETRIEVING THE CLEBSCH-GORDON COEFFICIENTS.  THE
C     INDEX IS SET TO 1 IF THE COEFFICIENT IS FOUND TO BE ZERO BY
C     'CHKCGC'.  THE INDEXING ASSUMES THAT LL GE L1 GE L2, THAT
C     MM GE M1 FOR LL=L1 AND THAT M1 GE M2 FOR L1=L2, BUT THIS IS
C     CHECKED FOR AND TEMPORARILY CORRECTED FOR BY 'ORDER' (AND ENTRY
C     POINT 'REORDR') IF THE CONDITION IS NOT IS SATISFIED.  NDXCGC
C     IS SET TO 2 FOR THE SPECIAL CASE OF C(M1,M1,0;L1,L1,0) AND
C     TO 3 IF THE INDEX IS GREATER THAN NCGC.
C
C     ENTRY:  NDXCG...
C
C        THIS ENTRY POINT MERELY SKIPS THE CHECK FOR ZERO COEFFICIENTS.
C
C     ENTRY:  NDXMAX...
C
C        THIS ENTRY POINT RETURNS THE MAXIMUM VALUE OF THE INDEX FOR A
C     GIVEN MAXIMUM VALUE OF L (LMAX).
C
C     VARIABLE DEFINITIONS:  (SEE GENCGC)
C
C     ROUTINES CALLED:  CHKCGC, ORDER(REORDR); IABS
C
C     COMMON USAGE:
C
C        /PARMS/   USES - IPARM(24)(=NCGC)
C
C-----------------------------------------------------------------------
      COMMON /PARMS/ APARM(20),IPARM(50),QPARM(50)
      EQUIVALENCE (IPARM(24),NCGC)
      NDXCGC=1
      CALL CHKCGC(MM,M1,M2,LL,L1,L2,QZERO)
      IF (QZERO) RETURN
C
C -->
      ENTRY NDXCG(MM,M1,M2,LL,L1,L2)
      QMXNDX=.FALSE.
      CALL ORDER(MM,M1,M2,LL,L1,L2,QREORD)
      IF ((M2.NE.0).OR.(L2.NE.0)) GO TO 10
      NDXCGC=2
      GO TO 20
	  
   10 CONTINUE
      IABM1=IABS(M1)
      IABM2=IABS(M2)
      N=L2+3
      IF (MM.EQ.(IABM1+IABM2)) N=2
      NDXCGC=L1*(L1+1)*(2*L1+1)*(L1*L1+L1+3)/15
     X      +IABM1*(L1+1)*(2*L1*L1+L1+IABM1+3)/3
     X      +L2*(L2+1)*(2*L2+1)/3
     X      +2*IABM2*(L2+1)
     X      +(LL-L1-L2)/2
     X      +N
      IF (NDXCGC.GT.NCGC .AND. NCGC.GT.0) NDXCGC=3
   20 IF (QREORD) CALL REORDR(MM,M1,M2,LL,L1,L2)
      RETURN
      END
C
C -->
      FUNCTION NDXMAX(LMAX)
      N=LMAX+3
      IF (LMAX.EQ.0) N=2
      NDXMAX=LMAX*(LMAX+1)*(2*LMAX+1)*(LMAX*LMAX+LMAX+3)/15
     X      +LMAX*(LMAX+1)*(2*LMAX*LMAX+LMAX+LMAX+3)/3
     X      +LMAX*(LMAX+1)*(2*LMAX+1)/3
     X      +2*LMAX*(LMAX+1)
     X      +N
      RETURN
      END
