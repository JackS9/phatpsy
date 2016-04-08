      SUBROUTINE GENVAB(STVAB,N,L,ETA,RAB,LV,MV,VEXP,OMEGAI,OMEGAJ,A,B,
     X                  SCRPLM,INDEX,ANORM,YLMNRM,D,CGC,QREALY)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
C-----------------------------------------------------------------------
C
C     GENVAB...
C
C        THIS ROUTINE GENERATES ALL THE MATRIX ELEMENTS BETWEEN STO'S
C     ON A SINGLE CENTER FOR A POTENTIAL ON ANOTHER CENTER OF THE
C     FOLLOWING FORM:
C
C        V(R) = SUM(K) C(K)*V(K;R)
C
C        LV(K).GE.0
C                   V(K;R) = R**(LV(K)-1)*EXP(-VEXP(K)*R)
C        LV(K).LT.0
C                   V(K;R) = R**-LV(K)*EXP(-VEXP(K)*R)*Y(-LV(K),MV(K))
C
C     WHERE
C
C        STVAB(IJ,K) = <I/V(K;R)/J>
C
C        IJ = MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
C
C        !I> = ANORM(I)*RA**(N(I)-1)*EXP(-ETA(I)*RA)*Y(L(I),M(I))
C
C     VARIABLE DEFINITIONS:
C
C        RAB............ DISTANCE BETWEEN THE TWO CENTERS.
C        NVTERM......... THE NUMBER OF TERMS IN THE POTENTIAL SUM.
C        OMEGAI(*,*).... AN ARRAY CONTAINING THE EXPANSION COEFFICIENTS
C                        OF THE I-TH STO IN ELLIPTICAL COORDINATES.
C        OMEGAJ(*,*).... SIMILAR TO OMEGAI(*,*) BUT FOR THE J-TH STO.
C        A(*)........... AN ARRAY FOR HOLDING THE A-FUNCTION VALUES
C                        USED IN THE EXPANSION OF AN STO.
C        B(*)........... SIMILAR TO A(*) BUT FOR THE B-FUNCTIONS.
C        INDEX(N)....... =N*(N-1)/2, INDEX ARRAY FOR SYMMETRY PACKING.
C        SCRPLM(*)...... A SCRATCH ARRAY FOR 'OGEN' USED TO HOLD PLM'S.
C        EXPMIN......... MINIMUM EXPONENT ALLOWED (TO AVOID UNDERFLOW).
C        NNMXM1......... =2*NMAX-1, WHERE NMAX IS THE MAXIMUM VALUE OF
C                        N ENCOUNTERED.
C        M2STO.......... =MSTO*(MSTO+1)/2, WHERE MSTO IS THE NUMBER
C                        OF STO'S (INCLUDING ML-VALUES).
C        NABDIM......... =4*NMAX-2, DIMESNION OF A(*) AND B(*).
C        LMXP1.......... =LMAX+1, WHERE LMAX IS THE MAXIMUM VALUE OF L
C                        ENCOUNTERED.
C        NYLM........... =LLMXP1*(LLXMP1+1)/2, DIMENSION OF YLMNRM(*).
C                        WHERE LLMXP1=2*LMAX+1.
C        NSTO........... NUMBER OF STO'S ON THIS ATOM (NOT INCLUDING
C                        DIFFERENT ML-VALUES).
C        FACT(N)........ =(N-1)-FACTORIAL.
C        BINOM(INDEX(N)+M)... =BINOMIAL COEFFICIENT OF X**N*Y**M/X*Y.
C        YLMNRM(INDEX(L+1)+/M/+1)... =NORMALIZATION CONSTANT FOR Y(L,M).
C        D(L3MX,2)...... D-COEFFICIENTS.
C        CGC(NCGC)...... CLEBSCH-GORDON COEFFICIENTS.
C
C     ROUTINES CALLED:  OGEN, ANMBNM, ABSUM, DCOEF, CGCOEF;
C                       IABS, DSQRT, DEXP, DABS, MOD, MAX0, MIN0
C
C     COMMON USAGE:
C
C        /PARMS/  USES - IPARM(19)(=NNMXM1),IPARM(21)(=LMXP1),
C                        IPARM(27)(=M2STO), IPARM(30)(=NABDIM),
C                        IPARM(32)(=NSTO),  IPARM(34)(=NYLM)
C                        APARM(1)(=EXPMIN)
C
C
C        /BPARMS/ USES - NVTERM
C
C        /TABLES/ USES - FACT(*), BINOM(*), ZERO, HALF, ROOTPI,
C                        HALFPI(2)(=PI)
C
C     RESTRICTIONS:
C
C        THE FACTORIALS, FACT(*), HAVE TO BE DEFINED UP TO 4*LMAX+1.
C        THE DIMENSION HERE IS SUFFICIENT FOR LMAX=5. - FACT(22)
C
C        THE BINOMIAL COEFFICIENT ARRAY, BINOM(*), MUST BE DIMENSIONED
C        UP TO NNMXM1*(NNMXM1+1)/2.  THE DIMENSION HERE IS SUF-
C        FICIENT FOR NMAX=7. - BINOM(91)
C
C.......................................................................
C
C     LAST REVISION:  OCTOBER 12, 1977
C
C     WRITTEN BY:  JACK A. SMITH
C                  QUANTUM THEORY PROJECT
C                  UNIVERSITY OF FLORIDA
C                  GAINESVILLE, FLORIDA
C
C     REFERENCE:   COMPUTATION METHODS OF QUANTUM CHEMISTRY. PT I.
C                  BY FRANK HARRIS (UNIV. OF UTAH).
C
C     SUBORDINATE ROUTINES:  GENCGC, NDXCGC, CGCOEF,CHKCGC,
C                            ORDER, GENDC, DCOEF, NDXD,
C                            OGEN, ANMBNM, ASCALE, BSCALE, ABSUM
C
C-----------------------------------------------------------------------
      COMMON /PARMS/ APARM(20),IPARM(50),QPARM(50)
      EQUIVALENCE (APARM(1),EXPMIN)
      EQUIVALENCE (IPARM(19),NNMXM1), (IPARM(21),LMXP1),
     X            (IPARM(24),NCGC),
     X            (IPARM(27),M2STO),  (IPARM(28),L3MX),
     X            (IPARM(30),NABDIM),
     X            (IPARM(32),NSTO),   (IPARM(34),NYLM)
      COMMON /BPARMS/ IBPARM(13),NVTERM
      COMMON /TABLES/ REALS(10),FACT(22),FFAC(19),BINOM(91),HALFPI(8),
     X                ZERO,HALF,ROOT(10),ROOTPI,CONST(10),CONVRT(10)
      EQUIVALENCE (HALFPI(2),PI)
      DIMENSION STVAB(M2STO,NVTERM),N(NSTO),L(NSTO),ETA(NSTO),
     X          LV(NVTERM),MV(NVTERM),
     X          OMEGAI(NNMXM1,NNMXM1),OMEGAJ(NNMXM1,NNMXM1),
     X          A(NABDIM),B(NABDIM),SCRPLM(LMXP1),INDEX(NNMXM1),
     X          VEXP(NVTERM),YLMNRM(NYLM),ANORM(NSTO),D(L3MX,2),
     X          CGC(NCGC)
	 
      RHALF=HALF*RAB
      CNORM=PI*ROOTPI*RAB*RAB
	  
      IJ=0
	  
      DO 100 I=1,NSTO
      ETAI=ETA(I)
      NI=N(I)
      LI=L(I)
      LLIP1=2*LI+1
	  
      DO 90 MIPLP1=1,LLIP1
      MI=MIPLP1-LI-1
	  
      DO 80 J=1,I
      ETAJ=ETA(J)
      RETA=RHALF*(ETAI+ETAJ)
      VNORM=ANORM(I)*ANORM(J)*CNORM
	  
      NJ=N(J)
      LJ=L(J)
	  NA=NI+NJ-1
	  
      LDIFP1=IABS(LI-LJ)+1
      LSUMP1=LI+LJ+1
      LLJP1=2*LJ+1
      MJHI=LLJP1
      IF (I.EQ.J) MJHI=MIPLP1
	  
      DO 70 MJPLP1=1,MJHI
      MJ=MJPLP1-LJ-1
      MDIF=IABS(MI-MJ)
      MSUM=IABS(MI+MJ)
	  
      MSIGN=1
      IF (MSUM-MDIF) 20,10,30
   10 IF (MI+MJ) 20,30,30
   20 MSIGN=-1
   30 CONTINUE
      MADIF=MSIGN*MDIF
      MASUM=MSIGN*MSUM
	  
      ILAMN=MAX0(LDIFP1,MIN0(MSUM,MDIF)+1)
      ILAMN=ILAMN+MOD(ILAMN+LSUMP1,2)
	  ILAMX=LSUMP1
	  
      IJ=IJ+1
	  
      DO 60 K=1,NVTERM
      STVAB(IJ,K)=ZERO
	  
      RVEXPK=RHALF*VEXP(K)
      RETAI=RETA+RVEXPK
      RETAJ=RETA-RVEXPK
      EXPON=DABS(RETAJ)-RETAI
      IF (EXPON.LT.EXPMIN) GO TO 60
	  
      VNORMK=VNORM*DEXP(EXPON)
	  
      LB=LV(K)
      IF (LB.LT.0) GO TO 35
	  
      NB=LB
      LB=0
      MB=0
      GO TO 37
	  
   35 LB=-LB
      NB=LB+1
      MB=MV(K)
	  
   37 CONTINUE
   
C...  Up N by 1 and use the Nuclear Attraction form of ABSUM
C     (Is the 2*ROOTPI Norm factor valid for NB>0 ??)

	  NB=NB+1

      L2B=INDEX(LB+1)
      SUMLA=ZERO
	  
      DO 50 ILA=ILAMN,ILAMX,2
      IF (ILA.GT.ILAMX) GO TO 50
      LA=ILA-1
	  L2A=INDEX(LA+1)
      LMIN=MIN0(LB,LA)
      IMLMX=2*LMIN+1
	  
      DO 40 IML=1,IMLMX
      ML=IML-LMIN-1
      IABSML=IABS(ML)
	  MLP1=IABSML+1
	  
      LIMI=NA-IABSML
      LIMJ=NB-IABSML
      LIMIJ=LIMI+LIMJ
	  
      CALL OGEN( NA,LA,ML,RHALF,OMEGAI,NNMXM1,BINOM,FACT,INDEX,SCRPLM)
      CALL OGEN(-NB,LB,ML,RHALF,OMEGAJ,NNMXM1,BINOM,FACT,INDEX,SCRPLM)
      CALL ANMBNM(A,RETAI,B,RETAJ,LIMIJ,IABSML,NABDIM)
      CALL ABSUM(SUM,A,B,OMEGAI,OMEGAJ,LIMI,LIMJ,1)
	  
      IF (ML.NE.0) SUM=HALF*SUM
	  
      SUM=SUM*DCOEF(D,LB,MB,ML)*YLMNRM(L2A+MLP1)*YLMNRM(L2B+MLP1)
	  
      SUMLA=SUMLA+SUM*DCOEF(D,LA,MADIF,ML)
     X               *CGCOEF(CGC,MADIF,MI,MJ,LA,LI,LJ,QREALY)
	 
      IF (MASUM.EQ.MADIF) GO TO 40
      SUMLA=SUMLA+SUM*DCOEF(D,LA,MASUM,ML)
     X               *CGCOEF(CGC,MASUM,MI,MJ,LA,LI,LJ,QREALY)
	 
   40 CONTINUE
   
   50 CONTINUE
   
      STVAB(IJ,K)=VNORMK*SUMLA
   60 CONTINUE
   
   70 CONTINUE
   80 CONTINUE
   
   90 CONTINUE
  100 CONTINUE
  
      RETURN
      END
