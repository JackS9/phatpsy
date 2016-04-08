      SUBROUTINE ESCA(EVEC,EVAL,N,L,ETA,NORB,MSTO,NSTO,IW)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C
C     ESCA...
C
C        THIS ROUTINE GENERATES ESCA INTENSITIES USING A PLANE WAVE
C     FOR THE PHOTOELECTRON.  THE FROZEN ORBITAL (KOOPMANS) MODEL IS
C     ASSUMED.  THE MOLECULAR ORBITALS ARE APPROXIMATED BY THE
C     ONE-CENTER PARTS.
C
C     VARIABLE DEFININTIONS:
C
C        EVEC(*,*)....... MOLECULAR ORBITAL COEFFICIENTS (ONE-CENTER
C                         CONTRIBUTIONS ONLY).
C        EVAL(*)......... MOLECULAR ORBITAL ENERGIES.
C        N(*)............ N-QUANTUM NUMBERS IN BASIS.
C        L(*)............ L-QUANTUM NUMBERS IN BASIS.
C        ETA(*).......... SLATER EXPONENTS IN BASIS.
C        NORB............ NUMBER OF MOLECULAR ORBITALS.
C        MSTO............ NUMBER OF BASIS FUNCTIONS (ONE CENTER).
C        NSTO............ NUMBER OF BASIS FUNCTIONS NOT INCLUDING
C                         ML-VALUES.
C        IW.............. FORTRAN I/O UNIT FOR WRITING.
C
C     NOTE:  LIMITED TO 4F-FUNCTIONS.
C
C-----------------------------------------------------------------------
      DIMENSION EVEC(MSTO,NORB),EVAL(NORB),N(NSTO),L(NSTO),ETA(NSTO)
      DATA ZERO/0.0D0/,TWO/2.0D0/,THREE/3.0D0/,FOUR/4.0D0/,FIVE/5.0D0/,
     X     SIX/6.0D0/,SEVEN/7.0D0/,EIGHT/8.0D0/,ANINE/9.0D0/,TEN/1.0D1/
      DATA ESOURC/46.0D0/,SCALE/1.0D3/,EVPERH/27.2107D0/
      EVSRC=ESOURC*EVPERH
      WRITE(IW,1000) EVSRC
 1000 FORMAT(//'      EV        INTENSITY      (SOURCE =',F8.2,')'/
     X         '      __        _________'/)
      DO 190 KORB=1,NORB
      AMP=ZERO
      IF (EVAL(KORB).LT.-ESOURC) GO TO 185
      EKAPPA=DSQRT(TWO*(EVAL(KORB)+ESOURC))
      IMSTO=0
      DO 180 ISTO=1,NSTO
      NI=N(ISTO)
      LI=L(ISTO)
      LLIP1=2*LI+1
      NLI=(NI*(NI-1))/2+LI+1
      IF (NLI.GT.10) GO TO 180
      ETAI=ETA(ISTO)
      GO TO (10,20,30,40,50,60,61,62,63,64),NLI
   10 AMPI=(EIGHT/DSQRT(TWO))*ETAI**(FIVE/TWO)*EKAPPA/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**2
      GO TO 65
   20 AMPI=(EIGHT/DSQRT(SIX))*ETAI**(FIVE/TWO)*EKAPPA*
     X     (THREE*ETAI*ETAI-EKAPPA*EKAPPA)/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**3
      GO TO 65
   30 AMPI=(FOUR*EIGHT/DSQRT(SIX))*ETAI**(SEVEN/TWO)*EKAPPA*EKAPPA/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**3
      GO TO 65
   40 AMPI=(FOUR*EIGHT/DSQRT(FIVE))*ETAI**(ANINE/TWO)*EKAPPA*
     X     (ETAI*ETAI-EKAPPA*EKAPPA)/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**4
      GO TO 65
   50 AMPI=(FOUR*EIGHT/DSQRT(FIVE*ANINE))*ETAI**(SEVEN/TWO)*EKAPPA*
     X     EKAPPA*(FIVE*ETAI*ETAI-EKAPPA*EKAPPA)/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**4
      GO TO 65
   60 AMPI=(FOUR*EIGHT/DSQRT(FIVE))*ETAI**(ANINE/TWO)*EKAPPA**3/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**4
      GO TO 65
   61 AMPI=(FOUR*EIGHT/DSQRT(SEVEN*TEN))*ETAI**(ANINE/TWO)*EKAPPA*
     X     (FIVE*ETAI**4-TEN*(ETAI*EKAPPA)**2+EKAPPA**4)/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**5
      GO TO 65
   62 AMPI=(EIGHT*EIGHT/DSQRT(SEVEN*TEN))*ETAI**(ANINE/TWO)*EKAPPA*
     X     ETAI*EKAPPA*(FIVE*ETAI*ETAI-THREE*EKAPPA*EKAPPA)/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**5
      GO TO 65
   63 AMPI=(EIGHT*EIGHT/DSQRT(SEVEN*TEN))*ETAI**(ANINE/TWO)*EKAPPA*
     X     EKAPPA*EKAPPA*(SEVEN*ETAI*ETAI-EKAPPA*EKAPPA)/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**5
      GO TO 65
   64 AMPI=(EIGHT*EIGHT*EIGHT/DSQRT(SEVEN*TEN))*ETAI**(ANINE/TWO)*
     X     ETAI*EKAPPA**4/
     X     (ETAI*ETAI+EKAPPA*EKAPPA)**5
   65 CONTINUE
      DO 170 IML=1,LLIP1
      IMSTO=IMSTO+1
      JMSTO=0
      DO 167 JSTO=1,NSTO
      NJ=N(JSTO)
      LJ=L(JSTO)
      LLJP1=2*LJ+1
      AMPJ=ZERO
      IF (LJ.NE.LI) GO TO 165
      NLJ=(NJ*(NJ-1))/2+LJ+1
      IF (NLJ.GT.10) GO TO 167
      ETAJ=ETA(JSTO)
      GO TO (110,120,130,140,150,160,161,162,163,164),NLJ
  110 AMPJ=(EIGHT/DSQRT(TWO))*ETAJ**(FIVE/TWO)*EKAPPA/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**2
      GO TO 165
  120 AMPJ=(EIGHT/DSQRT(SIX))*ETAJ**(FIVE/TWO)*EKAPPA*
     X     (THREE*ETAJ*ETAJ-EKAPPA*EKAPPA)/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**3
      GO TO 165
  130 AMPJ=(FOUR*EIGHT/DSQRT(SIX))*ETAJ**(SEVEN/TWO)*EKAPPA*EKAPPA/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**3
      GO TO 165
  140 AMPJ=(FOUR*EIGHT/DSQRT(FIVE))*ETAJ**(ANINE/TWO)*EKAPPA*
     X     (ETAJ*ETAJ-EKAPPA*EKAPPA)/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**4
      GO TO 165
  150 AMPJ=(FOUR*EIGHT/DSQRT(FIVE*ANINE))*ETAJ**(SEVEN/TWO)*EKAPPA**2*
     X     (FIVE*ETAJ*ETAJ-EKAPPA*EKAPPA)/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**4
      GO TO 165
  160 AMPJ=(FOUR*EIGHT/DSQRT(FIVE))*ETAJ**(ANINE/TWO)*EKAPPA**3/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**4
      GO TO 165
  161 AMPJ=(FOUR*EIGHT/DSQRT(SEVEN*TEN))*ETAJ**(ANINE/TWO)*EKAPPA*
     X     (FIVE*ETAJ**4-TEN*(ETAJ*EKAPPA)**2+EKAPPA**4)/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**5
      GO TO 165
  162 AMPJ=(EIGHT*EIGHT/DSQRT(SEVEN*TEN))*ETAJ**(ANINE/TWO)*EKAPPA*
     X     ETAJ*EKAPPA*(FIVE*ETAJ*ETAJ-THREE*EKAPPA*EKAPPA)/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**5
      GO TO 165
  163 AMPJ=(EIGHT*EIGHT/DSQRT(SEVEN*TEN))*ETAJ**(ANINE/TWO)*EKAPPA*
     X     EKAPPA*EKAPPA*(SEVEN*ETAJ*ETAJ-EKAPPA*EKAPPA)/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**5
      GO TO 165
  164 AMPJ=(EIGHT*EIGHT*EIGHT/DSQRT(SEVEN*TEN))*ETAJ**(ANINE/TWO)*
     X     ETAJ*EKAPPA**4/
     X     (ETAJ*ETAJ+EKAPPA*EKAPPA)**5
  165 CONTINUE
      DO 166 JML=1,LLJP1
      JMSTO=JMSTO+1
      IF (JML.NE.IML) GO TO 166
      AMP=AMP+EVEC(IMSTO,KORB)*EVEC(JMSTO,KORB)*AMPI*AMPJ
  166 CONTINUE
  167 CONTINUE
  170 CONTINUE
  180 CONTINUE
      AMP=SCALE*EKAPPA*AMP/ESOURC
  185 ENERGY=-EVPERH*EVAL(KORB)
      WRITE(IW,2000) ENERGY,AMP
 2000 FORMAT(' ',F9.3,F13.3)
  190 CONTINUE
      RETURN
      END