cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   The subroutine sptanh is a front end for the routine vinkiter and  c
c   vink                                                               c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SPTANH ( LMAX, S, DSAE, DSBE, KASE, NLAST )

c      implicit double precision (a-h,o-z)

      DIMENSION S(LMAX)
      INTEGER MXIT

      SMIN = S(1)
      SMAX = S(LMAX)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c.....VINKITER
c.....Vinokur's function creates a distribution of grid points which 
c.....satisfy a specified derivative condition, but we require a delta-s
c.....constraint instead. These two values are equivalent only to first 
c.....order, and hence, we resort to an iterative procedure to obtain  
c.....more accurate delta-s's. Up to ten iterative sweeps are
c.....made. The first guess sets ds/dxi = delta-s. The next guess 
c.....recalculates ds/dxi using the leading term in the truncation error
c.....(d2s/d(xi)**2). The next eight iterations use a 2-d secant 
c.....algorithm to home in on the ds/dxi's at both ends which will give 
c.....the correct delta-s.  In the cases where a single-sided 
c.....stretching function is required, (kase = 1 or 2) a secant 
c.....algorithm in 1-d is applied instead. 
c.....Inputs:
c.....   kase = 0 -> delta-s specified on both ends so use a 2D secant 
c.....               method to arrive at the values of dsa and dsb 
c.....               that will satisfy ds1e and ds2e within roundoff.
c.....   kase = 1 -> delta-s specified at the last endpoint only, so 
c.....               use a 1-D secant method to arrive at the value 
c.....               of dsb that will satisfy dsbe within roundoff.
c.....   kase = 2 -> delta-s specified at the first endpoint only, so 
c.....               use a 1-D secant method to arrive at the value 
c.....               of dsa which will satisfy ds1e within roundoff.
c
c      subroutine vinkiter( s, lmax, smin, smax, dsae, dsbe, 
c     &                     kase, nlast )
c
c#ifdef DP
c      implicit double precision (a-h,o-z)
c#endif

c      dimension s(*)
c      integer mxit

c#ifdef DP
c      tolmin2 = 1.0d-14
c#else
      TOLMIN2 = 1.0E-07
c#endif
      TOLMIN  = TOLMIN2 * TOLMIN2

      MXIT = 20

      IF ( KASE .EQ. 0 ) THEN
  
c........Initial guess, an = dsae, bn = dsbe.

         AN2 = DSAE 
         BN2 = DSBE  
         CALL VINK(S,LMAX,SMIN,SMAX,AN2,BN2,ESA,ESB,KASE )
         FN2 = ESA / DSAE - 1
         GN2 = ESB / DSBE - 1

c........Second guess, calculate ds1 and ds2 from a truncated 
c........Taylor Series.

         DSSA = 2.*S(   1)-5.*S(     2)+4.*S(     3)   -S(     4)
         DSSB = 2.*S(LMAX)-5.*S(LMAX-1)+4.*S(LMAX-2)   -S(LMAX-3)
         AN1 = DSAE-0.5*DSSA
         BN1 = DSBE+0.5*DSSB
         CALL VINK(S,LMAX,SMIN,SMAX,AN1,BN1,ESA,ESB,KASE )
         FN1 = ESA/DSAE - 1
         GN1 = ESB/DSBE - 1
         AN  = AN1
         BN  = BN1

c........3rd thru nth guesses, use 2-d secant method.

         DO 10 N = 3, MXIT 
c...........Calculate offset derivatives.
            CALL VINK(S,LMAX,SMIN,SMAX,AN2,BN1,ESA21,ESB21,KASE )
            CALL VINK(S,LMAX,SMIN,SMAX,AN1,BN2,ESA12,ESB12,KASE )
            FA = ( ESA - ESA21 )/DSAE
            FB = ( ESA - ESA12 )/DSAE
            GA = ( ESB - ESB21 )/DSBE
            GB = ( ESB - ESB12 )/DSBE
            DEN = FA*GB - FB*GA
            DELA = -(AN1 - AN2)
            DELB = -(BN1 - BN2)
c...........Stick with last guess if approaching roundoff.
            IF ( ABS(DEN) .LT. TOLMIN)THEN
               CALL VINK(S,LMAX,SMIN,SMAX,AN,BN,ESA,ESB,KASE )
               RETURN 
            END IF
c...........Calculate next distribution.
            AN = AN1 + DELA*(  GB*FN1 - FB*GN1 )/DEN
            BN = BN1 + DELB*(- GA*FN1 + FA*GN1 )/DEN
            CALL VINK(S,LMAX,SMIN,SMAX,AN,BN,ESA,ESB,KASE )
            FN = ESA/DSAE - 1
            GN = ESB/DSBE - 1
c...........Update n, n-1, n-2 and continue.
            AN2 = AN1
            BN2 = BN1
            AN1 = AN
            BN1 = BN
            FN1 = FN
            GN1 = GN
   10    CONTINUE

      ELSE IF ( KASE .EQ. 1 ) THEN
  
c........Initial guess, bn = dsbe.

         BN2 = DSBE 
         CALL VINK(S,LMAX,SMIN,SMAX,DSAE,BN2,ESA,ESB,KASE )
         FN2 = ESB/DSBE - 1

c........Second guess, calculate ds1 and ds2 from a truncated 
c........Taylor Series.

         DSSB = 2.*S(LMAX)-5.*S(LMAX-1)+4.*S(LMAX-2)   -S(LMAX-3)
         BN1 = DSBE-0.5*DSSB
         CALL VINK(S,LMAX,SMIN,SMAX,DSAE,BN1,ESA,ESB,KASE )
         FN1 = ESB/DSBE - 1
         BN = BN1

c........3rd thru nth guesses, use 1-d secant method.

         DO 20 N = 3, MXIT 
c...........Stick with last guess if approaching roundoff.
            DEN = FN1-FN2
            IF ( ABS(DEN) .LT. TOLMIN2)THEN
               CALL VINK(S,LMAX,SMIN,SMAX,DSAE,BN,ESA,ESB,KASE )
               RETURN
            END IF
c...........Calculate next distribution.
            BN = BN1 - FN1*(BN1-BN2)/DEN 
            CALL VINK(S,LMAX,SMIN,SMAX,DSAE,BN,ESA,ESB,KASE )
            FN = ESB/DSBE - 1
c...........Update n, n-1, n-2 and continue.
            BN2 = BN1
            BN1 = BN
            FN2 = FN1
            FN1 = FN
   20    CONTINUE

      ELSE IF ( KASE .EQ. 2 ) THEN
  
c........Initial guess, an = dsae.

         AN2 = DSAE 
         CALL VINK(S,LMAX,SMIN,SMAX,AN2,DSBE,ESA,ESB,KASE )
         FN2 = ESA/DSAE - 1

c........Second guess, calculate ds1 and ds2 from a truncated 
c........Taylor Series.

         DSSA = 2.*S(   1)-5.*S(     2)+4.*S(     3)   -S(     4)
         AN1 = DSAE-0.5*DSSA
         CALL VINK(S,LMAX,SMIN,SMAX,AN1,DSBE,ESA,ESB,KASE )
         FN1 = ESA/DSAE - 1
         AN = AN1

c........3rd thru nth guesses, use 1-d secant method.

         DO 30 N = 3, MXIT 
c...........Stick with last guess if approaching roundoff.
            DEN = FN1-FN2
            IF(ABS(DEN) .LT. TOLMIN2)THEN
               CALL VINK(S,LMAX,SMIN,SMAX,AN,DSBE,ESA,ESB,KASE )
               RETURN
            END IF
c...........Calculate next distribution.
            AN = AN1 - FN1*(AN1-AN2)/DEN 
            CALL VINK(S,LMAX,SMIN,SMAX,AN,DSBE,ESA,ESB,KASE )
            FN = ESA/DSAE - 1
c...........Update n, n-1, n-2 and continue.
            AN2 = AN1
            AN1 = AN
            FN2 = FN1
            FN1 = FN
   30    CONTINUE

      END IF

      RETURN 
      END 