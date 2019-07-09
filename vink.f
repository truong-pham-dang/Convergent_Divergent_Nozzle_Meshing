cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
c.....VINK
c.....Distributes points on a curve so that specified derivatives
c.....at the edges are satisfied (from NASA CR 3313 by Vinokur, 1980).
c.....Vinokur's algorithm is designed to distribute a specified 
c.....number of points along a curve, given the number of points, 
c.....the length of the curve, and the deriviative conditions at 
c.....both ends of the curve. In CFD applications, the user would
c.....usually rather specify the delta-s's at the ends of the curve, 
c.....which are equivalent to the derivatives only to first order. 
c.....Therefore, the user may wish to apply this algorithm iteratively
c.....to obtain an exact delta-s specification. Subroutine 
c.....vinkiter will iterate on this scheme until the proper delta-s 
c.....constraints are satisfied. 
c.....Inputs:
c.....   lmax -> number of points on the curve
c.....   smin, smax -> beginning and end values of s
c.....   ds1,  ds2  -> the derivative end conditions input into 
c.....                 Vinokur's function
c.....   kase = 0 -> satisfy delta-s on both ends
c.....        = 1 -> satisfy delta-s only at xi=ximax     
c.....        = 2 -> satisfy delta-s only at xi=ximin
c.....Outputs:
c.....   s(xi) -> resulting s distribution from Vinokur's function
c.....   es1 -> ( s(ximin+1)-s(ximin) ) <-  calculated delta-s
c.....   es2 -> ( s(ximax)-s(ximax-1) ) <-
c.....Note:
c.....   Additionally, this version uses the approximate inverse 
c.....   solution for y=sin(x)/x and y=sinh(x)/x rather than a Newton 
c.....   iteration. The approximate solution was also taken from 
c.....   NASA CR 3313.

      SUBROUTINE VINK( S, LMAX, SMIN, SMAX, DS1, DS2, ES1, ES2, KASE )

c#ifdef DP
c      implicit double precision (a-h,o-z)
c#endif

      DIMENSION S(*), D1(4,2), D2(4,2)

      PI = 4.0 * ATAN(1.0)

c.....Initialization.

      SDEL = SMAX - SMIN
      S0 = SDEL / REAL( LMAX-1 ) / DS1
      S1 = SDEL / REAL( LMAX-1 ) / DS2
      B = SQRT( S0 * S1 )
      A = SQRT( S0 / S1 )
      IF ( KASE .EQ. 1 ) THEN
         B = S1
      ELSE IF ( KASE .EQ. 2 ) THEN
         B = S0
      ENDIF

c.....Calculate x based on value of B.

      IF ( B .LT. 1.0 ) THEN

c........x is real.

         IF ( B .LT. 0.26938972 ) THEN
            X = ( 1.0 
     &          - B 
     &          + B**2 
     &          - ( 1.0 + PI ** 2.0 / 6.0 ) * B**3
     &          +  6.794732 * B**4 
     &          - 13.205501 * B**5  
     &          + 11.726095 * B**6 ) 
     &        * PI
         ELSE
            C = 1.0 - B
            X = ( 1.0
     &          + 0.15 * C 
     &          + 0.057321429 * C**2 
     &          + 0.048774238 * C**3
     &          - 0.053337753 * C**4 
     &          + 0.075845134 * C**5 )
     &        * SQRT( 6.0 * C )
         ENDIF

      ELSE IF ( B .EQ. 1.0 ) THEN

c........x is zero.

         X = 0.0

      ELSE

c........x is imaginary.

         IF ( B .LT. 2.7829681 ) THEN
            C = B - 1.0
            X = ( 1.0
     &          - 0.15 * C 
     &          + 0.0573214290 * C**2 
     &          - 0.0249072950 * C**3
     &          + 0.0077424461 * C**4 
     &          - 0.0010794123 * C**5 )
     &        * SQRT( 6.0 * C )
         ELSE
            V = LOG(B)
            W = 1.0 / B - 0.028527431
            X = V + ( 1.0 + 1.0 / V ) * LOG( 2.0*V ) - 0.02041793
     &        + 0.24902722 * W 
     &        + 1.9496443 * W**2
     &        - 2.6294547 * W**3 
     &        + 8.56795911 * W**4
         ENDIF

      ENDIF

c.....Distribute points along edge.

      IF ( KASE .EQ. 1 .OR. KASE .EQ. 2 ) THEN
         S(1   ) = 0.0
         S(LMAX) = SDEL
         DO 9 I = 2, LMAX-1
            J = LMAX + 1 - I
            XI = REAL(I-1) / (LMAX-1)
            IF ( B .GT. 1.0001 ) THEN
               U1 = 1.0 + TANH( X/2.0 * ( XI - 1.0 ) )
     &                  / TANH( X/2.0 )
            ELSE IF ( B .LT. 0.9999 ) THEN
               U1 = 1.0 + TAN ( X/2.0 * ( XI - 1.0 ) )
     &                  / TAN ( X/2.0 )
            ELSE
               U1 = XI * ( 1.0 - 0.5 * ( B - 1.0 ) * ( 1.0 - XI ) 
     &            * ( 2.0 - XI ) )
            ENDIF
            U2 = SINH( XI*X ) / SINH( X )
            IF ( KASE .EQ. 1 ) THEN
               FACT = ABS( DS1 )
               S(J) = ( ( 1.0 - FACT ) * ( 1.0 - U1 ) 
     &              + FACT * ( 1.0 - U2 ) ) * SDEL
            ELSE IF ( KASE .EQ. 2 ) THEN
               FACT = ABS( DS2 )
               S(I) = ( ( 1.0 - FACT ) * U1 + FACT * U2 ) * SDEL
            ENDIF
    9    CONTINUE
      ELSE
         DO 5 I = 1, LMAX
            XI = REAL(I-1) / REAL(LMAX-1)
            CNUM = X * ( XI-0.5 )
            CDEN = X / 2.0
            IF ( B .LT. 0.9999 ) THEN
               CC = TAN(CNUM) / TAN(CDEN)
               U = 0.5 * ( 1.0 + CC )
            ELSE IF ( B .GE. 0.9999 .AND. B .LE. 1.0001 ) THEN
               U = XI * ( 1.0 + 2.0 * ( B - 1.0 ) * ( XI - 0.5 ) 
     &           * ( 1.0 - XI ) )
            ELSE IF ( B .GT. 1.0001 ) THEN
               CC = TANH(CNUM) / TANH(CDEN)
               U = 0.5 * ( 1.0 + CC )
            ENDIF
            S(I) = U * SDEL / ( A + ( 1.0 - A ) * U )
    5    CONTINUE
      ENDIF

      DO 8 L = 1, LMAX
         S(L) = S(L) + SMIN
    8 CONTINUE
      ES1 = S(   2) - S(     1)
      ES2 = S(LMAX) - S(LMAX-1)

      RETURN
      END

