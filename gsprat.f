cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...file  gsprat.f  ( grid spacing ratio )                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GSPRAT ( S, N, SPRMAX )

c     Computes the maximum ratio between two grid points for a curve.
c   
c-----------------------------------------------------------------------

      INTEGER I, N
      REAL S(*), SPA, SPB, SPR, SPRMAX

c-----------------------------------------------------------------------

      SPA    = S(2) - S(1)
      SPB    = S(3) - S(2)
      SPRMAX = MAX( SPA, SPB ) / MIN( SPA, SPB ) - 1.0

      DO  I = 3, N - 1
        SPA    = S(I) - S(I-1)
        SPB    = S(I+1) - S(I)
        SPR    = MAX( SPA, SPB ) / MIN( SPA, SPB ) - 1.0
        SPRMAX = MAX( SPRMAX, SPR )
      ENDDO

c-----------------------------------------------------------------------

      RETURN
      END