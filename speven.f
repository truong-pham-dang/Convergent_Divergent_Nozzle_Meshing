cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...file  speven.f  ( evenly spaced grid )                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SPEVEN ( N, S )

c     Distributes s evenly given number of points, n, and the boundary 
c     values, s(1) and s(n).
c
c     n       Number of points
c     s       Parametric variable
c
c-----------------------------------------------------------------------

      IMPLICIT NONE

c.....Argument variables.

      INTEGER N

      REAL S (N)

c.....Local variables.

      INTEGER I

c-----------------------------------------------------------------------

c.....Perform a uniform distribution.

      DO  I = 2, N-1
        S(I) = S(1) + ( I - 1. ) * ( S(N) - S(1) ) / ( N - 1. )
      ENDDO

c-----------------------------------------------------------------------

      RETURN
      END