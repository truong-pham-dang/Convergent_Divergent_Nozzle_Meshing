cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...File cdnozg.f                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      PROGRAM CDNOZG

c        Computes the grid for a converging-diverging nozzle (CD Nozzle)
c        as described in AIAA-87-0355 by M.S. Liou.  Geometry and grid 
c        are in units of inches.
c
c-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NIM, NJM, NKM

      PARAMETER ( NIM  = 501 ) 
      PARAMETER ( NJM  = 151 ) 
      PARAMETER ( NKM  =  31 ) 

      INTEGER IDIM
      INTEGER I, J, K
      INTEGER NI, NJ, NK

      REAL FAREA
      REAL PI
      REAL XX, AX, RX 
      REAL SPWALL
      REAL THETNZ, THET
      REAL SPRMI, SPRMAX

      REAL X (NIM,NJM,NKM)
      REAL Y (NIM,NJM,NKM)
      REAL Z (NIM,NJM,NKM)
      REAL S (NJM)

c-----------------------------------------------------------------------

c...Header. 

      WRITE(*,*) ' '
      WRITE(*,*) '--- cdnozg ---' 

c...Constants.

      PI = 4.0 * ATAN(1.0)

c...Specify dimension of domain.

      WRITE(*,*) ' '
      WRITE(*,*) 'Enter domain dimensions'
      WRITE(*,*) '(1) 2D, (2) Axisymmetric, (3) Quasi-3D, (4) 3D'
      READ (*,*) IDIM 

c...Specify the dimensions of the grid.  For 3D (idim=4), a 5 degree 
c...section of the nozzle is swept out.

      WRITE(*,*) ' '

      IF ( IDIM .EQ. 4 ) THEN
        WRITE(*,*) 'Enter ni, nj, nk'
        READ (*,*) NI, NJ, NK
        THETNZ = 5.0 * PI / 180.0 
      ELSE
        NK = 1
        WRITE(*,*) 'Enter ni, nj'
        READ (*,*) NI, NJ
      ENDIF

c...Specify the wall spacing.

      WRITE(*,*) ' '
      WRITE(*,*) 'Enter the wall spacing (inches)'
      READ (*,*) SPWALL

c...Generate the grid. 

      OPEN ( UNIT=7, FILE='rx.d', FORM='FORMATTED', STATUS='UNKNOWN' )

      SPRMAX = 0.0

      DO  I = 1, NI

        XX = 10.0 * ( I - 1.0 ) / ( NI - 1.0 )
        AX = FAREA( XX )
        RX = SQRT( AX / PI )
        WRITE(7,'(2(2X,F12.5))') XX, RX

c...    Compute the distribution of radial grid points.

        S(1)  = 0.0
        IF     ( IDIM .EQ. 1 ) THEN
          S(NJ) = AX
          CALL SPTANH ( NJ, S, SPWALL, SPWALL, 0, 0 )
        ELSEIF ( IDIM .EQ. 2 ) THEN
          S(NJ) = RX
          CALL SPTANH ( NJ, S, SPWALL, SPWALL, 0, 0 )
        ELSEIF ( IDIM .EQ. 3 ) THEN
          S(NJ) = 1.0
          CALL SPEVEN ( NJ, S )
        ELSEIF ( IDIM .EQ. 4 ) THEN
          S(NJ) = RX
          CALL SPTANH ( NJ, S, SPWALL, SPWALL, 0, 0 )
        ENDIF
    
        CALL GSPRAT ( S, NJ, SPRMI )
        IF ( SPRMI .GT. SPRMAX )  SPRMAX = SPRMI

c...    Load grid.

        IF     ( IDIM .EQ. 1  .OR.  IDIM .EQ. 2 ) THEN
          DO  J = 1, NJ
            X(I,J,1) = XX
            Y(I,J,1) = S(J)
            Z(I,J,1) = 1.0
          ENDDO
        ELSEIF ( IDIM .EQ. 3 ) THEN
          DO  J = 1, NJ
            X(I,J,1) = XX
            Y(I,J,1) = S(J)
            Z(I,J,1) = RX
          ENDDO
        ELSEIF ( IDIM .EQ. 4 ) THEN
          DO  J = 1, NJ
            DO  K = 1, NK
              THET     = THETNZ * ( K - 1.0 ) / ( NK - 1.0 )
              X(I,J,K) = XX
              Y(I,J,K) = S(J) * COS(THET)
              Z(I,J,K) = S(J) * SIN(THET)
            ENDDO
          ENDDO
        ENDIF

      enddo 

c...Write out grid quality statements.

      WRITE (*,*) ' '
      WRITE (*,*) 'sprmax = ', SPRMAX

c...Write out the grid as a Plot3d file.

      WRITE (*,*) ' '
      WRITE (*,*) 'Writing out grid (cdnoz.x)'

      OPEN ( unit=8, file='cdnoz.x', form='unformatted', 
     &       status='unknown' )

      WRITE (8) NI, NJ, NK
      WRITE (8) ((( X(I,J,K), I=1,NI), J=1,NJ), K=1,NK),
     &          ((( Y(I,J,K), I=1,NI), J=1,NJ), K=1,NK),
     &          ((( Z(I,J,K), I=1,NI), J=1,NJ), K=1,NK)

c...Write out the grid as a mesh file which is read by project Euler2D_Nozzle.
	IF     ( IDIM .EQ. 1  .OR.  IDIM .EQ. 2 ) THEN
		OPEN ( UNIT=9, FILE='cdnoz.mesh', STATUS='REPLACE' )
		WRITE (9,*) NJ, NI
		DO J = 1, NJ
			DO I = 1, NI
				WRITE(9,*) X(I,J,1), Y(I,J,1)
			ENDDO
		ENDDO
	ENDIF
c...End Statement.

      WRITE (*,*) '     '
      WRITE (*,*) '--- End of cdnozg ---'
      WRITE (*,*) '     '

c-----------------------------------------------------------------------

      STOP
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION FAREA ( X )

c...Function for area of converging-divirging duct from aiaa-87-0355, 
c...M.S. liou (reference 75).
c
c-----------------------------------------------------------------------

      IMPLICIT NONE
      REAL FAREA, X, PI

c-----------------------------------------------------------------------

      PI = 4.0 * ATAN(1.0)

      IF ( X .LT. 5.0 ) THEN
        FAREA = 1.75 - 0.75 * COS( ( 0.2 * X - 1.0 ) * PI )
      ELSE
        FAREA = 1.25 - 0.25 * COS( ( 0.2 * X - 1.0 ) * PI )
      ENDIF

c-----------------------------------------------------------------------

      RETURN
      END



