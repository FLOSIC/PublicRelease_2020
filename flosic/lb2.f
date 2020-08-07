C UTEP Electronic Structure Lab (2020)
C     ******
C
C
C   ----------------------------------------------------------
C     DATA 
C   ----------------------------------------------------------
C
      BLOCK DATA LB2
      INTEGER LP,MP
      DOUBLE PRECISION GLTOL,STPMIN,STPMAX
      COMMON /LB3/MP,LP,GLTOL,STPMIN,STPMAX
      DATA MP,LP,GLTOL,STPMIN,STPMAX/6,6,9.0D-01,1.0D-20,1.0D+20/
      END
