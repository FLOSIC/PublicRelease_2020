C UTEP Electronic Structure Lab (2020)
       FUNCTION DISTANCE(A,B)
C
C      RETURNS THE DISTANCE BETWEN VECTORS A AND B
C
C      BY ULISES REVELES, JULY 2013.
C
C      -----------------------------------------------------------------
C
       IMPLICIT NONE
C
       REAL*8 A(*),B(*),DISTANCE
C
C      -----------------------------------------------------------------
C
       DISTANCE = 0.0
       DISTANCE = (A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2
       DISTANCE = SQRT(DISTANCE)
C
       RETURN
C
C      -----------------------------------------------------------------
C
       END
