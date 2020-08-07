C UTEP Electronic Structure Lab (2020)
c
c ************************************************************
c
c swap - swap two real values
c
       SUBROUTINE SWAP(A,B)
        IMPLICIT REAL*8 (A-H,O-Z)
        SAVE
        X= A
        A= B
        B= X
        RETURN
       END
