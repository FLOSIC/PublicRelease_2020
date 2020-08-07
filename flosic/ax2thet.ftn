C UTEP Electronic Structure Lab (2020)
C
       SUBROUTINE AX2THET(AX,THETA,BETA)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION AX(3)
       ANORM=0.0D0
       DO IX=1,3
        ANORM=ANORM+AX(IX)**2 
       ENDDO
       ANORM=SQRT(ANORM)
       IF(ANORM.LT.1D-6) THEN
        THETA=0.0D0
        BETA=0.0D0
        RETURN
       ENDIF
       DO IX=1,3
        AX(IX)=AX(IX)/ANORM
       ENDDO 
       THETA=ACOS(AX(3))
       ST=SIN(THETA)
       IF(ABS(ST).LT.1D-6) THEN
        BETA=0.0D0
        RETURN
       ENDIF
       BETA=ACOS(AX(1)/ST)
       RETURN
       END 
