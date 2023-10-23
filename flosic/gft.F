C UTEP Electronic Structure Lab (2020)
C
        SUBROUTINE GFT(NDMN,NMAX,ALPHA,RK,F,G)
C F(N)= INT X**(2*N+1)*SINF(KX)*EXPNL(-ALPHA*X*X)DX (-INF TO INF)
C G(N)= INT X**(2*N  )*COSF(KX)*EXPNL(-ALPHA*X*X)DX (-INF TO INF)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION F(0:NDMN),G(0:NDMN)
        DATA PI/3.14159265358979312D0/
             IF(NDMN.LT.NMAX)THEN
             write(6,*)'PARAMETER ERROR IN GFT'
             STOP
             END IF
            F3=1.0D0/ALPHA
            F2=F3/2 
            F1=RK*F2
        G(0)=SQRT(PI/ALPHA)*EXPNL(-RK*F1/2)
        F(0)=F1*G(0)
C         write(6,*)ALPHA,RK
C         write(6,*)G(0),F(0)
          DO I=1,NMAX
          G(I)=(2*I-1)*F2*G(I-1)-F1*  F(I-1)
          F(I)=F1        *G(I-1)+F3*I*F(I-1)
C         write(6,*)G(I),F(I)
          END DO
        RETURN 
        END
