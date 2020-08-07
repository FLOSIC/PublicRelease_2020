C UTEP Electronic Structure Lab (2020)
 
        SUBROUTINE BRTFRC(N,M,A,B,ALP,BET,RK,SSR,SSI)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION F(0:4),G(0:4)
              C=(ALP*A+BET*B)/(ALP+BET)
              AA=C-A
              BB=C-B
              ARG=ALP*BET*(A-B)*(A-B)/(ALP+BET)
C             ARG=ARG+RK*RK/(ALP+BET)/4.0D0
                   EXPN=EXPNL(-ARG)
              ARGC=EXPN*COSF(C*RK)
              ARGS=EXPN*SINF(C*RK)
                CALL GFT(4,4,ALP+BET,RK,F,G)
C F(N)= INT X**(2*N+1)*SINF(KX)*EXPNL(-ALPHA*X*X)DX (-INF TO INF)
C G(N)= INT X**(2*N )*COSF(KX)*EXPNL(-ALPHA*X*X)DX (-INF TO INF)
                SR=0.0D0
                SI=0.0D0
              DO J=0,N
              DO L=0,M
                COEF=                                 
     &                (FCTRL(N)/FCTRL(N-J)/FCTRL(J))*
     &                (FCTRL(M)/FCTRL(M-L)/FCTRL(L))*
     &                (AA**J)*(BB**L)
                NT=N+M-L-J
                IF((NT/2)*2.EQ.NT)THEN
                   SR=SR+COEF*G(NT/2)
                ELSE
                   SI=SI+COEF*F(NT/2)
                END IF
              END DO
              END DO
              SSR=SR*ARGC-SI*ARGS
              SSI=SI*ARGC+SR*ARGS
          RETURN
          END
