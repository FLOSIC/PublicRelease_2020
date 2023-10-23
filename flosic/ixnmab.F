C UTEP Electronic Structure Lab (2020)
C
        SUBROUTINE IXNMAB(RK,ALP,BET,AA,BB,NM,MM,CT,ST)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDMN=4)
        LOGICAL FIRST
        DIMENSION F(0:NDMN),G(0:NDMN),BIN(0:NDMN,0:NDMN)
        DIMENSION CT(0:NDMN,0:NDMN),ST(0:NDMN,0:NDMN)
        DATA FIRST/.TRUE./
           IF(FIRST)THEN
           FIRST=.FALSE.
C CALCULATE BINOMIAL COEFFICIENTS:
               DO NTOP=0,NDMN
               DO NBOT=0,NTOP
                 TOP=FCTRL(NTOP)
                 BOT=FCTRL(NBOT)*FCTRL(NTOP-NBOT)
                 BIN(NTOP,NBOT)=TOP/BOT
               END DO
               END DO
           END IF
            IF(NM+MM.GT.NDMN)STOP'NDMN'
            GAM=ALP+BET
            C=(ALP*AA+BET*BB)/GAM
            A=C-AA
            B=C-BB
              ARG=ALP*BET*(A-B)*(A-B)/GAM
C    &           +RK*RK/4.0D0/GAM
            ARG=EXPNL(-ARG)
C           ARG=EXPNL(-ALP*BET*(A-B)*(A-B)/GAM,2)
C           ARG=ARG*EXPNL(I*RK*C) (FIXED DOWN BELOW)
        CALL GFT(NDMN,NM+MM,GAM,RK,F,G)
              DO I=0,NM+MM
                G(I)=G(I)*ARG
                F(I)=F(I)*ARG
              END DO 
        DO N=0,NM
        DO M=0,MM
          CT(N,M)=0.0D0
          ST(N,M)=0.0D0
        END DO
        END DO
        DO N=0,NM
         AP=1.0D0
         DO J=0,N 
         IF(J.NE.0)AP=AP*A
        DO M=0,MM
         BP=AP*BIN(N,J)
         DO L=0,M
         IF(L.NE.0)BP=BP*B
           COEF=BIN(N,J)*BIN(M,L)*(B**L)*(A**J)
C          COEF=BIN(M,L)*BP
           NT=N+M-J-L
            IF(2*(NT/2).EQ.NT)THEN
            CT(N,M)=CT(N,M)+G(NT/2)*COEF
            ELSE
            ST(N,M)=ST(N,M)+F(NT/2)*COEF
            END IF
          END DO
          END DO
        END DO
        END DO
                COSD=COSF(RK*C)
                SIND=SINF(RK*C)
            DO N=0,NM
            DO M=0,MM
             PR=CT(N,M)*COSD-ST(N,M)*SIND 
             PI=CT(N,M)*SIND+ST(N,M)*COSD
             CT(N,M)=PR
             ST(N,M)=PI
            END DO
            END DO
        RETURN
        END
