C UTEP Electronic Structure Lab (2020)

        SUBROUTINE DOINT(II,JJ,ALP,BET,A,B,RK,DR,DI)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION N(3),M(3),A(3),B(3),SSR(3),SSI(3)
        DIMENSION RK(3)
        DIMENSION INDEX(3,10)
         DO I=1,10
            INDEX(1,I)=0
            INDEX(2,I)=0
            INDEX(3,I)=0
         END DO
            DO J=1,3
            INDEX(J,J+1)=1
            INDEX(J,J+4)=2
            END DO
            INDEX(1,8)=1
            INDEX(2,8)=1
            INDEX(1,9)=1
            INDEX(3,9)=1
            INDEX(2,10)=1
            INDEX(3,10)=1
            DO I=1,3
            N(I)=INDEX(I,II)
            M(I)=INDEX(I,JJ)
            END DO
            DO I=1,3
            CALL BRTFRC(N(I),M(I),A(I),B(I),ALP,BET,RK(I),SSR(I),SSI(I))
            END DO
            DR= SSR(1)*SSR(2)*SSR(3)-
     &          SSI(1)*SSI(2)*SSR(3)-
     &          SSI(1)*SSR(2)*SSI(3)-
     &          SSR(1)*SSI(2)*SSI(3)
            DI=-SSI(1)*SSI(2)*SSI(3)+
     &          SSI(1)*SSR(2)*SSR(3)+
     &          SSR(1)*SSI(2)*SSR(3)+
     &          SSR(1)*SSR(2)*SSI(3)
        RETURN
        END 
