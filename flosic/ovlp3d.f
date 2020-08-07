C UTEP Electronic Structure Lab (2020)
         SUBROUTINE OVLP3D(NDIM,ALPI,ALPJ,AI,AJ,SS)
c m.r. pederson 21-june 2000
c replacement for ovmxsf 
c note that dimension of ss must be passed...
c ndim=10  s,p and d states.
c ndim=20  s,p,d and f states....
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION AI(3),AJ(3),SS(NDIM,NDIM)
         DIMENSION OV(0:3,0:3,3)
         DIMENSION N(3,20)
         DATA N/0,0,0,  1,0,0, 0,1,0,  0,0,1,  
     &                  2,0,0, 0,2,0,  0,0,2,  
     &                  1,1,0, 1,0,1,  0,1,1,
     &                  3,0,0, 0,3,0,  0,0,3, 
     &                  2,1,0, 2,0,1,  0,2,1,
     &                  1,2,0, 1,1,1,  1,0,2,
     &                  0,1,2 /
           IF(NDIM.NE.10.AND.NDIM.NE.20)THEN
           write(6,*)'NDIM:',NDIM,' MUST BE 10 OR 20...'
           write(6,*)'INCORRECT USAGE OF OVLP3D'
           STOP
           END IF
                 DO IX=1,3
                 CALL OVLP1D(ALPI,ALPJ,AI(IX),AJ(IX),3,OV(0,0,IX))
                 END DO 
             DO I=1,NDIM
             DO J=1,NDIM
             SS(I,J)=OV(N(1,I),N(1,J),1)*
     &               OV(N(2,I),N(2,J),2)*
     &               OV(N(3,I),N(3,J),3) 
             END DO
             END DO
         RETURN
         END
