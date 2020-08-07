C UTEP Electronic Structure Lab (2020)
         SUBROUTINE KINM3D(NDIM,ALPI,ALPJ,AI,AJ,SS)
c m.r. pederson 21-june 2000
         IMPLICIT REAL*8 (A-H,O-Z)
c         INCLUDE 'const.inc'
         LOGICAL FIRST,EXIST
         PARAMETER (NFMAX=3)
         PARAMETER (NRELMAX=6)
         PARAMETER (NDMAX=2*NRELMAX+2)
         LOGICAL REL(NRELMAX)
         DIMENSION T1(0:NFMAX,0:NFMAX,0:NDMAX,3)
         DIMENSION RL(20,20)
         DATA FIRST/.TRUE./
       DATA HA2EV/27.2116D0/
       DATA CLIGHT/137.088D0/
       DATA HA2KEL/315891.1826D0/        !27.2116*1.602/1.38E-04   
       DATA PI /3.141592654D0/
       DATA TEMP/1.0D-4/
c replacement for knmxsf 
c note that dimension of ss must be passed...
c ndim=10  s,p and d states.
c ndim=20  s,p,d and f states....
c rel1 = .true. include del^4 scalar relativistic term...
c rel2 = .true. include del^6 scalar relativistic term....
c rel3 = .true. include del^8 scalar relativistic term....
c rel4 = .true. include del^10 scalar relativistic term....
         DIMENSION AI(3),AJ(3),SS(NDIM,NDIM)
         DIMENSION N(3,20)
         DATA N/0,0,0,  1,0,0, 0,1,0,  0,0,1,  
     &                  2,0,0, 0,2,0,  0,0,2,  
     &                  1,1,0, 1,0,1,  0,1,1,
     &                  3,0,0, 0,3,0,  0,0,3, 
     &                  2,1,0, 2,0,1,  0,2,1,
     &                  1,2,0, 1,1,1,  1,0,2,
     &                  0,1,2 /
           SAVE
           NRELA=0
           IF(FIRST) THEN
           DO I=1,NRELMAX
            REL(I)=.FALSE.
           ENDDO
           INQUIRE(FILE='RELA',EXIST=EXIST)
           IF (EXIST) THEN
             OPEN(71,FILE='RELA',STATUS='OLD')
             READ(71,*) NRELA
             write(6,*)'NRELA : ',NRELA
             IF(NRELA.GT.NRELMAX) NRELA=NRELMAX
             CLOSE(71)
             FIRST=.FALSE.
           ENDIF 
           DO I=1,NRELA
            REL(I)=.TRUE.
           ENDDO
           ENDIF 
C
           IF(NDIM.NE.10.AND.NDIM.NE.20)THEN
             write(6,*)'NDIM:',NDIM,' MUST BE 10 OR 20...'
             write(6,*)'INCORRECT USAGE OF KINM3D'
             CALL STOPIT
           END IF
             MXD=2
             IF(REL(1))MXD=4
             IF(REL(2))MXD=6
             IF(REL(3))MXD=8
             IF(REL(4))MXD=10
             IF(REL(5))MXD=12
             IF(REL(6))MXD=14
             DO IX=1,3
             CALL FONEDM(ALPI,ALPJ,AI(IX),AJ(IX),3,MXD,T1(0,0,0,IX))
             END DO 
c second order term
             DO I=1,NDIM
             DO J=1,NDIM
             SS(I,J)=T1(N(1,I),N(1,J),2,1)*
     &               T1(N(2,I),N(2,J),0,2)*
     &               T1(N(3,I),N(3,J),0,3) 
             SS(I,J)=T1(N(1,I),N(1,J),0,1)*
     &               T1(N(2,I),N(2,J),2,2)*
     &               T1(N(3,I),N(3,J),0,3) +SS(I,J)
             SS(I,J)=T1(N(1,I),N(1,J),0,1)*
     &               T1(N(2,I),N(2,J),0,2)*
     &               T1(N(3,I),N(3,J),2,3) +SS(I,J)
             END DO
             END DO
             JREL=1
             DO I=1,NDIM
             DO J=1,NDIM
             IF(JREL.EQ.1) THEN
               SS(J,I)=-0.5D0*SS(J,I) 
             ELSE
               SS(I,J)=-0.5D0/(1.0D0+SS(I,J)/(2.0D0*CLIGHT)**2)*SS(I,J)
             ENDIF
             END DO
             END DO
c fourth order term:
             IF(REL(1))THEN
             CSQR=CLIGHT**2
              write(6,*)' P^4 PART : ', REL(1)
             DO I=1,NDIM
             DO J=1,NDIM
             RL(I,J)=      T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),0,3) 
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=2.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=2.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=2.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             SS(I,J)=SS(I,J)-RL(I,J)/(8.0*CSQR)
             END DO
             END DO
             END IF
c sixth order term:
             IF(REL(2))THEN
              write(6,*)' P^6 PART : ', REL(2)
             DO I=1,NDIM
             DO J=1,NDIM
             RL(I,J)=      T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),0,3) 
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=3.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=3.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=3.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=3.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=3.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=3.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             SS(I,J)=SS(I,J)-RL(I,J)/(16.0*CSQR*CSQR)
             END DO
             END DO
             END IF
c eigth order term:
             IF(REL(3))THEN
              write(6,*)' P^8 PART : ', REL(3)
             DO I=1,NDIM
             DO J=1,NDIM
             RL(I,J)=      T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),0,3) 
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)
             RL(I,J)=4.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=4.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=4.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=4.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=4.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=4.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)      
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)      
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)      
             RL(I,J)=12.D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=12.D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=12.D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             SS(I,J)=SS(I,J)-5.0D0*RL(I,J)/(128.D0*CSQR**3)
             END DO
             END DO
             END IF
c tenth order term:
             IF(REL(4))THEN
              write(6,*)' P^10 PART : ', REL(4)
             DO I=1,NDIM
             DO J=1,NDIM
             RL(I,J)=      T1(N(1,I),N(1,J),10,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),0,3) 
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),10,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),10,3) +RL(I,J)
             RL(I,J)=5.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=5.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=5.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=5.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=5.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)
             RL(I,J)=5.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)
             RL(I,J)=10.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)      
             RL(I,J)=10.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)      
             RL(I,J)=10.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)      
             RL(I,J)=10.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)      
             RL(I,J)=10.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)      
             RL(I,J)=10.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)      
             RL(I,J)=20.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=20.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=20.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=30.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=30.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=30.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             SS(I,J)=SS(I,J)-7.0D0*RL(I,J)/(256.0D0*CSQR**4)
             END DO
             END DO
             END IF
c 12'th order term:
             IF(REL(5))THEN
              write(6,*)' P^12 PART : ', REL(5)
             DO I=1,NDIM
             DO J=1,NDIM
             RL(I,J)=      T1(N(1,I),N(1,J),12,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),0,3) 
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),12,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),12,3) +RL(I,J)
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),10,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),10,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),10,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),10,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),10,3) +RL(I,J)
             RL(I,J)=6.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),10,3) +RL(I,J)
             RL(I,J)=15.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)      
             RL(I,J)=15.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)      
             RL(I,J)=15.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)      
             RL(I,J)=15.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)      
             RL(I,J)=15.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)      
             RL(I,J)=15.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)      
             RL(I,J)=30.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=30.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=30.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)
             RL(I,J)=20.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=20.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=20.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=60.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=60.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=60.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=60.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=60.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=60.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=90.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             SS(I,J)=SS(I,J)-21.0D0*RL(I,J)/(1024.0D0*CSQR**5)
             END DO
             END DO
             ENDIF
c 14'th order term:
             IF(REL(6))THEN
              write(6,*)' P^14 PART : ', REL(6)
             DO I=1,NDIM
             DO J=1,NDIM
             RL(I,J)=      T1(N(1,I),N(1,J),14,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),0,3) 
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),14,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=      T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),14,3) +RL(I,J)
             RL(I,J)=7.0D0*T1(N(1,I),N(1,J),12,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=7.0D0*T1(N(1,I),N(1,J),12,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=7.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),12,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=7.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),12,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=7.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),12,3) +RL(I,J)
             RL(I,J)=7.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),12,3) +RL(I,J)
             RL(I,J)=21.0D0*T1(N(1,I),N(1,J),10,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)      
             RL(I,J)=21.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),10,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)      
             RL(I,J)=21.0D0*T1(N(1,I),N(1,J),10,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)      
             RL(I,J)=21.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),10,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)      
             RL(I,J)=21.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),10,3) +RL(I,J)      
             RL(I,J)=21.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),10,3) +RL(I,J)      
             RL(I,J)=35.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=35.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),0,3) +RL(I,J)
             RL(I,J)=35.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=35.0D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),0,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)
             RL(I,J)=35.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=35.0D0*T1(N(1,I),N(1,J),0,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)
             RL(I,J)=42.0D0*T1(N(1,I),N(1,J),10,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=42.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),10,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=42.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),10,3) +RL(I,J)
             RL(I,J)=105.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=105.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=105.0D0*T1(N(1,I),N(1,J),8,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=105.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),8,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=105.0D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)
             RL(I,J)=105.0D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),8,3) +RL(I,J)
             RL(I,J)=140.D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),2,3) +RL(I,J)
             RL(I,J)=140.D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),2,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=140.D0*T1(N(1,I),N(1,J),2,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=210.D0*T1(N(1,I),N(1,J),6,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             RL(I,J)=210.D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),4,2)*
     &                     T1(N(3,I),N(3,J),6,3) +RL(I,J)
             RL(I,J)=210.D0*T1(N(1,I),N(1,J),4,1)*
     &                     T1(N(2,I),N(2,J),6,2)*
     &                     T1(N(3,I),N(3,J),4,3) +RL(I,J)
             SS(I,J)=SS(I,J)-33.0D0*RL(I,J)/(2048.0D0*CSQR**6)
             END DO
             END DO
             ENDIF


         RETURN
         END
