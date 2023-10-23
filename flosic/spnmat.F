C UTEP Electronic Structure Lab (2020)
C
C
C
       SUBROUTINE SPNMAT(SPN, SPT, SPB, THETA, BETA )
       IMPLICIT REAL*8 (A-H,O-Z)
       COMPLEX*16 SPN(2,2,3),SPT(2,2,3),SPB(2,2,3)
       COMPLEX*16 PRD(2,2),ROT(2,2,3),SQM1
C INDICES S1,S2,IX  CORRESPOND TO +1/2 -1/2 AND SX, SY, AND SZ 
C
C |MAJ> =  COS(THETA/2) |+1/2> + EXP(+I*BETA) SIN(THETA/2) |-1/2>
C |MIN> = -EXP(-I*BETA)*SIN(THETA/2) |+1/2> + COS(THETA/2) |-1/2>
C
C THIS SUBROUTINE CALCULATES:
C <MAJ| S_X  | MAJ >    SPN(1,1,IX)  IX=1,2, OR 3 FOR X Y OR Z...
C <MAJ| S_X  | MIN >    SPN(1,2,IX) 
C <MIN| S_X  | MAJ >    SPN(2,1,IX)
C <MIN| S_X  | MIN >    SPN(2,2,IX)
C
C SPB = (D SPN) / (D THETA)
C SPT = (D SPN) / (D BETA )
C      write(6,*)'THETA,BETA:',THETA,BETA
       THE=THETA/2.0D0
       BET=BETA
       ROT(1,1,1)=DCMPLX(COS(THE),0.0D0)
       ROT(1,2,1)=DCMPLX(SIN(THE)*COS(BET) ,SIN(THE)*SIN(BET))
       ROT(2,1,1)=DCMPLX(-SIN(THE)*COS(BET),SIN(THE)*SIN(BET))
       ROT(2,2,1)=DCMPLX(COS(THE),0.0D0)
C DERIVATIVES WRT THETA:
       ROT(1,1,2)=0.5D0*DCMPLX(-SIN(THE),0.0D0)
       ROT(1,2,2)=0.5D0*DCMPLX( COS(THE)*COS(BET), COS(THE)*SIN(BET))
       ROT(2,1,2)=0.5D0*DCMPLX(-COS(THE)*COS(BET), COS(THE)*SIN(BET))
       ROT(2,2,2)=0.5D0*DCMPLX(-SIN(THE),0.0D0)
C DERIVATIVES WRT BETA: 
       ROT(1,1,3)=DCMPLX(0.0D0   ,0.0D0)
       ROT(1,2,3)=DCMPLX(-SIN(THE)*SIN(BET),SIN(THE)*COS(BET))
       ROT(2,1,3)=DCMPLX( SIN(THE)*SIN(BET),SIN(THE)*COS(BET))
       ROT(2,2,3)=DCMPLX(0.0D0   ,0.0D0)
C CHECK FOR UNITARY MATRIX:
       ERROR=0.0D0
       DO I=1,2
       DO J=1,2
           IF(I.NE.J)DIJ=0.0D0
           IF(I.EQ.J)DIJ=1.0D0
           SQM1=0.0D0
           DO K=1,2
           SQM1=SQM1+DCONJG(ROT(K,I,1))*ROT(K,J,1)
           END DO
           ERROR=ERROR+ABS(SQM1-DIJ)
       END DO
       END DO
C      write(6,*)'UNITARITY ERROR:',ERROR
       SPN(1,1,1)= 0.0D0
       SPN(2,2,1)= 0.0D0
       SPN(1,2,1)= 0.5D0
       SPN(2,1,1)= 0.5D0
C
       SPN(1,1,2)= 0.0D0
       SPN(2,2,2)= 0.0D0
       SPN(1,2,2)= DCMPLX(0.0,-0.5)
       SPN(2,1,2)= DCMPLX(0.0, 0.5)
C
       SPN(1,1,3)= 0.5D0
       SPN(2,2,3)=-0.5D0
       SPN(1,2,3)= 0.0D0
       SPN(2,1,3)= 0.0D0
         DO IX=1,3
C        write(6,*)IX
         DO I=1,2
C        PRINT 10,(SPN(I,J,IX),J=1,2)
         END DO
         END DO
C CHECK COMMUTATION RELATIONS:
        SQM1=DCMPLX(0.0D0,1.0D0)
         DO IT=1,3
            IF(IT.EQ.1)THEN
                IX=1
                JX=2
                KX=3
            ELSE IF(IT.EQ.2)THEN
                IX=2
                JX=3
                KX=1
            ELSE
                IX=3
                JX=1
                KX=2
            END IF
            DO K=1,2
            DO L=1,2
              PRD(K,L)=DCMPLX(0.0D0,0.0D0)
              DO J=1,2
              PRD(K,L)=PRD(K,L)+SPN(K,J,IX)*SPN(J,L,JX)/SQM1
              PRD(K,L)=PRD(K,L)-SPN(K,J,JX)*SPN(J,L,IX)/SQM1
              END DO
            END DO
            END DO
C         write(6,*)IX,' * ',JX,' ='
C           DO I=1,2
C           PRINT 10,(PRD(I,J   ),J=1,2)
C           END DO
         END DO
C CREATE GENERALIZED SPIN-MATRICES:
         DO IX=1,3
           DO I=1,2
           DO J=1,2
           PRD(I,J)=DCMPLX(0.0D0,0.0D0)
             DO K=1,2
             DO L=1,2
             PRD(I,J)=PRD(I,J)+DCONJG(ROT(K,I,1))*SPN(K,L,IX)*ROT(L,J,1) 
             END DO
             END DO
           END DO
           END DO
           DO I=1,2
           DO J=1,2
           SPN(J,I,IX)=PRD(J,I)         
           END DO
           END DO
         END DO
C CHECK NORMS:
            DO IX=1,3
            ERROR=0.0D0
            DO K=1,2
            DO L=1,2
                    IF(K.EQ.L)DLK=0.25D0
                    IF(K.NE.L)DLK=0.00D0
              PRD(K,L)= -DLK                
              DO J=1,2
              PRD(K,L)=PRD(K,L)+SPN(K,J,IX)*SPN(J,L,IX)
              END DO
            ERROR=ERROR+ABS(PRD(K,L))
            END DO
            END DO
C           write(6,*)'IX:',IX,' NORM ERROR:',ERROR
            END DO
C CHECK COMMUTATION RELATIONS FOR GENERALIZED SPIN-MATRICES:
        SQM1=DCMPLX(0.0D0,1.0D0)
          DO IT=1,3
            IF(IT.EQ.1)THEN
                IX=1
                JX=2
                KX=3
            ELSE IF(IT.EQ.2)THEN
                IX=2
                JX=3
                KX=1
            ELSE
                IX=3
                JX=1
                KX=2
            END IF
            ERROR=0.0D0
            DO K=1,2
            DO L=1,2
              PRD(K,L)=-SQM1*SPN(K,L,KX)
              DO J=1,2
              PRD(K,L)=PRD(K,L)+SPN(K,J,IX)*SPN(J,L,JX)
              PRD(K,L)=PRD(K,L)-SPN(K,J,JX)*SPN(J,L,IX)
              END DO
            ERROR=ERROR+ABS(PRD(K,L))
            END DO
            END DO
C         write(6,*)'|',IX,' * ',JX,' -I',KX,'| =',ERROR
C           DO J=1,2
C           PRINT 10,(PRD(I,J   ),I=1,2)
C           END DO
          END DO
C CALCULATE DERIVITIVES:
C D <MU |S| NU>/D(THETA)
C D <MU |S| NU>/D(BETA )
C THETA DERIVATIVES:
         DO IX=1,3
           DO I=1,2
           DO J=1,2
           PRD(I,J)=DCMPLX(0.0D0,0.0D0)
             DO K=1,2
             DO L=1,2
             PRD(I,J)=PRD(I,J)+DCONJG(ROT(K,I,2))*SPN(K,L,IX)*ROT(L,J,1) 
             PRD(I,J)=PRD(I,J)+DCONJG(ROT(K,I,1))*SPN(K,L,IX)*ROT(L,J,2) 
             END DO
             END DO
           END DO
           END DO
           DO I=1,2
           DO J=1,2
           SPT(J,I,IX)=PRD(J,I)         
           END DO
           END DO
         END DO
C BETA  DERIVATIVES:
         DO IX=1,3
           DO I=1,2
           DO J=1,2
           PRD(I,J)=DCMPLX(0.0D0,0.0D0)
             DO K=1,2
             DO L=1,2
             PRD(I,J)=PRD(I,J)+DCONJG(ROT(K,I,3))*SPN(K,L,IX)*ROT(L,J,1) 
             PRD(I,J)=PRD(I,J)+DCONJG(ROT(K,I,1))*SPN(K,L,IX)*ROT(L,J,3) 
             END DO
             END DO
           END DO
           END DO
           DO I=1,2
           DO J=1,2
           SPB(J,I,IX)=PRD(J,I)         
           END DO
           END DO
         END DO
 10      FORMAT(2G15.6,5X,2G15.6)
       END 
