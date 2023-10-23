C UTEP Electronic Structure Lab (2020)
C
C ******************************************************
C
      SUBROUTINE FFFMTC(MX,NPTS,XX,CHGD,S1) 
C ORIGINAL VERSION BY SUNGGI CHUNG
C MODIFIED BY MARK R PEDERSON (1990)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL FIRST
      DIMENSION F(10,1601),ONEOVRI(5)
      DIMENSION XX(MX),S1(MX,5),CHGD(MX)
      DIMENSION G(20)
      SAVE
      DATA FIRST/.TRUE./
      DATA F(1,1)/0.0D0/
      DATA ROOTPI/1.77245385090551602728D0/
C
C        DIFFERENCE BETWEEN DOUBLE $ SINGLE  ARE
C        REAL*8 CARD  $  FORMAT 55  $  COPOUT WORD SIZE
C
      IF(NPTS.GT.MX)THEN
       write(6,*)'FFFMTC: MX MUST BE AT LEAST: ',NPTS
       CALL STOPIT
      END IF
      IF (FIRST) THEN
       FIRST=.FALSE.
       FACT=1.0D0
       ONEOVER7=1.0D0/7.0D0
       ONEOVER5=1.0D0/5.0D0
       ONEOVER3=1.0D0/3.0D0
       DO 10 I=1,5
        AI=I
        ONEOVRI(I)=1.0D0/AI
        FACT=FACT*ONEOVRI(I)
   10  CONTINUE
       EPS=1.0D-12
       DO 20 I=1,10
        F(I,1)=1.0D0/(2*I-1)
   20  CONTINUE
       T=0.0D0
       DO 200 L=2,1601
        T=(1.0D-2)*(L-1)
        X=2*T
        DEN=39.0D0
        TERM=1.0D0/DEN
        SUM=TERM
        DO 130 I=2,100
         DEN=DEN+2.0D0
         TERM=TERM*X/DEN
         SUM=SUM+TERM
         Q=TERM/SUM
         IF(Q-EPS) 140,140,130
  130   CONTINUE
  140   EX=EXP(-T)
        G(20)=EX*SUM
C
C  USE DOWNWARD RECURSION
C
        DO 150 I=1,19
         K=21-I
         KK=K-1
         G(KK)=(X*G(K)+EX)/(2*KK-1)
  150   CONTINUE
        DO 160 I=1,10
         F(I,L)=G(I)
  160   CONTINUE
  200  CONTINUE
       DO 215 L=1,1601
        F(7,L)=F(7,L)/2.0D0
        F(8,L)=F(8,L)/6.0D0
  215  CONTINUE
      END IF
C
      DO IPTS=1,NPTS
       IF(XX(IPTS) .GT. 16.0D0)THEN
C       POL1= -EXP(-XX(IPTS))
        POL1= -CHGD(IPTS)
        RT=SQRT(XX(IPTS))*ROOTPI
        RECT1=0.5D0/XX(IPTS)
        S1(IPTS,1)=RECT1*(RT+POL1)
        S1(IPTS,2)=RECT1*(S1(IPTS,1)+POL1)
        S1(IPTS,3)=RECT1*(3*S1(IPTS,2)+POL1)
        S1(IPTS,4)=RECT1*(5*S1(IPTS,3)+POL1)
        S1(IPTS,5)=RECT1*(7*S1(IPTS,4)+POL1)
       ELSE
        IT=INT(100*XX(IPTS)+1.5D0)
C       POL1=EXP(-XX(IPTS))
        POL1=CHGD(IPTS)
        RECT1=2*XX(IPTS)
        DT=(1.0D-2)*(IT-1)-XX(IPTS)
        S1(IPTS,5)=F(5,IT)+DT*(F(6,IT)+DT*(F(7,IT)+DT*F(8,IT)))
        S1(IPTS,4)=(S1(IPTS,5)*RECT1+POL1)*ONEOVER7
        S1(IPTS,3)=(S1(IPTS,4)*RECT1+POL1)*ONEOVER5
        S1(IPTS,2)=(S1(IPTS,3)*RECT1+POL1)*ONEOVER3
        S1(IPTS,1)=(S1(IPTS,2)*RECT1+POL1)
       END IF
      END DO
      RETURN
      END 
