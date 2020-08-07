C UTEP Electronic Structure Lab (2020)
C
       SUBROUTINE GETGRP(NAME)
C ORIGINAL VERSION BY MARK R PEDERSON 1995-1996
       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*3 NAME
       DIMENSION R(3,3,121),PXYZ(3,3,6),RXYZ(3,3,3)
       SAVE
C
       IF (NAME .EQ. 'GRP') RETURN
C
C PXYZ(1): XYZ=>YXZ
C PXYZ(2): XYZ=>ZYX
C PXYZ(3): XYZ=>XZY
C PXYZ(4): XYZ=>YZX
C PXYZ(5): XYZ=>ZXY
C PXYZ(6): XYZ=>XYZ
C
       DO IP=1,3
        DO I=1,3
         DO J=1,3
          PXYZ(J,I,IP  )=0.0D0
          PXYZ(J,I,IP+3)=0.0D0
          RXYZ(J,I,IP  )=0.0D0
         END DO
         RXYZ(I,I,IP  )=1.0D0
        END DO
        RXYZ(IP,IP,IP)=-1.0D0
        IF (IP.EQ.1) THEN
         IX=1
         IY=2
         IZ=3
        ELSE IF (IP.EQ.2) THEN
         IX=3
         IY=1
         IZ=2
        ELSE IF (IP.EQ.3) THEN
         IX=2
         IY=3
         IZ=1 
        END IF
        PXYZ(IY,IX,IP)=1.0D0
        PXYZ(IX,IY,IP)=1.0D0
        PXYZ(IZ,IZ,IP)=1.0D0
       END DO
C
       PXYZ(1,2,4)=1.0D0
       PXYZ(2,3,4)=1.0D0
       PXYZ(3,1,4)=1.0D0
       PXYZ(1,3,5)=1.0D0
       PXYZ(2,1,5)=1.0D0
       PXYZ(3,2,5)=1.0D0
       PXYZ(1,1,6)=1.0D0
       PXYZ(2,2,6)=1.0D0
       PXYZ(3,3,6)=1.0D0
       MGP=1
       DO I=1,3
        DO J=1,3
         R(J,I,1)=PXYZ(J,I,6)
        END DO
       END DO
C
C START WITH GROUPS
C
       IF(NAME.EQ.'C2V')THEN
        MGP=6
        DO IGP=1,4
         DO I=1,3
          DO J=1,3
           R(J,I,IGP)=PXYZ(J,I,MGP)
          END DO
         END DO
        END DO
        R(1,1,2)=-R(1,1,2)
        R(2,2,3)=-R(2,2,3)
        R(1,1,4)=-R(1,1,4)
        R(2,2,4)=-R(2,2,4)
        MGP=4
       END IF


       IF(NAME.EQ.'D2D')THEN
        MGP=6
        DO IGP=1,2
         DO I=1,3
          DO J=1,3
           R(J,I,IGP)=PXYZ(J,I,MGP)
          END DO
         END DO
         MGP=MGP-5
        END DO
        DO IGP=3,5
         DO I=1,3
          DO J=1,3
           R(J,I,IGP)=-PXYZ(J,I,6)
          END DO
         END DO
        END DO
        R(2,2,3)=-R(2,2,3)
        R(3,3,4)=-R(3,3,4)
        R(1,1,5)=-R(1,1,5)
        MGP=5
       END IF

       IF(NAME.EQ.'D2H')THEN
        DO IGP=1,4
         DO I=1,3
          DO J=1,3
           R(J,I,IGP)=PXYZ(J,I,6)
          END DO
         END DO
        END DO
         DO I=1,3
           R(I,I,I+1)=-R(I,I,I+1)
         END DO

        DO IGP=5,8
         DO I=1,3
          DO J=1,3
           R(J,I,IGP)=-PXYZ(J,I,6)
          END DO
         END DO
        END DO

        R(3,3,5)=-R(3,3,5)
        R(2,2,6)=-R(2,2,6)
        R(1,1,7)=-R(1,1,7)
        MGP=8
       END IF

       IF(NAME.EQ.'C3V')THEN
        MGP=6
        DO IGP=1,6
         DO I=1,3
          DO J=1,3
           R(J,I,IGP)=PXYZ(J,I,MGP)
          END DO
         END DO
         MGP=MGP-1
        END DO
        MGP=6
       END IF
C
       IF (NAME.EQ.'IH') THEN
        DO IGP=1,3
         MGP=MGP+1
         DO I=1,3
          DO J=1,3
           R(J,I,MGP)=RXYZ(J,I,IGP)
          END DO
         END DO
        END DO
        CALL CLSGRP(MGP,R) 
        MGP=MGP+1
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=PXYZ(J,I,4)
         END DO
        END DO
        CALL CLSGRP(MGP,R) 
        A= 0.3090169943749474D0     
        B= 0.5000000000000000D0
        C= 0.8090169943749474D0  
        MGP=MGP+1
        R(1,1,MGP)= A
        R(1,2,MGP)= B
        R(1,3,MGP)=-C
        R(2,1,MGP)=-B
        R(2,2,MGP)= C
        R(2,3,MGP)= A
        R(3,1,MGP)= C
        R(3,2,MGP)= A
        R(3,3,MGP)= B 
       END IF
C
       IF ((NAME.EQ.'TD').OR.(NAME.EQ.'OH')) THEN
        DO I=1,3
         DO J=1,3
          R(J,I,2)=PXYZ(J,I,1)
         END DO
        END DO
        DO I=1,3
         DO J=1,3
          R(J,I,3)=PXYZ(J,I,4)
         END DO
        END DO
        DO I=1,3
         DO J=1,3
          R(J,I,4)=0.0D0       
         END DO
         R(I,I,4)=-1.0D0
        END DO
        R(1,1,4)=1.0D0
        MGP=4
C
        IF (NAME.EQ.'OH') THEN
         DO I=1,3
          DO J=1,3
           R(J,I,5)=0.0D0
          END DO
          R(I,I,5)=-1.0D0
         END DO
         MGP=5
        END IF
       END IF 
C
C CHECK FOR REFLECTIONS:
C
       IF ((NAME.EQ.'X').OR.(NAME.EQ.'Y').OR.(NAME.EQ.'Z')) THEN
        IF(NAME.EQ.'X')IC=1
        IF(NAME.EQ.'Y')IC=2
        IF(NAME.EQ.'Z')IC=3
        MGP=MGP+1
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=0.0D0
         END DO
         R(I,I,MGP)=1.0D0
        END DO
        R(IC,IC,MGP)=-1.0D0
       END IF
C
       IF ((NAME(1:2).EQ.'XY').OR.(NAME(1:2).EQ.'YX').OR.
     &     (NAME(1:2).EQ.'XZ').OR.(NAME(1:2).EQ.'ZX').OR.
     &     (NAME(1:2).EQ.'YZ').OR.(NAME(1:2).EQ.'ZY')) THEN
        IF(NAME(1:1).EQ.'X')I1=1
        IF(NAME(1:1).EQ.'Y')I1=2
        IF(NAME(1:1).EQ.'Z')I1=3
        IF(NAME(2:2).EQ.'X')I2=1
        IF(NAME(2:2).EQ.'Y')I2=2
        IF(NAME(2:2).EQ.'Z')I2=3
        MGP=MGP+1
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=0.0D0
         END DO
         R(I,I,MGP)=1.0D0
        END DO
        R(I1,I1,MGP)=-1.0D0
        MGP=MGP+1
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=0.0D0
         END DO
         R(I,I,MGP)=1.0D0
        END DO
        R(I2,I2,MGP)=-1.0D0
        IF ((NAME(3:3).EQ.'X').OR.(NAME(3:3).EQ.'Y').OR.
     &      (NAME(3:3).EQ.'Z')) THEN
         IF(NAME.EQ.'X')IC=1
         IF(NAME.EQ.'Y')IC=2
         IF(NAME.EQ.'Z')IC=3
         MGP=MGP+1
         DO I=1,3
          DO J=1,3
           R(J,I,MGP)=0.0D0
          END DO
          R(I,I,MGP)=1.0D0
         END DO
         R(IC,IC,MGP)=-1.0D0
        END IF
       END IF
C 
       IF ((NAME.EQ.'D4H').OR.(NAME.EQ.'C4V')) THEN
        DO MGP=1,3
         DO I=1,3
          DO J=1,3
           R(J,I,MGP)=0.0D0
          END DO
          R(I,I,MGP)=1.0D0
         END DO
         IF (MGP.GT.1) R(MGP-1,MGP-1,MGP)= -1.0D0
        END DO
        MGP=4
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=PXYZ(J,I,1)
         END DO
        END DO
        IF (NAME.EQ.'D4H') THEN
         MGP=MGP+1
         DO I=1,3
          DO J=1,3
           R(J,I,MGP)=0.0D0
          END DO
          R(I,I,MGP)=1.0D0
         END DO
         R(3,3,MGP)= -1.0D0
        END IF
       END IF
C
       CALL CLSGRP(MGP,R)
       OPEN(60,FILE='GRPMAT',FORM='FORMATTED',STATUS='UNKNOWN')
       WRITE(60,*) MGP,' ',NAME
       DO IGP=1,MGP
        DO I=1,3
         WRITE(60,60)(R(J,I,IGP),J=1,3)
        END DO
        WRITE(60,*)' '
       END DO
 60    FORMAT(' ',3G25.16)
       CLOSE(60)
       RETURN
       END
