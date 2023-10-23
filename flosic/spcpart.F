C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE SPCPART(MAXANG,NANG,RANG,WANG,RRAD,WRAD,RATM,IFATM,
     &                    INATM,NNTB,RNTB,IFNTB,RCDTAB,
     &                    MODPART,NDECAY,STOMSH)
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:58 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MAXANG, NANG, IFATM, INATM, NNTB, IFNTB, MODPART,
     & NDECAY, I, IANG, INUC, JNUC
       REAL*8 :: RANG , WANG, RRAD, WRAD, RATM, RNTB, RCDTAB, STOMSH,
     & POTF, POTMX, R1, R1RC, R2, R2RC, RATIO, RHOF, RSAV, RVEC, SUM,
     & W, WPRD, WPROD, X
        SAVE
        DIMENSION RANG(3,MAXANG),WANG(MAXANG),RATM(3)
        DIMENSION RCDTAB(MX_CNT,MX_CNT)
        DIMENSION RNTB(3,MX_CNT),STOMSH(6,MAXANG,2),IFNTB(MX_CNT)
        DIMENSION RVEC(3),WPROD(MX_CNT),RSAV(MX_CNT)
C
C CALCULATION OF PARTITION FUNCTION 
C
        DO 400 IANG=1,NANG
         RVEC(1)= RATM(1)+RRAD*RANG(1,IANG)
         RVEC(2)= RATM(2)+RRAD*RANG(2,IANG)
         RVEC(3)= RATM(3)+RRAD*RANG(3,IANG)
         R1= RRAD
         POTMX= 0.0D0
         IF (MODPART .EQ. 1) THEN
C
C RHO*POT/R**NDECAY PARTITIONING
C
          SUM= 0.0D0
          DO INUC= 1,NNTB
           R2= (RVEC(1)-RNTB(1,INUC))**2
     &        +(RVEC(2)-RNTB(2,INUC))**2
     &        +(RVEC(3)-RNTB(3,INUC))**2
           R2= SQRT(R2)
           R2RC= 1.0D0/R2
           CALL RPFIT(IFNTB(INUC),R2,R2RC,RHOF,POTF)
           POTF= ABS(POTF)
           POTMX= MAX(POTMX,POTF)
           DO I= 1,NDECAY
            POTF= POTF*R2RC
           END DO
           SUM= SUM+RHOF*POTF
          END DO
          R1RC= 1.0D0/R1
          CALL RPFIT(IFATM,R1,R1RC,RHOF,POTF)
          RATIO= RHOF*POTF
          DO I= 1,NDECAY
           RATIO= RATIO*R1RC
          END DO
          RATIO= ABS(RATIO)*WANG(IANG)*WRAD/SUM
         ELSE IF (MODPART .EQ. 2) THEN
C
C RHO/R**NDECAY PARTITIONING
C
          SUM= 0.0D0
          DO INUC= 1,NNTB
           R2= (RVEC(1)-RNTB(1,INUC))**2
     &        +(RVEC(2)-RNTB(2,INUC))**2
     &        +(RVEC(3)-RNTB(3,INUC))**2
           R2= SQRT(R2)
           R2RC= 1.0D0/R2
           CALL RPFIT(IFNTB(INUC),R2,R2RC,RHOF,POTF)
           POTF= ABS(POTF)
           POTMX= MAX(POTMX,POTF)
           DO I= 1,NDECAY
            RHOF= RHOF*R2RC
           END DO
           SUM= SUM+RHOF
          END DO
          R1RC= 1.0D0/R1
          CALL RPFIT(IFATM,R1,R1RC,RHOF,POTF)
          RATIO= RHOF
          DO I= 1,NDECAY
           RATIO= RATIO*R1RC
          END DO
          RATIO= ABS(RATIO)*WANG(IANG)*WRAD/SUM
         ELSE 
C
C BECKE PARTITIONING 
C
          DO INUC= 1,NNTB
           R2= (RVEC(1)-RNTB(1,INUC))**2
     &        +(RVEC(2)-RNTB(2,INUC))**2
     &        +(RVEC(3)-RNTB(3,INUC))**2
           RSAV(INUC)= R2
          END DO
          WPRD= 1.0D0
          DO 320 INUC= 1,NNTB
           WPROD(INUC)= 1.0D0
           RSAV(INUC)= SQRT(RSAV(INUC))
           R2= RSAV(INUC)
           CALL RPFIT(IFNTB(INUC),R2,1.0D0/R2,RHOF,POTF)
           POTMX= MAX(POTMX,ABS(POTF))
           IF (INUC .EQ. INATM) GOTO 320
           X= (R1-R2)*RCDTAB(INUC,INATM)
           DO I= 1,NDECAY
            X= X*(1.5D0-0.5D0*X*X)
           END DO
           WPRD= WPRD*0.5D0*(1.0D0-X)
  320     CONTINUE
C
          DO 350 INUC=1,NNTB
           R1= RSAV(INUC)
           DO 340 JNUC=INUC+1,NNTB
            IF ((WPROD(INUC) .LE. 0.0D0) .AND.
     &          (WPROD(JNUC) .LE. 0.0D0)) GOTO 340
            X= (R1-RSAV(JNUC))*RCDTAB(JNUC,INUC)
            DO I= 1,NDECAY
             X= X*(1.5D0-0.5D0*X*X)
            END DO
            W= 0.5D0*(1.0D0-X)
            WPROD(INUC)= WPROD(INUC)*W
            WPROD(JNUC)= WPROD(JNUC)*(1.0D0-W)
  340      CONTINUE
  350     CONTINUE
          SUM= 0.0D0
          DO INUC= 1,NNTB
           SUM= SUM+WPROD(INUC)
          END DO
          RATIO= WANG(IANG)*WRAD*WPRD/SUM
         END IF
C
C STORE POINT
C
         STOMSH(1,IANG,2)= RVEC(1)
         STOMSH(2,IANG,2)= RVEC(2)
         STOMSH(3,IANG,2)= RVEC(3)
         STOMSH(4,IANG,2)= RATIO
         STOMSH(5,IANG,2)= POTMX 
         STOMSH(6,IANG,2)= 0.0D0
  400   CONTINUE
        RETURN
       END
