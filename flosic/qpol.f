C UTEP Electronic Structure Lab (2020)
C ****************************************************************
C
       SUBROUTINE QPOL
C   CALCULATES THE DIPOLE MOMENT ONLY IN THE INEQUIVALENT PART OF THE 
C   MOLECULE.  
C
       use common2,only : RIDT, ZNUC, N_POS, NFNCT, IGGA, NSPN, DIPOLE
       use common3,only : RMAT
       use common7,only : MODDEN
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:59 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I1, I2, I_POS, ICNT, IFNCT, IGRP, IPOS, IPTS, IQ,
     & ISHELLA, ISIT, IX, LCNT, LSIT, MSITES, MXPT, NGRAD, NMSH, NSIT,
     & NSUM
       REAL*8 :: SYMBOL , CMIN, COULOMB, DIS, DIST, DMIN, FAC, PHIG,
     & PI, QPOLE, RDS, RHOG, RLOCA, TIMEGRB
       SAVE
       PARAMETER (MAXSPH=500)
       PARAMETER (MAXRAD=1000)
       PARAMETER (MAXANG=200)
C
       LOGICAL ICOUNT,EXIST
       CHARACTER*20 FNAME(4),FORMSTR
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       DIMENSION MSITES(1),RDS(3,MAX_IDENT),QPOLE(3,2),RLOCA(3,MX_GRP)
       DIMENSION DIST(MAX_IDENT)
       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,KRHOG,MXSPN)
     &  ,PHIG(MAX_PTS,2)
C
C
C FOR DENSITY EVALUATIONS, USE OLD SCHEME FOR DENSITY CALCULATION
C OTHERWISE, GET DENSITY FROM COUPOT
C
       PI=4.0D0*DATAN(1.0D0)
       PRINT '(A)',' '
       PRINT '(A)',' DIPOLE FOR THE INEQUIVALENT PART'
C
C READ IN NECESSARY INPUT DATA
C
         FNAME(1)='QPOLE'
         INQUIRE(FILE=FNAME(1),EXIST=EXIST)
C
C DETERMINE IF THERE IS ANYTHING TO DO 
C
        FORMSTR= ' '
        IF (.NOT.EXIST) FORMSTR= ' --> NOTHING TO DO'
        PRINT '(2A)','QPOLE :  ',FORMSTR
        IF (.NOT.EXIST) RETURN
C
C READ INPUT DATA
C
        OPEN(72,FILE=FNAME(1),FORM='FORMATTED',STATUS='OLD')
        REWIND(72)
C
        DO IGRP=1,MX_GRP
         DO IX=1,3
           RLOCA(IX,IGRP)=0.0D0
         END DO
        END DO
C
C NOW: CALCULATE ELECTRONIC DENSITY
C DENSITY WILL BE STORED IN RHOG
C
         NGRAD=1
         MODDEN=2
c
c
          I1=IGGA(1)
          I2=IGGA(2)
          IGGA(1)=0
          IGGA(2)=0
          CALL DENSOLD(TIMEGRB)
          IGGA(1)=I1
          IGGA(2)=I2

C
C UPDATE DATA IN RHOG
C
          DO IPTS=1,NMSH
           RHOG(IPTS,1,1)=RHOG(IPTS,1,1)+RHOG(IPTS,1,NSPN)
          END DO 

          DO IQ=1,2
           DO IX=1,3
             QPOLE(IX,IQ)=0.0D0
           END DO
          END DO

C       GET THE CORRECT WEIGHT AND MESH POINT

        NSUM=0
        DO IPTS=1,NMSH
         CALL GASITES(1,RMSH(1,IPTS),MXPT,RLOCA,MSITES)
         NSIT=MXPT
         NSUM=NSUM+NSIT

         ISHELLA=0
         DO IFNCT=1,NFNCT
          DO I_POS=1,N_POS(IFNCT)
           ISHELLA=ISHELLA+1

            DMIN=1.0D+30
            DO ISIT=1,NSIT
             DIS=(RLOCA(1,ISIT)-RIDT(1,ISHELLA))**2
     &          +(RLOCA(2,ISIT)-RIDT(2,ISHELLA))**2
     &          +(RLOCA(3,ISIT)-RIDT(3,ISHELLA))**2
             DIS=DSQRT(DIS)
             IF (DIS.LT.DMIN) THEN
              LSIT=ISIT 
              DMIN=DIS
             END IF
            END DO
            DIST(ISHELLA)=DMIN
            RDS(1,ISHELLA)=RLOCA(1,LSIT)
            RDS(2,ISHELLA)=RLOCA(2,LSIT)
            RDS(3,ISHELLA)=RLOCA(3,LSIT)
          END DO
         END DO     ! IFNCT
         CMIN=1.0D+30
         LCNT=0
          DO ICNT=1,ISHELLA 
           IF(CMIN.GT.DIST(ICNT)) THEN 
            CMIN=DIST(ICNT)
            LCNT=ICNT 
           END IF
          END DO
          RHOG(IPTS,1,1)=RHOG(IPTS,1,1)*WMSH(IPTS)/NSIT
          IF(IPTS.EQ.1) THEN
             WRITE(72,*)'ISHELLA = ', ISHELLA
             WRITE(72,*)'LCNT = ', LCNT
             WRITE(72,*)'DIST = ', (DIST(I),I=1,ISHELLA)
             WRITE(72,*)(RDS(IX,LCNT),IX=1,3)
          END IF
          DO IX=1,3
           QPOLE(IX,1)=QPOLE(IX,1)+RHOG(IPTS,1,1)
           QPOLE(IX,2)=QPOLE(IX,2)+RDS(IX,LCNT)*RHOG(IPTS,1,1)
          END DO
        END DO
C 
C      CONTRIBUTIONS FROM ATOM CENTERS   
        ISHELLA=0
        DO IFNCT=1,NFNCT
         DO IPOS=1,N_POS(IFNCT)
          ISHELLA=ISHELLA+1
          FAC=ABS(ZNUC(IFNCT))
          DO IX=1,3
           QPOLE(IX,2)=QPOLE(IX,2)-FAC*RIDT(IX,ISHELLA)
          END DO
         END DO
        END DO
       WRITE(72, *)'NSIT = ',NSUM
       WRITE(72, 9010)(QPOLE(IX,1),IX=1,3)
       WRITE(72, 9020)(QPOLE(IX,2),IX=1,3)

 9010  FORMAT('CHARGE :',3(1X,E20.12))
 9020  FORMAT('DIPOLE :',3(1X,F20.12))
       CLOSE(72)
       RETURN
       END
