C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE OPTRMSH(NRSIZ,IFNCT,RA,RB,OVLTAB,ERRMAX,
     &                    NRAD,RRAD,WRAD)
C
C CREATE OPTIMIZED RADIAL MESH FOR INTEGRATION IN THE INTERVAL (RA,RB)
C
       use common2,only : N_CON, LSYMMAX
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NRSIZ, IFNCT, NRAD, I, ICON, IMN, IOFS, IPT, ITYPE,
     & L, L1, MPT, MUL, NCOMP, NSEARCH, NTRUE
       REAL*8 :: SYMBOL , RA, RB, OVLTAB, ERRMAX, RRAD, WRAD, AA,
     & ALOPT, ALORC, BB, EBST, ERRG, ERRO, ERT1, ERT2, RMAX, RMIN,
     & SUMTAB, TINY, WC, WOVL, WT, XC, XT
        SAVE
        PARAMETER (NTRUE=96)
        PARAMETER (NCOMP=48)
        PARAMETER (NSEARCH=50)
        DIMENSION OVLTAB(MAX_CON,LDIM)
        DIMENSION RRAD(NRSIZ),WRAD(NRSIZ)
        DIMENSION XT(NTRUE),WT(NTRUE),XC(NCOMP),WC(NCOMP)
        DIMENSION SUMTAB(2,MAX_CON,LDIM,2),WOVL(MAX_CON,LDIM)
        LOGICAL OK
        DATA TINY /1.0D-10/
C
C CHECK INPUT PARAMETERS
C
        IF (RA .GE. RB)    RETURN
        IF (RA .LT. 0.0D0) RETURN
        IF (ERRMAX .LT. TINY) THEN
         write(6,*)'OPTRMSH: DESIRED ACCURACY IS TOO SMALL: ',ERRMAX
         write(6,*)'         MUST NOT BE SMALLER THAN: ',TINY
         CALL STOPIT
        END IF
C
C SETUP OF WOVL 
C
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          WOVL(ICON,L1)= 1.0D0/OVLTAB(ICON,L1)
         END DO
        END DO
C
C START OF MESH GENERATION 
C
        NRAD= 0
        RMIN= RA
        RMAX= RB
  100   CONTINUE
C
C CHECK IF THE INTEGRAL IS THE SAME WITH A NCOMP AND A NTRUE POINT MESH
C IF NOT, REDUCE RMAX AND TRY AGAIN
C
         CALL GAUSSP(RMIN,RMAX,NTRUE,XT,WT)
         CALL GAUSSP(RMIN,RMAX,NCOMP,XC,WC)
         CALL OPTRIGR(1,IFNCT,NTRUE,XT,WT,SUMTAB(1,1,1,1))
         CALL OPTRIGR(1,IFNCT,NCOMP,XC,WC,SUMTAB(1,1,1,2))
         DO L= 0,LSYMMAX(IFNCT)
          L1= L+1
          MUL= 4*L+2
          DO ICON= 1,N_CON(L1,IFNCT)
           DO I= 1,2
            ERRG= (SUMTAB(I,ICON,L1,1)-SUMTAB(I,ICON,L1,2))
     &            *WOVL(ICON,L1)*MUL
            IF (ABS(ERRG) .GT. ERRMAX) THEN
             RMAX= 0.5D0*(RMIN+RMAX)
             GOTO 100
            END IF
           END DO
          END DO
         END DO
C
C WE NOW KNOW THAT A NCOMP POINT MESH CAN'T FAIL ...
C NOW TRY MESHES BETWEEN 0 AND NCOMP POINTS
C IF (MPT .EQ. 0), CHECK ONLY SIMPLE GAUSS QUADRATURE
C
         OK= .FALSE.
         DO 400 MPT= 0,NCOMP
          ITYPE= 0
          EBST= ERRMAX
C
C DEAL WITH GAUSS MESH
C
          CALL GAUSSP(RMIN,RMAX,MPT,XC,WC)
          CALL OPTRIGR(1,IFNCT,MPT,XC,WC,SUMTAB(1,1,1,2))
          ERRG= 0.0D0
          DO L= 0,LSYMMAX(IFNCT)
           L1= L+1
           MUL= 4*L+2
           DO ICON= 1,N_CON(L1,IFNCT)
            DO I= 1,2
             AA= (SUMTAB(I,ICON,L1,1)-SUMTAB(I,ICON,L1,2))
     &           *WOVL(ICON,L1)*MUL
             ERRG= MAX(ERRG,ABS(AA))
            END DO
           END DO
          END DO
          IF (MPT .EQ. 0) THEN
           IF (ERRG .LT. ERRMAX) GOTO 500
           GOTO 400
          END IF
C
C CHECK IF GAUSS-LEGENDRE MESH IS ACCURATE ENOUGH
C
          IF (ERRG .LT. EBST) THEN
           OK= .TRUE.
           ITYPE= 1
           EBST= ERRG
           DO IPT= 1,MPT
            XT(IPT)= XC(IPT)
            WT(IPT)= WC(IPT)
           END DO
          END IF
C
C NOW DEAL WITH OPTIMIZED MESH
C NSEARCH OPTIMIZATION STEPS
C
          ERRO= 0.0D0
          ERT1= 1.0D30
          ERT2= 1.0D30
          DO 200 IMN= 1,NSEARCH
           CALL MINIMIZE(IMN,1.0D0,ALOPT,ERRO)
           ALORC= 1.0D0/ALOPT
C
C SWITCH LIMITS TO FOR EXPONENTIAL MESH
C CALL GAUSSP AND TRANSFORM POINTS/WEIGHTS
C
           AA= EXP(-ALOPT*RMAX)
           BB= EXP(-ALOPT*RMIN)
           CALL GAUSSP(AA,BB,MPT,XC,WC)
           DO I= 1,MPT
            WC(I)= WC(I)/(XC(I)*ALOPT)
            XC(I)= -ALORC*LOG(XC(I))
           END DO
C
C GET ERROR FOR TRANSFORMED MESH
C
           CALL OPTRIGR(1,IFNCT,MPT,XC,WC,SUMTAB(1,1,1,2))
           ERRO= 0.0D0
           DO L= 0,LSYMMAX(IFNCT)
            L1= L+1
            MUL= 4*L+2
            DO ICON= 1,N_CON(L1,IFNCT)
             DO I= 1,2
              AA= (SUMTAB(I,ICON,L1,1)-SUMTAB(I,ICON,L1,2))
     &            *WOVL(ICON,L1)*MUL
              ERRO= MAX(ERRO,ABS(AA))
             END DO
            END DO
           END DO
C
C IF WE HAVE A SUCCESSFUL MESH, STORE IT IN XT, WT
C EBST CONTAINS ERROR OF BEST MESH SO FAR
C
           IF (ERRO .LT. EBST) THEN
            OK= .TRUE.
            ITYPE= 2
            EBST= ERRO
            DO IPT= 1,MPT
             XT(IPT)= XC(IPT)
             WT(IPT)= WC(IPT)
            END DO
           END IF
C 
C IF THE ERROR FOR THE BEST MESH HAS CHANGED LESS THAN 1 PERCENT
C DURING THE LAST 5 ITERATIONS, STOP OPTIMIZING ALOPT
C
           ERT2= MIN(ERT2,ERRO) 
           IF (MOD(IMN,5) .EQ. 0) THEN 
            IF (ERT2 .GT. ERT1*0.99D0) GOTO 300
            ERT1= ERT2
           END IF
  200     CONTINUE
  300     IF (OK) GOTO 410
  400    CONTINUE
  410    CONTINUE
C
C STORE MESH POINTS IN RRAD, WRAD (MULTIPLY WEIGHTS BY R**2)
C
         IF (NRAD+MPT .GT. NRSIZ) THEN
          write(6,*)'OPTRMSH: NRSIZ EXCEEDED: ',NRSIZ,NRAD+MPT
          CALL STOPIT
         END IF
         IF (ITYPE .EQ. 1) THEN
          DO IPT= 1,MPT
           RRAD(NRAD+IPT)= XT(IPT)
           WRAD(NRAD+IPT)= WT(IPT)*XT(IPT)**2
          END DO
         ELSE
          IOFS=NRAD+MPT+1
          DO IPT=1,MPT
           RRAD(IOFS-IPT)= XT(IPT)
           WRAD(IOFS-IPT)= WT(IPT)*XT(IPT)**2
          END DO
         END IF
         NRAD=NRAD+MPT
  500    CONTINUE
C
C THE MESH WAS SUCCESSFUL, UPDATE RMIN AND RMAX
C
         IF (RMAX .NE. RB) THEN
          RMIN= RMAX
          RMAX= RB
          GOTO 100
         END IF
        CONTINUE
        RETURN
       END
