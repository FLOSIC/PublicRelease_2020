C UTEP Electronic Structure Lab (2020)
C> @param[in] NDERV
C> @param[out] IUPDAT
C> @param[out] IALLCTR
C> @param[in] NPTS
C> @param[in] PTS
C> @param[in] IFNCT
C> @param[out] GRAD
C
C *****************************************************************
C ORIGINAL VERSION BY MARK R PEDERSON (1987)
C GORBDRV VERSION DIRK POREZAG AUGUST 1994
C GORBDRV CALCULATES THE CONTRIBUTION OF ONE ATOM TO THE
C OCCUPIED WAVEFUNCTIONS UP TO THE SECOND DERIVATIVES
C ON A MESH OF POINTS
C NDERV DETERMINES THE MAXIMUM DERIVATIVE (0 TO 2) 
C
       SUBROUTINE GORBDRV(NDERV,IUPDAT,IALLCTR,NPTS,PTS,IFNCT,GRAD)
       use debug1
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX
C       use common7,only : GAUSS_CUT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:48 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NDERV, NPTS, IFNCT, I, I_BARE, IANG, IEND, IFC, IGR,
     & IPTS, ISHIFT, ISIZE, J, JGR, K, L1, LMAX, LMAX1, NB, NGRAD
       REAL*8 :: PTS , GRAD, ALP, ALP2, ALP2R, CONTMP, CONTR, GAUS,
     & GEPS, PANG, REXPON, RSQR
       SAVE
       LOGICAL IUPDAT,ICOUNT,FIRST,IALLCTR,IBARCTR,ITMPCTR
       DIMENSION IALLCTR(MAX_CON,3),IBARCTR(MAX_CON,3)
       DIMENSION PTS(NSPEED,3),GAUS(NSPEED),RSQR(NSPEED)
       DIMENSION GRAD(NSPEED,10,6,MAX_CON,3)
       DIMENSION PANG(NSPEED,10,10),ALP2R(NSPEED,3)
       DIMENSION ISIZE(3),IEND(3),ISHIFT(3)
C
       DATA ISIZE /1,3,6/
       DATA IEND  /1,4,10/
       DATA ISHIFT/0,1,4/
       DATA FIRST/.TRUE./
C
C      STORAGE IN GRAD:
C      
C      GRAD( A     ,     B     ,     C     ,      D     ,     E )
C            ^           ^           ^            ^
C          POINT   1= WAVEFNCT    FUNCTION    CONTRACTED   ANGULAR
C                  2= GRAD X      (SAME L)      ORBITAL    MOMENTUM
C                  3= GRAD Y                              
C                  4= GRAD Z
C                  5= GRAD XX
C                  6= GRAD YY
C                  7= GRAD ZZ
C                  8= GRAD XY
C                  9= GRAD XZ
C                 10= GRAD YZ
C
C START
C
       IF ((NDERV.LT.0).OR.(NDERV.GT.2)) THEN
        write(6,*)'GORBDRV: INVALID VALUE FOR NDERV: ',NDERV
        CALL STOPIT
       END IF
       LMAX=LSYMMAX(IFNCT)
       LMAX1=LMAX+1
       IF (FIRST) THEN
        GEPS= EXP(-CUTEXP)
        IF (DEBUG) write(6,*)'IN GORBDRV: GAUSS_CUT NOT IN USE' 
        FIRST=.FALSE.
       END IF
       IF(NPTS.GT.NSPEED)THEN
        write(6,*)'GORBDRV: NSPEED MUST BE AT LEAST ',NPTS
        CALL STOPIT
       END IF
C
C GET ANGULAR PART AND SQUARE OF RADIUS
C
       DO IPTS=1,NPTS
        PANG(IPTS,1,1) = 1.0D0
        PANG(IPTS,1,5) = PTS(IPTS,1)*PTS(IPTS,1)
        PANG(IPTS,1,6) = PTS(IPTS,2)*PTS(IPTS,2)
        PANG(IPTS,1,7) = PTS(IPTS,3)*PTS(IPTS,3)    
        RSQR(IPTS)=PANG(IPTS,1,5)+PANG(IPTS,1,6)+PANG(IPTS,1,7)
       END DO
       IF (LMAX.GT.0) THEN
        DO IPTS=1,NPTS
         PANG(IPTS,1,2) = PTS(IPTS,1)
         PANG(IPTS,1,3) = PTS(IPTS,2)
         PANG(IPTS,1,4) = PTS(IPTS,3)
        END DO  
        IF (LMAX.GT.1) THEN
         DO IPTS=1,NPTS
          PANG(IPTS,1,8) = PTS(IPTS,1)*PTS(IPTS,2)
          PANG(IPTS,1,9) = PTS(IPTS,1)*PTS(IPTS,3) 
          PANG(IPTS,1,10)= PTS(IPTS,2)*PTS(IPTS,3) 
         END DO 
        END IF
       END IF 
C
C INITIALIZE POLYNOMIAL DERIVATIVES
C
       NGRAD=1
       IF (NDERV.EQ.0) GOTO 50
       IF (NDERV.EQ.1) NGRAD=4 
       IF (NDERV.EQ.2) NGRAD=10
C
       DO IANG=1,IEND(LMAX1)
        DO IGR=2,NGRAD
         DO IPTS=1,NPTS
          PANG(IPTS,IGR,IANG)=0.0D0
         END DO        
        END DO        
       END DO        
C
C DEFINE NONZERO ELEMENTS OF POLYNOMIAL DERIVATIVE MATRIX
C
C FIRST DERIVATIVES OF P FUNCTIONS
C
       IF (LMAX.GT.0) THEN
        DO IPTS=1,NPTS
         PANG(IPTS,2,2)=   1.0D0
         PANG(IPTS,3,3)=   1.0D0
         PANG(IPTS,4,4)=   1.0D0
        END DO
C          
C FIRST DERIVATIVES OF D FUNCTIONS
C
        IF (LMAX.GT.1) THEN
         DO IPTS=1,NPTS
          PANG(IPTS,2,5)=   2*PTS(IPTS,1)
          PANG(IPTS,2,8)=   PTS(IPTS,2)
          PANG(IPTS,2,9)=   PTS(IPTS,3)
          PANG(IPTS,3,6)=   2*PTS(IPTS,2)
          PANG(IPTS,3,8)=   PTS(IPTS,1)
          PANG(IPTS,3,10)=  PTS(IPTS,3)
          PANG(IPTS,4,7)=   2*PTS(IPTS,3)
          PANG(IPTS,4,9)=   PTS(IPTS,1)
          PANG(IPTS,4,10)=  PTS(IPTS,2)
         END DO
         IF (NDERV.EQ.1) GOTO 50
C
C SECOND DERIVATIVES OF D FUNCTIONS
C
         DO IPTS=1,NPTS
          PANG(IPTS,5,5)=   2.0D0
          PANG(IPTS,6,6)=   2.0D0
          PANG(IPTS,7,7)=   2.0D0
          PANG(IPTS,8,8)=   1.0D0
          PANG(IPTS,9,9)=   1.0D0
          PANG(IPTS,10,10)= 1.0D0
         END DO
        END IF
       END IF
   50  CONTINUE
C
C INITIALIZE GRAD AND IALLCTR
C
       DO L1=1,LMAX1
        DO NB=1,N_CON(L1,IFNCT)
         DO IANG=1,ISIZE(L1)
          DO IGR= 1,NGRAD
           DO IPTS=1,NPTS
            GRAD(IPTS,IGR,IANG,NB,L1)= 0.0D0
           END DO
          END DO
         END DO
        END DO
       END DO
C
       DO L1=1,LMAX1
        DO NB=1,N_CON(L1,IFNCT)
         IALLCTR(NB,L1)= .FALSE.
        END DO
       END DO
C
C LOOP OVER ALL BARE GAUSSIANS
C
       IUPDAT=.FALSE.
       DO 100 I_BARE=1,N_BARE(IFNCT)
        ALP=BFALP(I_BARE,IFNCT)
        ALP2=2*ALP
C
C CREATE ICOUNT AND EXPONENTIALS, CHECK IF UPDATE NECESSARY
C
        ICOUNT=.FALSE.
        DO IPTS=1,NPTS
         REXPON=ALP*RSQR(IPTS)
         IF (REXPON.LT.CUTEXP) THEN
          ICOUNT=.TRUE.
          GAUS(IPTS)=EXP(-REXPON)
         ELSE
          GAUS(IPTS)=0.0D0
         END IF          
        END DO
        IF (.NOT.ICOUNT) GOTO 100
        IUPDAT=.TRUE.
C
C CREATE ARRAYS CONTAINING 2*ALPHA*(X,Y,Z)
C
        IF (NDERV.GE.1) THEN
         DO IGR=1,3
          DO IPTS=1,NPTS
           ALP2R(IPTS,IGR)=ALP2*PTS(IPTS,IGR)
          END DO
         END DO
        END IF
C
C CREATE IBARCTR
C
        DO L1=1,LMAX1
         DO NB=1,N_CON(L1,IFNCT)
          CONTMP= ABS(BFCON(I_BARE,NB,L1,IFNCT))
          ITMPCTR= .FALSE.
          DO IPTS=1,NPTS
           IF (CONTMP*GAUS(IPTS) .GT. GEPS) ITMPCTR= .TRUE.
          END DO
          IBARCTR(NB,L1)= ITMPCTR
          IF (ITMPCTR) IALLCTR(NB,L1)= .TRUE.
         END DO
        END DO
C
C WAVEFUNCTION
C
        DO L1=1,LMAX1
         DO NB=1,N_CON(L1,IFNCT)
          IF (IBARCTR(NB,L1)) THEN
           CONTR=BFCON(I_BARE,NB,L1,IFNCT)
           DO IANG=1,ISIZE(L1)
            IFC=IANG+ISHIFT(L1)
            DO IPTS=1,NPTS
             GRAD(IPTS,1,IANG,NB,L1)=GRAD(IPTS,1,IANG,NB,L1)
     &       +PANG(IPTS,1,IFC)*GAUS(IPTS)*CONTR
            END DO
           END DO
          END IF
         END DO
        END DO
        IF (NDERV.EQ.0) GOTO 100
C
C FIRST DERIVATIVES
C
        DO L1=1,LMAX1
         DO NB=1,N_CON(L1,IFNCT)
          IF (IBARCTR(NB,L1)) THEN
           CONTR=BFCON(I_BARE,NB,L1,IFNCT)
           DO IANG=1,ISIZE(L1)
            IFC=IANG+ISHIFT(L1)
            DO IGR=1,3
             I=IGR+1
             DO IPTS=1,NPTS
              GRAD(IPTS,I,IANG,NB,L1)=GRAD(IPTS,I,IANG,NB,L1)
     &        +(PANG(IPTS,I,IFC)
     &         -PANG(IPTS,1,IFC)*ALP2R(IPTS,IGR))
     &        *GAUS(IPTS)*CONTR
             END DO
            END DO
           END DO
          END IF
         END DO
        END DO
        IF (NDERV.LE.1) GOTO 100
C
C SECOND DERIVATIVES (XX,YY,ZZ)
C
        DO L1=1,LMAX1
         DO NB=1,N_CON(L1,IFNCT)
          IF (IBARCTR(NB,L1)) THEN
           CONTR=BFCON(I_BARE,NB,L1,IFNCT)
           DO IANG=1,ISIZE(L1)
            IFC=IANG+ISHIFT(L1)
            DO IGR=1,3
             I=IGR+1
             J=IGR+4
             DO IPTS=1,NPTS
              GRAD(IPTS,J,IANG,NB,L1)=GRAD(IPTS,J,IANG,NB,L1)
     &        +(PANG(IPTS,J,IFC)
     &         -PANG(IPTS,I,IFC)*ALP2R(IPTS,IGR)*2
     &         +PANG(IPTS,1,IFC)*(ALP2R(IPTS,IGR)**2-ALP2))
     &        *GAUS(IPTS)*CONTR
             END DO
            END DO
           END DO
          END IF
         END DO
        END DO
C
C SECOND DERIVATIVES (XY,XZ,YZ)
C
        DO L1=1,LMAX1
         DO NB=1,N_CON(L1,IFNCT)
          IF (IBARCTR(NB,L1)) THEN
           CONTR=BFCON(I_BARE,NB,L1,IFNCT)
           DO IANG=1,ISIZE(L1)
            IFC=IANG+ISHIFT(L1)
            DO IGR=1,2
             I=IGR+1
             DO JGR=I,3
              J=JGR+1
              K=I+J+3
              DO IPTS=1,NPTS
               GRAD(IPTS,K,IANG,NB,L1)=GRAD(IPTS,K,IANG,NB,L1)
     &         +(PANG(IPTS,K,IFC)
     &          -PANG(IPTS,J,IFC)*ALP2R(IPTS,IGR)
     &          -PANG(IPTS,I,IFC)*ALP2R(IPTS,JGR)
     &          +PANG(IPTS,1,IFC)*ALP2R(IPTS,IGR)*ALP2R(IPTS,JGR))
     &         *CONTR*GAUS(IPTS)
              END DO
             END DO
            END DO
           END DO
          END IF
         END DO
        END DO
  100  CONTINUE
       RETURN
       END
