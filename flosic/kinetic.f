C UTEP Electronic Structure Lab (2020)
       SUBROUTINE KINETIC 
C WRITTEN BY MARK R. PEDERSON
C
C MAGNETOANISOTROPY ENERGY BY MRP 18-MARCH 1999
C MAJOR CHANGE ON 2-AUGUST 1999
C CALCULATES SPIN-ORBIT COUPLING TERMS
C
       use common2,only : ISPN
       use common3,only : RMAT
       use common5,only : PSI
       use common8,only : REP, N_REP, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:51 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IBAS, IDEBUG, IREP, JBAS, KK, KND, NBAS, NLOW,
     & NMAX
       REAL*8 :: SYMBOL , CLIGHT, EV, EVL, GRAD, H, HA2EV, HA2KEL,
     & HAMKN, OVLP, PI, PSIG, PSIKN, PTS, SC1R, SUM, TEMP
       SAVE
       PARAMETER (MXSPEC=10000)
       PARAMETER (NMAX=MPBLOCK)
       PARAMETER (MDH=MAX_OCC)
C
C  MODE=1   CALCULATE OVERLAPS                -STORE IN HSTOR(I,1)
C  MODE=2   CALCULATE HAMILTONIAN             -STORE IN HSTOR(I,2)
C  MODE=3   CALCULATE KINETIC+NONLOCAL EOP    -STORE IN HSTOR(I,1)
C  MODE=4   CALCULATE KINETIC EOP             -STORE IN HSTOR(I,1)

       LOGICAL ICOUNT,EXIST,FAST,DOIT,RDIT
       LOGICAL LMOM,LDIR,DMOM,LSQMOM
       CHARACTER*20 FNAME
       COMMON/TMP2/PSIG(4,NMAX,MAX_OCC),H(NDH,NDH)
       COMMON/TMP1/PSIKN(NDH,NDH)
     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,ICOUNT(MAX_CON,3)
     &  ,HAMKN(NDH*(NDH+1)/2),EVL(NDH),SC1R(NDH)
     &  ,EV(NDH)
       DIMENSION OVLP(NDH,NDH)
       SAVE
       DATA HA2EV/27.2116D0/
       DATA CLIGHT/137.088D0/
       DATA HA2KEL/315891.1826D0/        !27.2116*1.602/1.38E-04   
       DATA PI /3.141592654D0/
       DATA TEMP/1.0D-4/
       DATA IDEBUG/0/
C
C

       DO I=1,NLOW
          HAMKN(I)=0.0D0
       END DO

C      K4=0
C      DO IREP=1,N_REP
C       NBAS=NS_TOT(IREP)
C        DO IBAS=1,NBAS
C         DO JBAS=IBAS,NBAS
C           K4=K4+1
C           HAMKN(K4)=HSTOR(K4,2)
C         END DO
C        END DO
C      END DO
C      
C
       ISPN=1
       CALL OVERLAP(1)
       CALL OVERLAP3D(2)
C      K4=0
       KIND=0 
       KND=0 
       KK=0 
       SUM=0.0D0
         DO IREP=1,N_REP 
         NBAS=NS_TOT(IREP) 
C
             DO IBAS=1,NBAS
               SC1R(IBAS)=0.0D0
                EVL(IBAS)=0.0D0
              DO JBAS=1,NBAS
               OVLP(JBAS,IBAS)=0.0D0
               H(JBAS,IBAS)   =0.0D0
              END DO
             END DO
C
C      FILL OVLP AND HAM:
C
C           DO IBAS=1,NBAS
C            DO JBAS=IBAS,NBAS
C              K4=K4+1
C               H(JBAS,IBAS)=HSTOR(K4,2)
C              H(JBAS,IBAS)=HAMKN(K4)
C              H(IBAS,JBAS)=H(JBAS,IBAS)
C            END DO
C           END DO
C
C
            DO IBAS=1,NBAS
             DO JBAS=IBAS,NBAS
               KND=KND+1
               H   (JBAS,IBAS)=HSTOR(KND,2)
               H   (IBAS,JBAS)=H(JBAS,IBAS)
               OVLP(JBAS,IBAS)=HSTOR(KND,1)
               OVLP(IBAS,JBAS)=OVLP(JBAS,IBAS)
             END DO
            END DO
C
C        
           CALL DIAGGE(NDH,NBAS,H,OVLP,EVL,SC1R,1)
           DO I=1,NBAS
             WRITE(6,100)I,EVL(I) 
           END DO
         
        END DO   ! IREP
  100  FORMAT(I5,F20.8)
  200  FORMAT(4G20.8)
  300  FORMAT(I5,2F20.8)
       STOP 
       END
