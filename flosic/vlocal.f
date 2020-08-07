C UTEP Electronic Structure Lab (2020)
C
C *******************************************************************
C
       SUBROUTINE VLOCAL(MODE,NPV,IFU,RDIS,VLOC)
C
C CALLED FROM APOTNL
C MODE=1: DETERMINE THE LOCAL PART OF THE NUCLEAR PSEUDOPOTENTIAL
C MODE=2: DETERMINE THE DERIVATIVE OF THE LOCAL PART OF THE NUCLEAR 
C         PSEUDOPOTENTIAL VERSUS THE DISTANCE, DIVIDED BY THE DISTANCE
C
       use debug1
       use common1,only : PSPSYM,BHSALP,BHSCOF,RRADTAB,VLRTAB,NRADTAB
       use common2,only : ZELC
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MODE, NPV, IFU, I, IALP, IPV, MXITER
       REAL*8 :: SYMBOL , RDIS, VLOC, A0, A1, ACCUM, AI, AIA, AIF, B0,
     & B1, FAC, G, GOLD, OLD, R, RDRC, RRC, RTPIRC, SUM, VDRV, X, X2
        SAVE
        PARAMETER(RTPIRC=0.56418958354775628695D0)
        PARAMETER(MXITER=100)
        DIMENSION RDIS(NSPEED),VLOC(NSPEED)
        DIMENSION RDRC(NSPEED),VDRV(2,NSPEED)
        DATA ACCUM/1.0D-14/
C
C ALL-ELECTRON
C
        IF (NPV .GT. NSPEED) THEN
         PRINT *,'VLOCAL: NSPEED MUST BE AT LEAST: ',NPV
         CALL STOPIT
        END IF
        IF ((MODE .NE. 1) .AND. (MODE .NE. 2)) THEN
         PRINT *,'VLOCAL: INVALID MODE'
         CALL STOPIT
        END IF
        DO IPV=1,NPV
         RDRC(IPV)= 1.0D0/RDIS(IPV)
        END DO
        IF (PSPSYM(IFU)(1:3) .EQ. 'ALL') THEN
         IF (MODE .EQ. 1) THEN
          DO IPV=1,NPV
           VLOC(IPV)= -ZELC(IFU)*RDRC(IPV)
          END DO
         ELSE 
          DO IPV=1,NPV
           VLOC(IPV)= +ZELC(IFU)*RDRC(IPV)*RDRC(IPV)*RDRC(IPV)
          END DO
         END IF
C
C BHS
C USE OWN ERROR FUNCTION SINCE NOT F77 STANDARD
C TAYLOR EXPANSION FOR SMALL ARGUMENTS, CONTINUED FRACTION FOR LARGE ONES
C
        ELSE IF ((PSPSYM(IFU)(1:3) .EQ. 'BHS')
     &      .OR. (PSPSYM(IFU)(1:3) .EQ. 'ECP'))THEN
         DO 30 IPV=1,NPV
          R= RDIS(IPV)
          VLOC(IPV)= 0.0D0
          DO IALP=1,2
           X= BHSALP(IALP,IFU)*R
           X2=X*X
           IF (X2 .LE. 4.0D0) THEN
            FAC= 2*RTPIRC*EXP(-X2)
            SUM= 0.0D0
            DO I=1,MXITER
             OLD= SUM
             SUM= SUM+FAC
             IF (SUM .EQ. OLD) GOTO 10
             FAC= FAC*X2/(0.5D0+I)
            END DO
            PRINT *,'VLOCAL: BHS ERROR(1): X2= ',X2
            CALL STOPIT
   10       CONTINUE
            SUM= SUM*X
           ELSE
            GOLD= 0.0D0
            A0= 1.0D0
            B0= 0.0D0
            A1= X2
            B1= 1.0D0
            FAC= 1.0D0
            DO I=1,MXITER
             AI=I
             AIA=AI-0.5D0
             A0=(A1+A0*AIA)*FAC
             B0=(B1+B0*AIA)*FAC
             AIF=AI*FAC
             A1= X2*A0+AIF*A1
             B1= X2*B0+AIF*B1
             IF (A1 .NE. 0.0D0) THEN
              FAC= 1.0D0/A1
              G= B1*FAC
              IF (ABS(G-GOLD) .LE. ACCUM*ABS(GOLD)) GOTO 20
              GOLD=G
             END IF
            END DO
            PRINT *,'VLOCAL: BHS ERROR(2): X2= ',X2
            CALL STOPIT
   20       CONTINUE
            SUM= 1.0D0-RTPIRC*X*EXP(-X2)*G
           END IF
           VLOC(IPV)= VLOC(IPV)+BHSCOF(IALP,IFU)*SUM
          END DO
          IF (MODE .EQ. 1) THEN
           VLOC(IPV)= -ZELC(IFU)*VLOC(IPV)*RDRC(IPV)
          ELSE
           RRC= RDRC(IPV)
           VLOC(IPV)= ZELC(IFU)*RRC*RRC*(VLOC(IPV)*RRC-2*RTPIRC
     &     *(BHSCOF(1,IFU)*BHSALP(1,IFU)*EXP(-(BHSALP(1,IFU)*R)**2)
     &      +BHSCOF(2,IFU)*BHSALP(2,IFU)*EXP(-(BHSALP(2,IFU)*R)**2)))
          END IF
   30    CONTINUE
C
C TAB
C
        ELSE IF (PSPSYM(IFU)(1:3) .EQ. 'TAB') THEN
         IF (MODE .EQ. 1) THEN
          CALL FINTPOL(8,NPV,RDIS,0.1D0,NRADTAB(IFU),2,1,RRADTAB(1,IFU),
     &                 VLRTAB(1,1,IFU),VDRV)
          DO IPV=1,NPV
           VLOC(IPV)= VDRV(1,IPV)*RDRC(IPV)
          END DO
         ELSE
          CALL FINTPOL(8,NPV,RDIS,0.1D0,NRADTAB(IFU),2,2,RRADTAB(1,IFU),
     &                 VLRTAB(1,1,IFU),VDRV)
          DO IPV=1,NPV
           VLOC(IPV)= (VDRV(2,IPV)-VDRV(1,IPV)*RDRC(IPV))*RDRC(IPV)**2
          END DO
         END IF
        END IF
        RETURN
       END
