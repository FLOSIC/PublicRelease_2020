C UTEP Electronic Structure Lab (2020)
       SUBROUTINE EFG
C
C ORIGINALLY WRITTEN BY ALAN JACKSON FALL/1999
C MINOR CHANGES JK 01/2000
C
       use common2,only : RIDT, IFUIDT, NIDENT, ZNUC, SYMATM
       use common3,only : RMAT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:43 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MAXSPH, I, IAT, IATOM, IEV, IH, IID, ISWITCH, J,
     & JATOM, K, KMAX, KMIN, KY, L, MSITES, MTOT, NSPHERE
       REAL*8 :: SYMBOL , CENTER, EFAK, EFGMAX, EFGMIN, EFGRAD, EFGVAL,
     & ETA, PI, R2, R5, RHOFERMI, RREL, RTOL, RVECA, SC, TIME1, TIME2,
     & VZZ, ZCHRG
       SAVE
       PARAMETER (MAXSPH=200) 
       PARAMETER (NMAX=MPBLOCK)
       PARAMETER (NRHOG=6*MXSPN-2)
       PARAMETER (EFAK=9.717366384D0)
C                  Ha/bohr^2 -> V/m^2
C
       LOGICAL LMKFIL,EXIST 
       CHARACTER*80 LINE       
       DIMENSION RREL(3),RVECA(3,MX_GRP)
       DIMENSION EFGRAD(3,3,MAXSPH)
       DIMENSION RHOFERMI(MAXSPH)
       DIMENSION EFGVAL(3,MAXSPH),SC(3)
       DIMENSION CENTER(3,MAXSPH)
C
C
C RETURN IF INPUT FILE DOES NOT EXIST
C
        PRINT '(A)','CALCULATING ELECTRIC FIELD GRADIENTS'
        INQUIRE(FILE='EFG',EXIST=EXIST)
        IF (.NOT.EXIST) THEN
         PRINT '(2A)','EFG: FILE EFG DOES NOT EXIST ',
     &                '--> NOTHING TO DO'
         RETURN
        END IF         
C
C CREATE A STANDARD INPUT FILE IF THE CURRENT INPUT FILE IS EMPTY
C
        CALL GTTIME(TIME1)
        LMKFIL=.TRUE.
        OPEN(74,FILE='EFG',FORM='FORMATTED',STATUS='OLD')
        REWIND(74)
        READ(74,*,END=5,ERR=5) ISWITCH
        IF (ISWITCH .NE. 0) LMKFIL=.FALSE.
    5   CLOSE(74)
C                            
        IF (LMKFIL) THEN
         OPEN(74,FILE='EFG',FORM='FORMATTED',STATUS='OLD')
         REWIND(74)
         WRITE(74,*) '0  auto=0, otherwise user-defined'
         WRITE(74,1010) NIDENT
 1010    FORMAT(' ',I5,' NUMBER OF CENTERS')
         DO IATOM=1,NIDENT 
          WRITE(74,1022)(RIDT(J,IATOM),J=1,3),
     &        (SYMATM(L,IATOM),L=1,10)
 1022     FORMAT(3(1X,F10.5),5X,10A)
         END DO
         CLOSE(74)
        END IF               
C
C READ INPUT FILE
C CENTER CONTAINS THE COORDINATES (X,Y,Z) 
C
        OPEN(74,FILE='EFG',FORM='FORMATTED',STATUS='OLD')
        REWIND(74)
        READ(74,*,END=10) ISWITCH
        READ(74,*,END=10) NSPHERE
        IF (NSPHERE.GT.MAXSPH) THEN
         PRINT *,'EFG: MAXSPH MUST BE AT LEAST: ',NSPHERE
         GOTO 20
        END IF                         
        DO I=1,NSPHERE
         READ(74,33) LINE
         READ(LINE,*)(CENTER(J,I),J=1,3)
         DO J=80,1,-1
          IF (LINE(J:J).NE.' ') THEN
           DO L=1,10
            IH=J-10+L
            SYMATM(L,I)=LINE(IH:IH)
           ENDDO
           GOTO 13
          ENDIF
         ENDDO
   13   CONTINUE
        END DO
        GOTO 30
   10   PRINT *,'EFG: INPUT FILE IS INVALID'
   20   CLOSE(74)
        GOTO 900
   30   CONTINUE
   33   FORMAT(A80)
C                              
C
        PI= 4*ATAN(1.0D0)
          CALL MOSS1(RHOFERMI,EFGRAD,CENTER,NSPHERE)
        RTOL = 0.0001D0
        DO IAT = 1,NSPHERE
         PRINT1022,(CENTER(J,IAT),J=1,3),(SYMATM(L,IAT),L=1,10)
C
C LOOP OVER ALL OTHER SITES
C
        DO IID=1,NIDENT
          ZCHRG = ZNUC(IFUIDT(IID))
          CALL GASITES(1,RIDT(1,IID),MTOT,RVECA,MSITES)
          DO JATOM=1,MTOT
             RREL(1) = CENTER(1,IAT) - RVECA(1,JATOM)
             RREL(2) = CENTER(2,IAT) - RVECA(2,JATOM)
             RREL(3) = CENTER(3,IAT) - RVECA(3,JATOM)
             R2 = RREL(1)**2 + RREL(2)**2 + RREL(3)**2 
             R5 = R2**2.5
             IF(R2.LT.RTOL) GO TO 100
             DO I=1,3
             DO J= 1,3
               EFGRAD(I,J,IAT) = EFGRAD(I,J,IAT) 
     &                - ZCHRG*3.0D0*RREL(I)*RREL(J)/R5
               IF(I.EQ.J) THEN
                  EFGRAD(I,I,IAT) = EFGRAD(I,I,IAT) + ZCHRG*R2/R5
               END IF
             END DO
             END DO
100          CONTINUE
          END DO
         END DO
       END DO
C
C
C PRINT RESULTS
C
        WRITE(74,*)
        WRITE(74,*)
        IEV=1
        DO J=1,NSPHERE
         WRITE(74,*) J,'  ',(SYMATM(L,J),L=1,10),('=',K=1,50)
         WRITE(74,*) RHOFERMI(J)
         DO K=1,3
           WRITE(74,90) (EFGRAD(I,K,J),I=1,3)
         END DO
         WRITE(74,*) ' '
         CALL DIAGSP(3,3,EFGRAD(1,1,J),EFGVAL(1,J),SC,IEV)
         WRITE(74,*) 'AFTER DIAG : '
         WRITE(74,*) (EFGVAL(I,J),I=1,3)
         WRITE(74,*) ' '
         DO K=1,3
           WRITE(74,90) (EFGRAD(I,K,J),I=1,3)
         END DO             
         WRITE(74,*) ' '
C EFG AND ETA
         EFGMAX=-100.d0
         EFGMIN=999999.9d0
         DO K=1,3
          IF(ABS(EFGVAL(K,J)).GT.EFGMAX) THEN
             KMAX=K
             EFGMAX=ABS(EFGVAL(K,J))
          ENDIF
          IF(ABS(EFGVAL(K,J)).LT.EFGMIN) THEN
             KMIN=K
             EFGMIN=ABS(EFGVAL(K,J))
          ENDIF
         ENDDO 
         KY=6-KMAX-KMIN
         VZZ=EFGVAL(KMAX,J)
         ETA=0.0D0
         IF(ABS(VZZ).GT.RTOL) ETA=(EFGVAL(KY,J)-EFGVAL(KMIN,J))/VZZ
         WRITE(74,91) (SYMATM(L,J),L=1,10),VZZ,VZZ*EFAK
         WRITE(74,92) (SYMATM(L,J),L=1,10),ETA
         WRITE(74,*)' '
        END DO
        CLOSE(74)
  90    FORMAT(3F12.6)
  91    FORMAT(1X,'EFG : ',10A,2X,'= ',F12.6,3X,'= ',
     &        F12.6,' *10**21 V/m**2')
  92    FORMAT(1X,'ETA : ',10A,2X,'= ',F12.6)
C
  900   CONTINUE
        CALL GTTIME(TIME2)
        CALL TIMOUT('ELECTRIC FIELD GRADIENTS:          ',TIME2-TIME1)    
        RETURN
       END
