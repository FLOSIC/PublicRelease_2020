C UTEP Electronic Structure Lab (2020)
       SUBROUTINE CALLDECOMP(LNEWWF,INFIL)
       use common2,only : RIDT, IFUIDT, NIDENT, NCNT, ZNUC, E_UP, E_DN,
     &   ISPN, NSPN, SYMATM
       use common3,only : RMAT
       use common4,only : ISPIDT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:36 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MAXSPH, I, IATOM, IDENT, IDUM, IFAIL, ISIZE, ISWITCH,
     & ITDUM, J, JATOM, KATOM, LMAX, MSITES, MTOT, NCNTX, NPOW, NSPHERE
       REAL*8 :: SYMBOL , AMIN, CENTER, DD, DUM, E1, E2, E_TOT,
     & EDN_SAVE, EMAX, EMIN, ERRMAX, EUP_SAVE, FWHM, QWRT, RVECA,
     & TIME1, TRDUM, Z
       SAVE         
       PARAMETER (MAXSPH=500)     
       CHARACTER*7 INFIL,FNAME
       CHARACTER*7 INFILN
       CHARACTER*4 OUTFILE(2)
       CHARACTER*3 PTYP
       LOGICAL LMKFIL,LNEWWF,EXIST,WINDOW
       DIMENSION CENTER(6,MAXSPH),RVECA(3,MX_GRP)
       DIMENSION ISIZE(3),MSITES(1)   
C
C RETURN IF INPUT FILE DOES NOT EXIST
C
        PRINT '(A)','DECOMPOSING WAVEFUNCTIONS'
C        INQUIRE(FILE=INFIL,EXIST=EXIST)
C        IF (.NOT. EXIST) THEN
C         PRINT '(4A)','DECOMP: FILE ',INFIL,' DOES NOT EXIST ',
C     &                '--> NOTHING TO DO'
C         RETURN
C        END IF
C
C CREATE A STANDARD INPUT FILE IF THE CURRENT INPUT FILE IS EMPTY
C
        CALL GTTIME(TIME1)
        LMKFIL=.TRUE.
        OPEN(74,FILE=INFIL,FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(74)
        READ(74,*,END=5,ERR=5) ISWITCH
        IF (ISWITCH .NE. 0) LMKFIL=.FALSE.
    5   CLOSE(74)
C
        IF (LMKFIL) THEN
         IF(NIDENT.GT.MAXSPH) THEN
         PRINT *,'DECOMP: MAXSPH MUST BE AT LEAST: ',NIDENT
         RETURN
         ENDIF
         OPEN(74,FILE=INFIL,FORM='FORMATTED',STATUS='UNKNOWN')
         WRITE(74,*) '0  auto=0, otherwise user-defined'
         WRITE(74,*) '1.0D-4 0.020 5 2  ERR,ALPMIN,LMAX,NPOW'
         WRITE(74,*) '0.005  0.001      FWHM,QWRT'
         WRITE(74,1010) NIDENT
 1010    FORMAT(' ',I5,' NUMBER OF CENTERS')
         DO IATOM=1,NIDENT
          DO J=1,3
           CENTER(J,IATOM)=RIDT(J,IATOM)
          END DO
          Z=ZNUC(IFUIDT(IATOM))
          CENTER(6,IATOM)=200*ABS(Z)**3
         END DO
         DO IATOM=1,NIDENT
          CENTER(4,IATOM)=0.0D0
          CENTER(5,IATOM)=50.0D0
          DO JATOM=1,NIDENT
           IF (JATOM.NE.IATOM) THEN
            CALL GASITES(1,CENTER(1,JATOM),MTOT,RVECA,MSITES)
            DO KATOM=1,MTOT
             DD=0.5D0*SQRT((CENTER(1,IATOM)-RVECA(1,KATOM))**2+
     &                     (CENTER(2,IATOM)-RVECA(2,KATOM))**2+
     &                     (CENTER(3,IATOM)-RVECA(3,KATOM))**2)
             IF (DD.LT.CENTER(5,IATOM)) CENTER(5,IATOM)=DD
            END DO
           END IF
          END DO
          WRITE(74,1020)(CENTER(J,IATOM),J=1,6)
 1020     FORMAT(3(1X,F10.5),2(1X,F6.2),1X,E10.2,' -0.6 0.0')
         END DO
         WRITE(74,*) 'POSITION,RMIN,RMAX,ALPMAX,ERGMIN,ERGMAX ',
     &               'FOR EACH CENTER'
         CLOSE(74)
        END IF
                OPEN(74,FILE=INFIL)
                WINDOW=.TRUE.
                REWIND(74)
                READ(74,*)IDUM
                READ(74,*) DUM
                READ(74,*) DUM
                READ(74,*)NCNTX
                        EMIN= 1.0D30
                        EMAX=-1.0D30
                        DO IATOM=1,NCNTX
                        READ(74,*)(CENTER(J,IATOM),J=1,6),E1,E2
                           EMIN=MIN(EMIN,E1)
                           EMAX=MAX(EMIN,E2)
                        END DO
                CLOSE(74)
C
       IF(NSPN.EQ.1) THEN
        OUTFILE(1)='DOSO'
        INFILN=INFIL
        INFILN(7:7)='1'
        CALL DECOMP(LNEWWF,INFIL,INFILN,OUTFILE(1))
       ELSE
        IF(E_UP.NE.E_DN) THEN
         OUTFILE(1)='MAJO' 
         OUTFILE(2)='MINO' 
        ELSE
         OUTFILE(1)='AF1O' 
         OUTFILE(2)='AF2O' 
        ENDIF
        EUP_SAVE=E_UP
        EDN_SAVE=E_DN
        E_TOT=E_UP+E_DN
C        LNEWWF=.TRUE.
        LNEWWF=.FALSE.
        E_UP=E_TOT
        E_DN=0
        write(6,*)'DEC: ',E_UP,E_DN
        ITDUM=888
        TRDUM=0.0D0
        CALL WFWIND(EMIN,EMAX,.TRUE.,.FALSE.,IFAIL)
        DO I=1,7
          IF(INFIL(I:I).NE.' ')THEN
            J=I
            INFILN(I:I)=INFIL(I:I)
          END IF
        END DO
        INFILN(J:J)='1'
        CALL DECOMP(LNEWWF,INFIL,INFILN,OUTFILE(1)) 
        E_UP=0
        E_DN=E_TOT
        TRDUM=0.0D0
        ITDUM=889
        CALL WFWIND(EMIN,EMAX,.FALSE.,.TRUE.,IFAIL)
              INFILN(J:J)='2'
        CALL DECOMP(LNEWWF,INFIL,INFILN,OUTFILE(2)) 
        E_UP=EUP_SAVE
        E_DN=EDN_SAVE
        OPEN(99,FILE='EVAL888',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE='EVAL889',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
c        OPEN(99,FILE='EVALUES',STATUS='UNKNOWN')
c        CLOSE(99,STATUS='DELETE')
       ENDIF
C
C READ THE NUMBER OF SPHERE FOR map.dat
C
       OPEN(74,FILE=INFIL,FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(74)
        READ(74,*,END=10) ISWITCH
        READ(74,*,END=10) ERRMAX,AMIN,LMAX,NPOW
        READ(74,*,END=10) FWHM,QWRT
        READ(74,*,END=10) NSPHERE
        GOTO 30
   10   PRINT *,'DECOMP: INPUT FILE IS INVALID'   
   30   CONTINUE
        CLOSE(74)
        IF (NSPHERE.GT.MAXSPH) THEN
         RETURN
        END IF                                         
C 
C
       OPEN(63,FILE='map.dat',STATUS='unknown')
       DO ISPN=1,NSPN
        DO IDENT=1,NSPHERE
        WRITE(FNAME,'(A,I3.3)') OUTFILE(ISPN),IDENT
        IF(ISPIDT(IDENT).EQ.0) PTYP='UPO'
        IF(ISPIDT(IDENT).EQ.1) PTYP='SUP'
        IF(ISPIDT(IDENT).EQ.-1) PTYP='SDN'
        WRITE(63,88) FNAME,(SYMATM(J,IDENT),J=1,10),PTYP
        ENDDO
       ENDDO
       CLOSE(63)
   88  FORMAT(A7,4X,10A,4X,A3)
       RETURN
       END
