C UTEP Electronic Structure Lab (2020)
       SUBROUTINE MOLDENDRV(OPTION)
C
C      BY MULLINS (2010)
C      UPDATED BY ULISES REVELES (JUNE 2013)
C
c      This subroutine does the writing to (and updating of) the
C      molden format file that tracks the history of the calculation 
C      and allows molden to graphically display the results. 
C
C      OPTION:
C        XYZ = WRITE ONLY CARTESIAN COORDINATES
C        OPT = WRITE COORDINATES, ENERGY, AND FORCES 
C        MOS = WRITE BASIS SETS AND MO's COEFFICIENTES
C
C     ------------------------------------------------------------------
C
C      --- COMMON VARIABLES ---
C
      use debug1
      use xmol, only : NUM_ATMS,XMOL_LIST,GET_LETTER,AU2ANG 
      use common2,only : ETOTAL,RIDT,NIDENT,NFNCT,NATOMS,N_POS,IFUIDT,
     &                   FTOT,ZNUC,GMAX
      use global_inputs,only : EXCITED1
C
C     --- PARAMETERS ---
C
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:54 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: INUC, M_NUC, MSITES, NNUC
       REAL*8 :: GRAD 
C
C     --- LOCAL VARIABLES ---
C
      LOGICAL EXIST
      CHARACTER*2 LETTER
      CHARACTER*3 OPTION
      CHARACTER*132 INLINE,STRCOMP
C
      INTEGER I,IID,IGEO,IGCNV,IFOR,ITEST,J,K
!      REAL*8 AU2ANG,GRAD
C
!      DATA AU2ANG/0.529177249D0/   !YY 3 more decimal added
C
C     --- LOCAL DYNAMICAL FIELDS ---
C
      INTEGER ALLOCATION
      INTEGER,ALLOCATABLE :: IDINDX(:)
      REAL*8,ALLOCATABLE  :: REQV(:,:),RNUC(:,:)
C
C     ------------------------------------------------------------------
C
C      --- ALLOCATE LOCAL FIELDS ---
C
       IF(NATOMS==0) THEN
         INQUIRE(FILE='XMOL.DAT',EXIST=EXIST)
         IF(EXIST)THEN
           OPEN(19,FILE='XMOL.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
           READ(19,*)NATOMS
           READ(19,*)
           IF(.NOT.ALLOCATED(XMOL_LIST)) ALLOCATE(XMOL_LIST(NATOMS))
           DO I=1,NATOMS
             READ(19,*) XMOL_LIST(I)%ANUM,XMOL_LIST(I)%RX,
     &                XMOL_LIST(I)%RY,  XMOL_LIST(I)%RZ
!YY XMOL.DAT is in Ang. These need to be in Bohr radius.
             XMOL_LIST(I)%RX=XMOL_LIST(I)%RX/AU2ANG
             XMOL_LIST(I)%RY=XMOL_LIST(I)%RY/AU2ANG
             XMOL_LIST(I)%RZ=XMOL_LIST(I)%RZ/AU2ANG
             WRITE(6,*) XMOL_LIST(I)%ANUM,XMOL_LIST(I)%RX,
     &                XMOL_LIST(I)%RY,  XMOL_LIST(I)%RZ
           
           END DO
           CLOSE(19)
         ELSE
           CALL CHECK_INPUTS
           IF(EXCITED1) THEN
             PRINT *,"ERROR: XMOL.DAT File is required for 
     &               Delta-SCF+MOLDEN"
             CALL STOPIT
           END IF
         ENDIF
       ENDIF
       NUM_ATMS=NATOMS
C       CALL TRACER('MOLDENDRV NATOMS',NATOMS)
C       CALL TRACER(OPTION,NATOMS)
       ALLOCATE(IDINDX(NATOMS),REQV(3,NATOMS),RNUC(3,NATOMS),
     &          STAT=ALLOCATION)
       IF (ALLOCATION.GT.0) THEN
         WRITE(*,*) 'ERROR: ALLOCATION FAILED IN ','MOLDENDRV'
       END IF
C
C CHECK FOR EXISTENCE OF MOLDEN FILE
       CALL CHECK_MOLDEN
C
C
C      --- WRITE ONLY COORDINATES INTO THE MOLDEN FILE ---
C
       IF (OPTION.EQ."XYZ") THEN
         OPEN(99,FILE='CLUSTER.MOLDEN',FORM='FORMATTED',STATUS='OLD')
         OPEN(98,FILE='TEMP.MOLDEN',FORM='FORMATTED',STATUS='NEW')
C
         ITEST = 0
         DO WHILE(ITEST.EQ.0)
            READ(99,'(132a)',END=100) INLINE
            DO J=1,20
               IF(INLINE(J:J+7).EQ.'[Atoms]') THEN
                 ITEST=1
                 EXIT
               EN DIF
            END DO
C
            DO J=132,1,-1
              IF(INLINE(J:J).ne.' ') EXIT
            END DO
            WRITE(98,'(A)') INLINE(1:J)
         END DO
C
C      --- IF FOUND ATOMS SECTION, REPLACE IT WITH NEW SET ---
C      --- FIRST GET FULL UNSYMMETRIZED CARTESIANS COORDINATES ---
C
         DO INUC=1,NUM_ATMS
           CALL GET_LETTER(XMOL_LIST(INUC)%ANUM,LETTER) 
           WRITE(98,'(1X,A4,2I7,3(1X,F16.6))') LETTER,INUC,
     &           XMOL_LIST(INUC)%ANUM,XMOL_LIST(INUC)%RX*AU2ANG,
     &           XMOL_LIST(INUC)%RY*AU2ANG,XMOL_LIST(INUC)%RZ*AU2ANG
         END DO
C
C      --- SKIP OVER OLD SECTION ---
C
         ITEST = 0
         DO WHILE(ITEST.EQ.0)
           READ(99,'(132a)',end=200) INLINE
           DO J=1,20
             IF(INLINE(J:J).EQ.'[') THEN
               ITEST=1
               EXIT
             END IF
           END DO
         END DO
C
C      --- NOW SPOOL THE REST OF THE FILE OVER ---
C
         DO J=132,1,-1
           IF(INLINE(J:J).ne.' ') EXIT
         END DO
         WRITE(98,'(A)') INLINE(1:J)
         ITEST = 0
         DO WHILE(ITEST.EQ.0)
           READ(99,'(132a)',END=200) INLINE
           DO J=132,1,-1
             IF(INLINE(J:J).NE.' ') EXIT
           END DO
           WRITE(98,'(A)') INLINE(1:J)
         END DO
C
  100  CONTINUE
C
C      --- NO ATOMS SECTION FOUND, WE WILL APPEND ONE ---
C
         WRITE(98,'(a)') '[Atoms] Angs'
C
C      --- WRITE FULL SET OF CARTESIAN COORDIANTES ---
C
         DO INUC=1,NUM_ATMS
           CALL GET_LETTER(XMOL_LIST(INUC)%ANUM,LETTER) 
           WRITE(98,'(1X,A4,2I7,3(1X,F16.6))') LETTER,INUC,
     &           XMOL_LIST(INUC)%ANUM,XMOL_LIST(INUC)%RX*AU2ANG,
     &           XMOL_LIST(INUC)%RY*AU2ANG,XMOL_LIST(INUC)%RZ*AU2ANG
         END DO
C
  200    CONTINUE
c
         CLOSE(99)
         CLOSE(98)
         CALL SYSTEM('rm CLUSTER.MOLDEN')
         CALL SYSTEM('mv TEMP.MOLDEN CLUSTER.MOLDEN')
C
C      --- WRITE COORDINATES, ENERGIES, AND FORCES INTO MOLDEN FILE ---
C
       ELSE IF (OPTION.EQ.'OPT') THEN
C
C      --- SET UP FLAGS FOR WRITING GEOMETRIES, CONVERGENCE INFO, ---
C      --- AND CARTESIAN FORCES ---
C 
         IGEO = 0
         IGCNV = 0
         IFOR = 0
C
         OPEN(99,FILE='CLUSTER.MOLDEN',FORM='FORMATTED',
     &        STATUS='UNKNOWN')
         OPEN(98,FILE='TEMP.MOLDEN',FORM='FORMATTED',STATUS='UNKNOWN')
         ENDFILE(98)
         REWIND(98)
C
C      --- READ FIRST LINE OF EXISTING MOLDEN FILE AND VERIFY THE ---
C      --- CORRECT HEADER, IF NOT, GENERATE AND PASS OVER THE      ---
C      --- CORRECT ONE FOR LATER ---
C
         READ(99,'(132a)',END=900) INLINE
         ITEST = 0
         DO J=1,20
           IF(INLINE(J:J+14).EQ.'[Molden Format]') THEN
             ITEST = 1
             EXIT
           END IF
         END DO
C
         IF(ITEST.NE.1) THEN
           INLINE(1:15)='[Molden Format]'
           J = 1
         END IF
         DO K=132,J,-1
           IF(INLINE(K:K).NE.' ') EXIT
         END DO
         WRITE(98,'(A)') INLINE(J:K)
C
C      --- SCAN TO THE NEXT SET OF MAJOR HEADERS: LINES BEGIN ---
C      --- WITH '[' AND WRITE OVER THE NEW DATA AS NEEDED     ---
C
  350    CONTINUE
C
C      --- FINISH READING IF ALL TASK HAVE BEEN COMPLETED ---
C
         IF ((IGEO.NE.0).AND.(IGCNV.NE.0).AND.(IFOR.NE.0)) THEN
           GO TO 900
         END IF 
C
         ITEST = 0
         DO WHILE(ITEST.EQ.0)
           READ(99,'(132a)',END=900) INLINE
           DO J=1,20
             IF(INLINE(J:J).EQ.'[') THEN
               ITEST=1
               EXIT
             END IF
           END DO
C
           DO K=132,J,-1
             IF(INLINE(K:K).NE.' ') EXIT
           END DO
           WRITE(98,'(A)') INLINE(1:K)
         END DO
         BACKSPACE(98)
         WRITE(98,'(A)') INLINE(J:K)
C
C      --- FOUND THE GEOMETRIES SECTION, NOW WE MOVE TO THE END ---
C      --- OF THE SECTION OF THE END OF THE FILE ---
C
         IF(INLINE(J:J+11).EQ.'[GEOMETRIES]') THEN
C
           IGEO  = 1
           ITEST = 0
           DO WHILE(ITEST.EQ.0)
             READ(99,'(132a)',END=400) INLINE
             DO J=1,20
               IF(INLINE(J:J).EQ.'[') THEN
                 ITEST=1
                 EXIT
               END IF
             END DO
C
             DO K=132,J,-1
               IF(INLINE(K:K).NE.' ') EXIT
             END DO
             WRITE(98,'(A)') INLINE(1:K)
           END DO
C
           BACKSPACE(98)
           BACKSPACE(99)
  400      CONTINUE
C
C      --- NOW APPEND THE NEW DATA ---
C      --- WRITE FULL SET OF CARTESIAN COORDIANTES ---
C
           WRITE(98,*) NUM_ATMS
           WRITE(98,*)

           DO INUC=1,NUM_ATMS
             CALL GET_LETTER(XMOL_LIST(INUC)%ANUM,LETTER) 
             WRITE(98,'(2X,A2,3(2X,F15.7))') LETTER,
     &           XMOL_LIST(INUC)%RX*AU2ANG,
     &           XMOL_LIST(INUC)%RY*AU2ANG,XMOL_LIST(INUC)%RZ*AU2ANG
         END DO
C
           GOTO 350
C
C      --- FOUND GEOCONV SECTION, NOW WE MOVE TO THE END OF THE ---
C      --- ENERGY SUBSECTION BY FINDING THE MAZ-FORCE LINE AND ---
C      --- BACK UP ONE ---
C
         ELSEIF(INLINE(J:J+8).EQ.'[GEOCONV]') THEN
C
           IGCNV = 1
           ITEST = 0
           DO WHILE(ITEST.EQ.0)
             READ(99,'(132a)',end=500) INLINE
             DO J=1,20
               IF(INLINE(J:J+8).EQ.'max-force') THEN
                 ITEST=1
                 EXIT
               END IF
             END DO
C
             DO K=132,1,-1
               IF(INLINE(K:K).ne.' ') EXIT
             END DO
             WRITE(98,'(A)') INLINE(1:K)
           END DO
C
           BACKSPACE(98)
           BACKSPACE(99)
  500      WRITE(98,'(5x,e15.7)') ETOTAL
C
C     --- NOW WE SKIP TO RMS-FORCE SECTION, AND BACK UP ONE ---
C
           ITEST = 0
           DO WHILE(ITEST.EQ.0)
             READ(99,'(132a)',END=600) INLINE
             DO J=1,20
               IF(INLINE(J:J+8).EQ.'rms-force') THEN
                 ITEST=1
                 EXIT
               END IF
             END DO
C
             DO K=132,1,-1
               IF(INLINE(K:K).ne.' ') EXIT
             END DO
             WRITE(98,'(A)') INLINE(1:K)
           END DO
C
C     --- CALCULATE MAX-FORCE ANS WRITE IT INTO MOLDEN FILE ---
C
           GMAX=0.0D0
           GRAD=0.0D0
           DO IID=1,NIDENT
             DO I=1,3
               GMAX=MAX(GMAX,ABS(FTOT(I,IID)))
               GRAD=GRAD+FTOT(I,IID)**2
             END DO
           END DO
           GRAD=SQRT(GRAD)
C
           BACKSPACE(98)
           BACKSPACE(99)
  600      WRITE(98,'(5x,e15.7)') GMAX
C
C     --- FINALLY, SKIP TO THE NEXT SECTION OR EOF ---
C
           ITEST = 0
           DO WHILE(ITEST.EQ.0)
             READ(99,'(132a)',end=700) INLINE
             DO J=1,20
               IF(INLINE(J:J).EQ.'[') THEN
                 ITEST=1
                 EXIT
               END IF
             END DO
C
             DO K=132,1,-1
               IF(INLINE(K:K).ne.' ') EXIT
             END DO
             WRITE(98,'(A)') INLINE(1:K)
           END DO
C
           BACKSPACE(98)
           BACKSPACE(99)
  700      WRITE(98,'(5x,e15.7)') GRAD
C
           GOTO 350
C
C     --- FOUND FORCES SECTION, NOW WE MOVE TO THE END OF THE ---
C     --- SECTION OF THE THE END OFTHE FILE ---
C
         ELSEIF (INLINE(J:J+7).EQ.'[FORCES]') THEN
C
           IFOR = 1
           ITEST = 0
           DO WHILE(ITEST.EQ.0)
             READ(99,'(132a)',end=800) INLINE
             DO J=1,20
               IF(INLINE(J:J).EQ.'[') THEN
                 ITEST = 1
                 EXIT
               END IF
               IF(INLINE(J:J+4).EQ.'point') THEN
                 IFOR = IFOR + 1
                 EXIT
               END IF
             ENDDO
             DO K=132,1,-1
               IF(INLINE(K:K).ne.' ') EXIT
             END DO
             WRITE(98,'(A)') INLINE(1:K)
           END DO
           BACKSPACE(98)
           BACKSPACE(99)
  800      CONTINUE
C
C     --- APPEND NEW DATA ---
C     --- GET FULL UNSYMMETRIZED CARTESIANS FORCES ---
C
         NNUC=0
         DO IID=1,NIDENT
           CALL GASITES(1,FTOT(1,IID),M_NUC,REQV,MSITES)
           DO INUC=1,M_NUC
             RNUC(1,NNUC+INUC)=REQV(1,INUC)
             RNUC(2,NNUC+INUC)=REQV(2,INUC)
             RNUC(3,NNUC+INUC)=REQV(3,INUC)
             IDINDX(NNUC+INUC)=IID
           END DO
           NNUC=NNUC+M_NUC
         END DO
C
C      --- WRITE ALL CARTESIAN FORCES ---
C
           WRITE(98,*) 'point', IFOR
           WRITE(98,*) NUM_ATMS
           DO J=1,NUM_ATMS
             WRITE(98,'(3(2X,F15.7))') (RNUC(I,J), I=1,3)
           END DO
           IFOR = 1
           GOTO 350
C
C     --- SOME OTHER TYPE OF SECTION FOUND, JUST MOVE TO THE ---
C     --- NEXT SECTION --- 
C
         ELSE
          GOTO 350
C
         END IF
C
  900    CONTINUE
C
C      --- END OF THE OLD FILE, CHECK THAT EVERYTHING WAS WRITTEN ---
C      --- IF NOT, APPEND IT TO THE NEW FILE ---
C 
         IF(IGEO.EQ.0) THEN
C
C      --- NO GEOMETRIES SECTION FOUND, WE'LL APPEND ONE ---
C
           WRITE(98,'(a)') '[GEOMETRIES] XYZ'
C
C      --- WRITE FULL SET OF CARTESIAN COORDIANTES ---
C
           WRITE(98,*) NUM_ATMS
           WRITE(98,*)

           DO INUC=1,NUM_ATMS
             CALL GET_LETTER(XMOL_LIST(INUC)%ANUM,LETTER) 
             WRITE(98,'(2X,A2,3(2X,F15.7))') LETTER,
     &           XMOL_LIST(INUC)%RX*AU2ANG,
     &           XMOL_LIST(INUC)%RY*AU2ANG,XMOL_LIST(INUC)%RZ*AU2ANG
           END DO
C
         ENDIF
C
C      --- NO GEOCONV SECTION FOUND, APPEND ONE ---
C
         IF(IGCNV.EQ.0) THEN
C
C      --- CALCULATE RMS_FORCE AND WRITE IT INTO MOLDEN FILE ---
C
           GMAX=0.0D0
           GRAD=0.0D0
           DO IID=1,NIDENT
             DO I=1,3
               GMAX=MAX(GMAX,ABS(FTOT(I,IID)))
               GRAD=GRAD+FTOT(I,IID)**2
             END DO
           END DO
           GRAD=SQRT(GRAD)
C
           WRITE(98,'(a)') '[GEOCONV]'
           WRITE(98,'(a)') 'energy'
           WRITE(98,'(5x,e15.7)') ETOTAL
           WRITE(98,'(a)') 'max-force'
           WRITE(98,'(5x,e15.7)') GMAX
           WRITE(98,'(a)') 'rms-force'
           WRITE(98,'(5x,e15.7)') GRAD
         ENDIF
C
C     ---- NO FORCES SECTION FOUND, WE'LL APPEND ONE ---
C
         IF(IFOR.EQ.0) THEN
C
           WRITE(98,'(a)') '[FORCES]'
C
C     --- GET FULL UNSYMMETRIZED CARTESIANS FORCES ---
C
         NNUC=0
         DO IID=1,NIDENT
           CALL GASITES(1,FTOT(1,IID),M_NUC,REQV,MSITES)
           DO INUC=1,M_NUC
             RNUC(1,NNUC+INUC)=REQV(1,INUC)
             RNUC(2,NNUC+INUC)=REQV(2,INUC)
             RNUC(3,NNUC+INUC)=REQV(3,INUC)
             IDINDX(NNUC+INUC)=IID
           END DO
           NNUC=NNUC+M_NUC
         END DO
C
C     --- WRITE ALL CARTESIAN FORCES ---
C
           WRITE(98,*) 'point', 1
           WRITE(98,*) NUM_ATMS 
           DO J=1,NUM_ATMS
             WRITE(98,'(3(2X,F15.7))') (RNUC(I,J), I=1,3)
           END DO
         END IF
C
C     --- DONE, CLOSE UP SHOP, AND CONTINUE ---
C
         CLOSE(99)
         CLOSE(98)
         CALL SYSTEM('rm CLUSTER.MOLDEN')
         CALL SYSTEM('mv TEMP.MOLDEN CLUSTER.MOLDEN')
C
C     --- WRITE BASIS SETS AND MOS INTO MOLDEN FILE ---
C
       ELSE IF (OPTION.EQ.'MOS') THEN
         CALL MOLDENMOS
         call wfxdrv
C
       ENDIF
C
C      --- DEALLOCATE LOCAL FIELDS ---
C
      DEALLOCATE(IDINDX,REQV,RNUC,STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR: DEALLOCATION FAILED IN ','MOLDENDRV'
      END IF
C
      RETURN
C
C     ------------------------------------------------------------------
C
      END
