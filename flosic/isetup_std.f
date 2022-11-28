C UTEP Electronic Structure Lab (2020)
      SUBROUTINE ISETUP_STD 
C
C     ORIGINAL VERSION BY MARK R PEDERSON (1996)
C     UPDATED BY ULISES REVELES, JULY 2013.
C     UPDATED BY LUIS BASURTO, JANUARY 2014.
C
C     ------------------------------------------------------------------
C
C      use debug1
      use global_inputs,only : CALCTYPE1,MOLDEN1,SPNPOL1
      use common1, only : PSPSYM
      use common2, only : ATOMMAP,ATOMSPIN,CHNET,OPTTYP,SPNNET
      use common2, only: natoms
      use common3, only : RMAT, NGRP
      use common9, only: old_mode
      !<LA: uncomment
      use inputvars, only: functype, symgrp, geometry, atoms, 
     & atoms_int, spinpseudo
      use FRM,only     : LFRM
!      use inp_module
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IATOMS, IATRD, IERR, IGP, IL, IOK, IREAD, ISET,
     & ITOT, IZE, IZN, JTOT, L, MMATOMS, MSITES, MXCHG, NALP, NATOM,
     & NBASF, NFIND, NNUC, NSETS
       REAL*8 :: SYMBOL , ALPSET, CONSET, ELECTRONS, ERROR, RNUC,
     & TOLER, V, ZTOT
       INTEGER :: IPSCHG
       external IPSCHG
      SAVE

      PARAMETER(MXCHG=56)

      INTEGER INEW,TOTNSETS,J,K
      LOGICAL   EXIST,PSEUDO,FIRST,AFERRO

      CHARACTER*60 NAMFUNCT
      CHARACTER*3  NAMGRP,TABLE(0:MXCHG)
      CHARACTER*6  CHARINF
      CHARACTER*3,ALLOCATABLE ::  NAMPSP(:),SPNDIR(:)
      INTEGER,ALLOCATABLE :: IZNUC(:),IZELC(:),IDXSET(:)
      REAL*8,ALLOCATABLE :: R(:,:)

      DIMENSION ALPSET(MAX_BARE),CONSET(MAX_BARE,MAX_CON,3),NBASF(2,3)
      DIMENSION RNUC(3,MX_GRP),MSITES(1)
      DIMENSION V(3),ELECTRONS(2)

      DATA TABLE/
     & 'PSE',
     & 'HYD','HEL',
     & 'LIT','BER','BOR','CAR','NIT','OXY','FLU','NEO',
     & 'SOD','MAG','ALU','SIL','PHO','SUL','CHL','ARG',
     & 'POT','CAL','SCA','TIT','VAN','CHR','MAN','IRO','COB','NIC',
     & 'COP','ZIN','GAL','GER','ARS','SEL','BRO','KRY',
     & 'RUB','STR','YTR','ZIR','NIO','MOL','TEC','RHU','RHO','PAL',
     & 'SLV','CAD','IND','TIN','ANT','TEL','IOD','XEN',
     & 'CES','BAR'/
C
      DATA TOLER/1.0D-5/
C
C     ------------------------------------------------------------------
C
C     --- INTIALIZATION ---
C
      PSEUDO = .FALSE.
      FIRST = .TRUE.
      TOTNSETS = 0
      AFERRO = .FALSE.
C
C     --- CHECK IF CLUSTER FILE EXIST, IF NOT USE DEFAULT  ---
C
      !<LA: old mode start
      if (old_mode) then
       INQUIRE(FILE='CLUSTER',EXIST=EXIST)
       IF (.NOT. EXIST) GOTO 800
C
C     --- OPEN CLUSTER FILE ---
C
       OPEN(90,FILE='CLUSTER',FORM='FORMATTED',STATUS='OLD')
       REWIND(90)
       READ(90,'(A)', END=200) NAMFUNCT

C     --- COORDINATES GIVEN IN ADDITIONAL FILE ---
C     
       IF(NAMFUNCT(1:1).EQ.'@')THEN
        CLOSE(90)
        CALL EZSTART
        WRITE(6,*)'XMOL FILE PROCESSED INTO CLUSTER FILE, RERUN PROGRAM'
        CALL STOPIT
       END IF
C     --- READ SYMMETRY POINT GROUP AND GET GROUP MATRIX IF DESIRED ---
C
       READ(90,'(A3)',END=200) NAMGRP
       CALL GETGRP(NAMGRP)

C     --- READ GROUP INFORMATION, CHECK GROUP, AND STORE IT IN COMMON /GROUP/ ---
C     --- LB: NOW CALLED COMMON3 IN modules2.f90 ---
C
       CALL FGMAT
C
C     --- NOW DEAL WITH ATOMS ---
C
C LB:     --- AND FRAGMENTS --- 
C
       READ(90,*,END=200) MMATOMS
       IF (MMATOMS .LE. 0) THEN
         PRINT *,'ISETUP: NUMBER OF ATOMS IS ZERO'
         CLOSE(90)
         CALL STOPIT
       ELSE IF (MMATOMS .GT. MAX_IDENT) THEN
         PRINT *,'ISETUP: MAX_IDENT SHOULD BE AT LEAST:',MMATOMS
!         GOTO 900
         CLOSE(90)
         CALL STOPIT
       END IF
      !<LA: new_mode
      else
      !<LA: uncomment
        namgrp = symgrp
        namfunct = functype
        call getgrp(namgrp)
        call fgmat
        mmatoms = natoms
        if (mmatoms .le. 0) then
          write(*,*) 'isetup: number of atoms is zero'
          call stopit
        elseif (mmatoms .gt. max_ident) then
          write(*,*) 'isetup: max_ident should be at least:',mmatoms
          call stopit
        endif
      endif
      !<LA: old_mode end
C
C     --- ALLOCATE LOCAL ARRAYS ---
C
      ALLOCATE(NAMPSP(MMATOMS),STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR ALLOCATING NAMPSP'
      ALLOCATE(SPNDIR(MMATOMS),STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR ALLOCATING SPNDIR'
      ALLOCATE(IZNUC(MMATOMS),STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR ALLOCATING IZNUC'
      ALLOCATE(IZELC(MMATOMS),STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR ALLOCATING IZELC'
      ALLOCATE(R(3,MMATOMS),STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR ALLOCATING R'
      ALLOCATE(IDXSET(MMATOMS),STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR ALLOCATING IDXSET'
C
C     --- LOOP OVER NUMBER OF NON-EQUIVALENT ATOMS ---
C
C     --- GENERATE INPUT FILES                   ---
C     --- FIRST, WRITE HEADER OF SYMBOL ---
C
       OPEN(92,FILE='SYMBOL', FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND(92)
C
C     --- CREATE/OPEN THE MOLDEN OUTPUT FILE AND ---
C     --- MAKE SURE THAT IS EMPTY ! ---
C
      if (old_mode) call check_inputs
      !<LA: this does nothing and is never read in, commenting out
!      IF(MOLDEN1)THEN
!        OPEN(42, FILE='CLUSTER.MOLDEN', STATUS='UNKNOWN')
!        REWIND(42)
!        ENDFILE(42)
!        REWIND(42)
!        WRITE(42,'(A)') '[Molden Format]'
!        CLOSE(42)
!      ENDIF
      SELECT CASE (CALCTYPE1)
        CASE (1)
          WRITE(92,'(A)') 'LBFGS'
        CASE (2)
          WRITE(92,'(A)') 'SCF-ONLY'
        CASE (3)
          WRITE(92,'(A)') 'VERLET'
        CASE (4)
          WRITE(92,'(A)') 'VIBRATIONAL'
        CASE (5)
          WRITE(92,'(A)') 'REDUNDANT'
      END SELECT

      WRITE(92,'(A)') NAMFUNCT
      WRITE(92,'(A)') 'OLDMESH'
      WRITE(92,'(A)') '1    NUMBER OF SYMBOLIC FILES'
      WRITE(92,'(A)') 'ISYMGEN = INPUT'
      WRITE(92,'(I3,A)') MMATOMS+2,'  NUMBER OF SYMBOLS IN LIST'
      WRITE(92,'(I3,A)') MMATOMS,  '  NUMBER OF NUCLEI'
      if (old_mode) natoms = mmatoms
C
C     --- CARTESIAN CONSTRAINS: ALL ATOMS ARE FREE TO MOVE ---
C       
      DO IATOMS=1,MMATOMS 
        WRITE(92,'(A)') '1.0  1 1 1'
      END DO
      WRITE(92,'(A)') '  1  CALCULATION SETUP BY ISETUP'

      PRINT '(A)','GENERATING SYMBOL FILE FROM DATA IN FILE CLUSTER'

      NATOM=0
      NSETS=0
      DO IATRD=1,MMATOMS
        NATOM=NATOM+1
        IDXSET(NATOM)=0
        NSETS=NSETS+1  
C
C     --- READ INFO FOR ONE ATOM, CHECK IF CHARGE IN BOUNDS ---
C
        !<LA: old_mode 
        if (old_mode) then
         CALL IGETATM(90,R(1,NATOM),IZNUC(NSETS),CHARINF,IREAD)
         NAMPSP(NSETS)(1:3)=CHARINF(1:3)
         SPNDIR(NATOM)(1:3)=CHARINF(4:6)
         IF (IREAD .NE. 0) GOTO 200
         IF (NAMPSP(NSETS) .EQ. 'BHS') PSEUDO= .TRUE.
         IF (IZNUC(NSETS) .GT. MXCHG) THEN
          PRINT *,'ISETUP: ATOMS WITH Z > ',MXCHG,' ARE NOT SUPPORTED'
          CLOSE(90)
          CALL STOPIT
         ELSE IF (IZNUC(NSETS) .LT. 0) THEN
          PRINT *,'ISETUP: ATOMS WITH Z < 0 ARE NOT SUPPORTED'
          CLOSE(90)
          CALL STOPIT
         END IF
        else
        !LA: new mode
      !<LA: uncomment
         r(:,natom) = geometry(:,natom)
         charinf = spinpseudo(nsets)
         nampsp(nsets)(1:3) = charinf(1:3)
         spndir(natom)(1:3) = charinf(4:6)
         iznuc(nsets) = atoms_int(nsets)
         if (nampsp(nsets) .eq. 'bhs') pseudo= .true.
         if (iznuc(nsets) .gt. mxchg) then
          print *,'isetup: atoms with z > ',mxchg,' are not supported'
          call stopit
         elseif (iznuc(nsets) .lt. 0) then
          print *,'isetup: atoms with z < 0 are not supported'
          call stopit
         endif
        endif

C
C     --- CHECK FOR EQUIVALENT ATOMS ---
C
        ITOT=0
        DO IATOMS=1,NATOM-1
          JTOT=0
          DO IGP=1,NGRP
            DO I=1,3
              V(I)=0.0D0
              DO J=1,3
                V(I)=V(I)+RMAT(I,J,IGP)*R(J,IATOMS)
              END DO
            END DO
            ERROR=ABS(V(1)-R(1,NATOM)) 
     &           +ABS(V(2)-R(2,NATOM)) 
     &           +ABS(V(3)-R(3,NATOM)) 
            IF (ERROR .LT. TOLER) JTOT=JTOT+1
          END DO
          IF ((JTOT .NE. 0) .AND. 
     &       ((IZNUC (IDXSET(IATOMS)) .NE. IZNUC (NSETS)) .OR.
     &        (NAMPSP(IDXSET(IATOMS)) .NE. NAMPSP(NSETS)))) THEN
            PRINT *,'ISETUP: FOUND ATOMS OF DIFFERENT TYPE AT THE'
            PRINT *,'SAME COORDINATES - ASK A NUCLEAR PHYSICIST'
!            GOTO 900
            if (old_mode) CLOSE(90)
            CALL STOPIT
          END IF
          ITOT=ITOT+JTOT
        END DO
C
C     --- FORGET THIS ATOM IF IT CAN BE CONSTRUCTED BY SYMMETRY ---
C
        IF (ITOT .NE. 0) THEN 
          NATOM=NATOM-1
          NSETS=NSETS-1
        ELSE
C
C     --- CHECK IF NEW ATOM TYPE BELOGS TO AN ALREDY DEFINED SET ---
C
          ITOT=NSETS
          DO ISET=1,NSETS-1 
            IF ((IZNUC(ISET) .EQ. IZNUC(NSETS)) .AND.
     &        (NAMPSP(ISET) .EQ. NAMPSP(NSETS))) ITOT=ISET
          END DO
          IF (ITOT .NE. NSETS) NSETS=NSETS-1
          IDXSET(NATOM)=ITOT
        END IF
C
C     --- END OF LOOP FOR ATOMS IN FRAGMENT ---
C
      END DO
C
C     --- NOW UPDATE SYMBOL AND PSPINP ---
      TOTNSETS=TOTNSETS+NSETS
C
C
C     --- GENERATE POINTER TO MAP COORDINATES:                  ---
C     --- ORIGINAL ORDER IN CLUSTER -> NEW ORDER IN SYMBOL FILE ---
C
      INEW = 0 
      ATOMMAP(1:NATOM,3:4) = 0
C
      DO ISET=1,NSETS  
        DO IATOMS=1,NATOM
          IF (IDXSET(IATOMS) .EQ. ISET) THEN
            INEW = INEW + 1
            ATOMMAP(IATOMS,3) = INEW
          END IF 
        END DO 
      END DO
C
C     --- STORE ATOMMAP ON TAPE ---
C
C YY 8.5.2020 commented out the next two lines
C     CALL OPTIO(ATOMMAP(1:NATOM,3),NATOM,'WRITE','ATOMMAP')
C     CALL OPTIO(ATOMMAP(1:NATOM,3),NATOM,'CLOSE','ATOMMAP')
C
C     --- GENERATE REVERSE POINTER ---
C
      DO IATOMS=1,NATOM
        ATOMMAP(ATOMMAP(IATOMS,3),4) = IATOMS
      END DO
CJUR
C     --- WRITE POINTER - USE FOR DEBUG --
C
C     WRITE(*,*)
C     WRITE(*,*)'POINTER MAPPING REORDERING OF COORDINATES'
C     WRITE(*,*)'ATOMMAP(I,3)'
C     DO IATOM=1,NATOM
C       WRITE(*,*)IATOM,ATOMMAP(IATOM,3)
C     END DO
C
C     WRITE(*,*)
C     WRITE(*,*)'ATOMMAP(I,4)'
C     DO IATOM=1,NATOM
C       WRITE(*,*)IATOM,ATOMMAP(IATOM,4)
C     END DO
CJUR
C     --- LOOP OVER SETS OF ATOMS ---
C
      DO ISET=1,NSETS  
        IZN=IZNUC(ISET)
        PSPSYM(1)(1:3)= NAMPSP(ISET)
        PSPSYM(1)(4:4)= '-'
        PSPSYM(1)(5:7)= TABLE(IZN)
C
C     --- FIND ATOMS THAT BELONG TO THIS SET ---
C
        NFIND=0
        DO IATOMS=1,NATOM
          IF (IDXSET(IATOMS) .EQ. ISET) THEN
            NFIND=NFIND+1
          END IF 
        END DO 
C
C     --- ADD TO SYMBOL ---
C
        NFIND=0
        DO IATOMS=1,NATOM
          IF (IDXSET(IATOMS) .EQ. ISET) THEN
            NFIND=NFIND+1
            IOK=0
            IF(SPNDIR(IATOMS).EQ.'SUP')IOK=1
            IF(SPNDIR(IATOMS).EQ.'SDN')IOK=1
            IF(SPNDIR(IATOMS).EQ.'UPO')IOK=1
            IF(IOK.EQ.0)SPNDIR(IATOMS)='UPO'
C
C     --- SET ATOMSPIN ---
C
            !YY copy spin from CLUSTER to SYMBOL
            if(IOK .eq. 1) ATOMSPIN(IATOMS) = SPNDIR(IATOMS)

            IF ((ATOMSPIN(IATOMS).EQ.'SUP').OR.
     &          (ATOMSPIN(IATOMS).EQ.'SDN')) THEN
              AFERRO = .TRUE.
            ELSE
              ATOMSPIN(IATOMS) = 'UPO' 
            END IF
            WRITE(92,1040) PSPSYM(1),NFIND,(R(J,IATOMS),J=1,3),
     &                     ATOMSPIN(IATOMS) 
 1040       FORMAT(A7,I3.3,' =',3(1X,F15.6),' ',A3)
          END IF
        END DO
      ENDDO
C
      ZTOT=0.0D0
      DO IATOMS=1,NATOM
        ISET=IDXSET(IATOMS)
        IZELC(ISET)=IPSCHG(NAMPSP(ISET),IZNUC(ISET))  
        CALL GASITES(1,R(1,IATOMS),NNUC,RNUC,MSITES)
        ZTOT=ZTOT+NNUC*IZELC(ISET)
      END DO
C
C     --- DEFINE NET CHARGE AND SPIN IPSCHG RETURNS THE ---
C     --- ACTUAL NUMBER OF ELECTRONS ON THE ATOM        ---
C
      if (old_mode) READ(90,*,END=200) CHNET,SPNNET
C
      ELECTRONS(1)= -(CHNET-SPNNET-ZTOT)/2.0D0
      ELECTRONS(2)= -(CHNET+SPNNET-ZTOT)/2.0D0
C
C     --- CLOSE CLUSTER FILE ---
C
      if (old_mode) CLOSE(90)
!YY. spin polarized calculation if SPNPOLV is set to Y
      if (old_mode) call check_inputs
      IF(SPNPOL1)  AFERRO = .true.
!If FRMIDT contains down FODs, make it spin polarized caculation
      IF(LFRM(2) .GT. 0) AFERRO = .true.

C
C     --- ADD SOME STUFF TO SYMBOL ---
C
      IF (AFERRO) THEN
        WRITE(92,1080) ELECTRONS(1),-ELECTRONS(2)
      ELSE
        WRITE(92,1080) ELECTRONS
      END IF
      WRITE(92,1090)

 1080 FORMAT('ELECTRONS  =',2(1X,F15.6))
 1090 FORMAT('EXTRABASIS = 0')
C
C     --- CLOSE SYMBOL FILE ---
C
      CLOSE(92)
C
C     --- OPEN CLUSTER FILE ---
C
!<LA: why is this done twice? the next sections are repeated from
!earlier

      if (old_mode) then
       OPEN(90,FILE='CLUSTER',FORM='FORMATTED',STATUS='OLD')
       REWIND(90)
       READ(90,'(A)', END=200) NAMFUNCT
C
C     --- COORDINATES GIVEN IN ADDITIONAL FILE ---
C     
       IF(NAMFUNCT(1:1).EQ.'@')THEN
        CLOSE(90)
        CALL EZSTART
        WRITE(6,*)'XMOL FILE PROCESSED INTO CLUSTER FILE, RERUN PROGRAM'
        CALL STOPIT
       END IF
C
C     --- READ SYMMETRY POINT GROUP AND GET GROUP MATRIX IF DESIRED ---
C
       READ(90,'(A3)',END=200) NAMGRP
C
C     --- NOW DEAL WITH ATOMS ---
C LB:     --- AND FRAGMENTS FOR ISYMGEN --- 
C
       READ(90,*,END=200) MMATOMS
       IF (MMATOMS .LE. 0) THEN
         PRINT *,'ISETUP: NUMBER OF ATOMS IS ZERO'
!         GOTO 900
         CLOSE(90)
         CALL STOPIT
       ELSE IF (MMATOMS .GT. MAX_IDENT) THEN
         PRINT *,'ISETUP: MAX_IDENT SHOULD BE AT LEAST:',MMATOMS
!         GOTO 900
         CLOSE(90)
         CALL STOPIT
       END IF
      endif
C
C     --- PROCESS ISYMGEN ---
C
       !<LA: The next few lines are unique 
      IF (PSEUDO) THEN
        OPEN(94,FILE='PSPINP' ,FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(94)
      END IF

      OPEN(92,FILE='ISYMGEN',FORM='FORMATTED',STATUS='UNKNOWN')
      REWIND(92)

      PRINT '(A)','GENERATING ISYMGEN FILE FROM DATA IN FILE CLUSTER'
      !<LA: END of unique lines

      NSETS=0
      NATOM=0
      DO IATRD=1,MMATOMS
        NATOM=NATOM+1
        IDXSET(NATOM)=0
        NSETS=NSETS+1  
C
C     --- READ INFO FOR ONE ATOM, CHECK IF CHARGE IN BOUNDS ---
C
        if (old_mode) then
         CALL IGETATM(90,R(1,NATOM),IZNUC(NSETS),CHARINF,IREAD)
         NAMPSP(NSETS)(1:3)=CHARINF(1:3)
         SPNDIR(NATOM)(1:3)=CHARINF(4:6)
         IF (IREAD .NE. 0) GOTO 200
         IF (NAMPSP(NSETS) .EQ. 'BHS') PSEUDO= .TRUE.
         IF (IZNUC(NSETS) .GT. MXCHG) THEN
          PRINT *,'ISETUP: ATOMS WITH Z > ',MXCHG,' ARE NOT SUPPORTED'
!          GOTO 900
          CLOSE(90)
          CALL STOPIT
         ELSE IF (IZNUC(NSETS) .LT. 0) THEN
          PRINT *,'ISETUP: ATOMS WITH Z < 0 ARE NOT SUPPORTED'
!          GOTO 900
          CLOSE(90)
          CALL STOPIT
         END IF
        else
        !LA: new mode
      !<LA: uncomment
         r(:,natom) = geometry(:,natom)
         charinf = spinpseudo(nsets)
         nampsp(nsets)(1:3) = charinf(1:3)
         spndir(natom)(1:3) = charinf(4:6)
         iznuc(nsets) = atoms_int(nsets)
         if (nampsp(nsets) .eq. 'bhs') pseudo= .true.
         if (iznuc(nsets) .gt. mxchg) then
          print *,'isetup: atoms with z > ',mxchg,' are not supported'
          call stopit
         elseif (iznuc(nsets) .lt. 0) then
          print *,'isetup: atoms with z < 0 are not supported'
          call stopit
         endif
        endif
C
C     --- CHECK FOR EQUIVALENT ATOMS WTHIN THIS FRAGMENT ---
C
        ITOT=0
        DO IATOMS=1,NATOM-1
          JTOT=0
          DO IGP=1,NGRP
            DO I=1,3
              V(I)=0.0D0
              DO J=1,3
                V(I)=V(I)+RMAT(I,J,IGP)*R(J,IATOMS)
              END DO
            END DO

            ERROR=ABS(V(1)-R(1,NATOM)) 
     &           +ABS(V(2)-R(2,NATOM)) 
     &           +ABS(V(3)-R(3,NATOM)) 
            IF (ERROR .LT. TOLER) JTOT=JTOT+1
          END DO
          IF ((JTOT .NE. 0) .AND. 
     &       ((IZNUC (IDXSET(IATOMS)) .NE. IZNUC (NSETS)) .OR.
     &        (NAMPSP(IDXSET(IATOMS)) .NE. NAMPSP(NSETS)))) THEN
            PRINT *,'ISETUP: FOUND ATOMS OF DIFFERENT TYPE AT THE'
            PRINT *,'SAME COORDINATES - ASK A NUCLEAR PHYSICIST'
!            GOTO 900
            if (old_mode) CLOSE(90)
            CALL STOPIT
          END IF
          ITOT=ITOT+JTOT
        END DO
C
C     --- FORGET THIS ATOM IF IT CAN BE CONSTRUCTED BY SYMMETRY ---
C
        IF (ITOT .NE. 0) THEN 
          NATOM=NATOM-1
          NSETS=NSETS-1
        ELSE
C
C     --- CHECK IF NEW ATOM TYPE BELOGS TO AN ALREDY DEFINED SET ---
C
          ITOT=NSETS
          DO ISET=1,NSETS-1 
            IF ((IZNUC (ISET) .EQ. IZNUC (NSETS)) .AND.
     &        (NAMPSP(ISET) .EQ. NAMPSP(NSETS))) ITOT=ISET
          END DO
          IF (ITOT .NE. NSETS) NSETS=NSETS-1
          IDXSET(NATOM)=ITOT
        END IF
C
C END OF LOOP OVER ATOMS IN FRAGMENT
C
      END DO
      if (old_mode) CLOSE(90)
      !<LA: end of first repeated section 

C     --- WRITE INTO ISYMGEN FILE ---
C
      IF(FIRST) THEN
        WRITE(92,'(1X,I3,10X,A)') TOTNSETS,'TOTAL NUMBER OF ATOM TYPES'
        FIRST=.FALSE.
      ENDIF
      !<LA: the following lines are repeated but not in the same order
      ! as earlier, seems to have to do with writing to isymgen, so the
      ! question is why did the prior section need to be repeated 
C
C     --- GENERATE POINTER TO MAP COORDINATES:                  ---
C     --- ORIGINAL ORDER IN CLUSTER -> NEW ORDER IN SYMBOL FILE ---
C
      !<LA: atommap is reset but not called again
      INEW = 0 
      ATOMMAP(1:NATOM,3:4) = 0
C
C     --- LOOP OVER SETS OF ATOMS ---
C
      DO ISET=1,NSETS  
        IZN=IZNUC(ISET)
        IZE=IZELC(ISET)
        PSPSYM(1)(1:3)= NAMPSP(ISET)
        PSPSYM(1)(4:4)= '-'
        PSPSYM(1)(5:7)= TABLE(IZN)

        CALL SETPSP(94,PSPSYM(1),IZN)
        CALL SETBAS(IZN,NAMPSP(ISET),ALPSET,CONSET,NALP,NBASF)
C
C     --- FIND ATOMS THAT BELONG TO THIS SET ---
C
        NFIND=0
        DO IATOMS=1,NATOM
          IF (IDXSET(IATOMS) .EQ. ISET) THEN
            NFIND=NFIND+1
          END IF 
        END DO 
C
C WRITE TO ISYMGEN
C
        PRINT 1020,ISET,PSPSYM(1)(1:7),IZN,NFIND
 1020   FORMAT('ATOM TYPE ',I3,' (',A,'): NUCLEAR CHARGE= ',I3,', ',
     &         I3,' ATOM(S)')
        WRITE(92,'(2(1X,I3),6X,A)') IZE,IZN,
     &           'ELECTRONIC AND NUCLEAR CHARGE'
        IF (NAMPSP(ISET) .EQ. 'ALL') THEN
          WRITE(92,'(A3,11X,A)') PSPSYM(1),'ALL-ELECTRON ATOM TYPE'
        ELSE
          WRITE(92,'(A7,7X,A)') PSPSYM(1),'PSP-SYMBOL'
        END IF
        WRITE(92,'(1X,I3,10X,A,A3)') NFIND,
     &           'NUMBER OF ATOMS OF TYPE ',TABLE(IZN)
C
C     --- ADD TO ISYMGEN ---
C
        NFIND=0
        DO IATOMS=1,NATOM
          IF (IDXSET(IATOMS) .EQ. ISET) THEN
            NFIND=NFIND+1
            WRITE(92,1030) PSPSYM(1),NFIND
 1030       FORMAT(A7,I3.3)
          END IF
        END DO
C
C     --- ADD BASIS SET TO ISYMGEN ---
C 
        WRITE(92,'(2A)') 'EXTRABASIS    CONTROLS USAGE OF ',
     &                   'SUPPLEMENTARY BASIS FUNCTIONS'
        WRITE(92,1050) NALP,'NUMBER OF BARE GAUSSIANS'
 1050   FORMAT(1X,I3,10X,A) 
        WRITE(92,1060)(NBASF(1,L), L=1,3),'NUMBER OF S,P,D FUNCTIONS'
        WRITE(92,1060)(NBASF(2,L), L=1,3),
     &                'SUPPLEMENTARY S,P,D FUNCTIONS'
 1060   FORMAT(3(1X,I3),2X,A) 
 1070 FORMAT(3(1X,D22.10))  !YY Increased precision by 2 digits
        WRITE(92,1070)(ALPSET(J), J=1,NALP)
        WRITE(92,*)
C
        if (old_mode) call check_inputs
        !write L=1:3 S,P,D types
        DO L=1,3
          !write functions plus supplementary functions         
          DO IL=1,NBASF(1,L)+NBASF(2,L)
            !write coefficients NALP times
            !YY If molden1 is on, enforce normalization of basis set
            !IF(MOLDEN1) call renorm(L,ALPSET(:),CONSET(:,IL,L),NALP)
            WRITE(92,1070)(CONSET(J,IL,L), J=1,NALP)
            WRITE(92,*)
          END DO
        END DO 
      ENDDO

      WRITE(92,'(A)') 'ELECTRONS'
      WRITE(92,'(A)') 'WFOUT'
C
C     --- CLOSE ISYMGEN FILE ---
C
      CLOSE(92)
      IF (PSEUDO) CLOSE(94)

      PRINT '(A)',' '
C
C     --- DEALLOCATE LOCAL ARRAYS ---
C
      DEALLOCATE(NAMPSP,STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR DEALLOCATING NAMPSP'
      DEALLOCATE(SPNDIR,STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR DEALLOCATING SPNDIR'
      DEALLOCATE(IZNUC,STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR DEALLOCATING IZNUC'
      DEALLOCATE(IZELC,STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR DEALLOCATING IZELC'
      DEALLOCATE(R,STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR DEALLOCATING R'
      DEALLOCATE(IDXSET,STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ISETUP: ERROR DEALLOCATING IDXSET'

      RETURN
C
C     --- SETUP DEFAULT CLUSTER FILE ---
C
  800 OPEN(90,FILE='CLUSTER',FORM='FORMATTED',STATUS='NEW')
      REWIND(90)
      WRITE(90,'(2A)') 'GGA-PBE*GGA-PBE          ',
     &                 '(DF TYPE EXCHANGE*CORRELATION)'
      WRITE(90,'(2A)') 'TD                       ',
     &                 '(TD, OH, IH, X, Y, XY, ... OR GRP)'
      WRITE(90,'(2A)') '2                        ',
     &                 '(NUMBER OF INEQUIV. ATOMS IN CH4)'
      WRITE(90,'(2A)') '0.00  0.00  0.00  6  BHS ',
     &                 '(R, Z, PSEUDOPOT. TYPE FOR CARBON)'
C      WRITE(90,'(2A)') '1.189042  1.189042  1.189042  1  ALL ',
C     &                 '(R, Z, PSEUDOPOT. TYPE FOR HYDROGEN)'
      WRITE(90,'(2A)') '1.20  1.20  1.20  1  ALL ',
     &                 '(R, Z, PSEUDOPOT. TYPE FOR HYDROGEN)'
      WRITE(90,'(2A)') '0.0 0.0                  ',
     &                 '(NET CHARGE AND NET SPIN)' 
      WRITE(90,*)'--------------OR-------------------'
      WRITE(90,'(1A)')'@XMOL.DAT'
      WRITE(90,*)'IF YOU WISH TO START FROM AN XYZ XMOL FILE'
      CLOSE(90)
      PRINT '(A)','LOOK AT AND EDIT FILE NAMED CLUSTER, THEN RERUN'
      RETURN
C
C     --- ERROR ---
C
 200  PRINT *,'ISETUP: FILE CLUSTER IS INVALID'
       CLOSE(90)
       CALL STOPIT
C
C     --- ERROR ---
C
!  900 CLOSE(90)
!      CALL STOPIT
C
C     ------------------------------------------------------------------
C
      END
