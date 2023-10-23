C UTEP Electronic Structure Lab (2020)
       SUBROUTINE DVPMESH(ITBEG,ICALC)

       use xmol,only : NUM_ATMS,XMOL_LIST,AU2ANG,GET_LETTER
       use mesh1,only : wmsh,rmsh,nmsh
       use common2,only : RIDT, IFUIDT, NIDENT, ZNUC,
     &   N_CON, LSYMMAX, NFNCT
       use common3,only : RMAT, NGRP
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:48 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: ITBEG, ICALC, I, IANG, IC, ICON, ICONA, ICONB, IERR,
     & IFNCT, IGRP, IID, INATM, INDX, INUC, IPV, IRAD, IS, ISA, ISB, J,
     & JANG, JCALC, JFNCT, JNUC, L, L1, L1A, L1B, LA, LB, LCUR,
     & LMAXSTO, LMXANG, LMXRAD, LS, LSA, LSB, LSIZ, MAX_CNT, MAX_GRP,
     & MAXANG, MAXRAD, MODPART, MSITES, NANG, NDECAY, NLW, NNUC, NPTAT,
     & NPV, NR1, NR2, NRAD, NSTO, NTOT, NUP
       REAL*8 :: SYMBOL , ANNFAC, BF, BF2, DIFF, ERRAD, ERRANG, ERRCOR,
     & ERROR, FAC, POTF, POTMX, RA, RASMAX, RCOR, RDIF, RDN, RHOF,
     & RNNFAC, RNNS, RRC, RRSMAX, RTMP, RTRAFO, RTRINIT, RUP, SREF,
     & SUM, THRD, TIME1, TIME2, TIME3, TIME4, TIMEA, TIMEB, TOLER
        SAVE
C
C MXRAD:  MAX. NUMBER OF POINTS FOR A RADIAL MESH
C MXANG:  MAX. NUMBER OF ANGULAR POINTS FOR AN ATOMIC SPHERE
C
        PARAMETER (MAXRAD=300)
        PARAMETER (MAXANG=1250)
        PARAMETER (LSIZ=(LDIM*(LDIM+1))/2)
        LOGICAL EXIST,LSYWT,ISUCC
C
C LOCAL STUFF 
C
        LOGICAL IUPDAT
C        COMMON/TMP1/STOMSH(6,MAXANG,2),PTS(NSPEED,3),RVEC(3)
C     &   ,GRAD(NSPEED,10,LSIZ,MAX_CON,LDIM)
C     &   ,OVLTAB(MAX_CON,LDIM,MAX_FUSET)
C     &   ,OVLDIA(2,LSIZ,MAX_CON,LDIM,MAX_FUSET)
C     &   ,ONSOVL(2,LSIZ,LSIZ,MAX_CON,LDIM,MAX_CON,LDIM)
C     &   ,OVLSUM(2,LSIZ,MAX_CON,LDIM,MX_CNT)
C     &   ,RRAD(MAXRAD),RANG(3,MAXANG),WRAD(MAXRAD),WANG(MAXANG)
C     &   ,RNUC(3,MX_CNT),REQV(3,MX_GRP),RSQR(MX_CNT)
C     &   ,RCDTAB(MX_CNT,MX_CNT)
C     &   ,NMULT(MAX_IDENT),IDINDX(MX_CNT),IFINDX(MX_CNT)
C     &   ,ICOUNT(MAX_CON,LDIM)
        REAL*8 :: RVEC(3)
        REAL*8,ALLOCATABLE :: STOMSH(:,:,:),PTS(:,:),GRAD(:,:,:,:,:),
     &   OVLTAB(:,:,:),OVLDIA(:,:,:,:,:),ONSOVL(:,:,:,:,:,:,:),
     &   OVLSUM(:,:,:,:,:),RRAD(:),RANG(:,:),WRAD(:),WANG(:),RNUC(:,:),
     &   REQV(:,:),RSQR(:),RCDTAB(:,:)
        INTEGER,ALLOCATABLE :: NMULT(:),IDINDX(:),IFINDX(:)
        LOGICAL,ALLOCATABLE :: ICOUNT(:,:)
        CHARACTER*2         :: LETTER
C
        DATA LSYWT /.FALSE./
        DATA MODPART,NDECAY /1, 2/
        DATA ERRCOR /1.0D-6/
        DATA LMXRAD,ERRAD,LMXANG,ERRANG /17, 1.0D-6, 29, 1.0D-6/
        DATA RNNFAC,ANNFAC /2.0D0, 2.0D0/
C
C TOLER IS A GENERAL ACCURACY MARGIN.
C
        DATA TOLER /1.0D-8/
C
C ALLOCATE LOCAL ARRAYS
C
       CALL GTTIME(TIMEA)
       ALLOCATE(STOMSH(6,MAXANG,2),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING STOMSH'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING PTS'
       ALLOCATE(GRAD(NSPEED,10,LSIZ,MAX_CON,LDIM),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING STOMSH'
       ALLOCATE(OVLTAB(MAX_CON,LDIM,MAX_FUSET),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING OVLTAB'
       ALLOCATE(OVLDIA(2,LSIZ,MAX_CON,LDIM,MAX_FUSET),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING OVLDIA'
       ALLOCATE(ONSOVL(2,LSIZ,LSIZ,MAX_CON,LDIM,MAX_CON,LDIM),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING ONSOVL' 
       ALLOCATE(OVLSUM(2,LSIZ,MAX_CON,LDIM,MX_CNT),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING OVLSUM' 
       ALLOCATE(RRAD(MAXRAD),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING RRAD' 
       ALLOCATE(RANG(3,MAXANG),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING RANG' 
       ALLOCATE(WRAD(MAXRAD),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING WRAD' 
       ALLOCATE(WANG(MAXANG),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING WANG' 
       ALLOCATE(RNUC(3,MAX_CNT),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING RNUC' 
       ALLOCATE(REQV(3,MAX_GRP),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING REQV' 
       ALLOCATE(RSQR(MX_CNT),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING RSQR' 
       ALLOCATE(RCDTAB(MX_CNT,MX_CNT),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING RCDTAB' 
       ALLOCATE(NMULT(MAX_IDENT),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING NMULT' 
       ALLOCATE(IDINDX(MX_CNT),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING IDINDX' 
       ALLOCATE(IFINDX(MX_CNT),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING IFINDX' 
       ALLOCATE(ICOUNT(MAX_CON,LDIM),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR ALLOCATING ICOUNT' 
       CALL GTTIME(TIMEB)
       WRITE(6,*)'DVPMESH: ALLOCATION TIME',TIMEB-TIMEA
C
C CHECK IF OLDMESH EXISTS AND READ OLD MESH IF DESIRED
C

        PRINT '(A)',' '
        PRINT '(A)','CREATING INTEGRATION MESH'
        CALL GTTIME(TIME1)
        TIME2=TIME1
        TIME3=TIME1
        IF (ITBEG.NE.0) THEN
         IF (ITBEG.LT.0) ITBEG=0
         INQUIRE(FILE='OLDMESH',EXIST=EXIST)
         IF (EXIST) THEN
          PRINT '(A)','TRYING TO READ OLD MESH'
          OPEN(99,FILE='OLDMESH',FORM='UNFORMATTED',STATUS='UNKNOWN')
          REWIND(99)
          READ(99) NMSH,JCALC
          PRINT '(A,I10)','TOTAL NUMBER OF MESH POINTS: ',NMSH
          IF (NMSH.GT.MAX_PTS) THEN
           write(6,*)'DVPMESH: MAX_PTS MUST BE AT LEAST: ',NMSH
           CALL STOPIT
          END IF
          READ(99)((RMSH(J,I),J=1,3),I=1,NMSH)
          READ(99)(WMSH(I),I=1,NMSH)
          CLOSE(99)
          IF (JCALC.EQ.ICALC) THEN
           CALL GTTIME(TIME2)
           CALL REORMSH
           CALL GTTIME(TIME3)
           GOTO 900
          END IF
          write(6,*)'DVPMESH: OLDMESH IS INCOMPATIBLE WITH CURRENT'
          write(6,*)'         GEOMETRY --> A NEW MESH WILL BE CREATED'
          write(6,*)'         OLDMESH WILL BE SAVED IN FILE SAVMESH'
          OPEN(99,FILE='SAVMESH',FORM='UNFORMATTED',STATUS='UNKNOWN')
          WRITE(99) NMSH,JCALC
          WRITE(99)((RMSH(J,I),J=1,3),I=1,NMSH)
          WRITE(99)(WMSH(I),I=1,NMSH)
          CLOSE(99)
         END IF
        END IF
C
C GET DEGREE OF HIGHEST LEBEDEV MESH
C
        LMAXSTO= -1
        CALL LEBEDEV(MAXANG,LMAXSTO,NANG,RANG,WANG)
C
C READ DVPMESH
C
        OPEN(99,FILE='DVPMESH',FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(99)
        READ(99,*,END=10) LSYWT
        READ(99,*,END=10) MODPART,NDECAY
        READ(99,*,END=10) ERRCOR
        READ(99,*,END=10) LMXRAD,ERRAD,RNNFAC
        READ(99,*,END=10) LMXANG,ERRANG,ANNFAC
        IF (NDECAY.LT.1) THEN
         write(6,*)'DVPMESH: NDECAY MUST BE >= 1'
         CALL STOPIT
        END IF
        IF (ERRCOR .LE. 0.0D0) THEN
         write(6,*)'DVPMESH: ERRCOR MUST BE > 0'
         CALL STOPIT
        END IF
        IF ((LMXRAD .LE. 0) .OR. (LMXANG .LE. 0)) THEN
         write(6,*)'DVPMESH: LMXRAD AND LMXANG MUST BE >= 0'
         CALL STOPIT
        END IF
        IF ((LMXRAD .GT. LMAXSTO) .OR. (LMXANG .GT. LMAXSTO)) THEN
         write(6,*)'DVPMESH: LMXRAD AND LMXANG MUST BE <= ',LMAXSTO
         CALL STOPIT
        END IF
        IF ((ERRAD .LE. 0.0D0) .OR. (ERRANG .LE. 0.0D0)) THEN
         write(6,*)'DVPMESH: ERRAD AND ERRANG MUST BE > 0'
         CALL STOPIT
        END IF
        IF ((RNNFAC .LT. 0.0D0) .OR. (ANNFAC .LT. 0.0D0)) THEN
         write(6,*)'DVPMESH: RNNFAC AND ANNFAC MUST BE >= 0'
         CALL STOPIT
        END IF
C
C WRITE DVPMESH
C
   10   REWIND(99)
        WRITE(99,'(1X,L5)') LSYWT
        WRITE(99,'(1X,I5,1X,I5)') MODPART,NDECAY
        WRITE(99,'(1X,E12.5)') ERRCOR 
        WRITE(99,'(1X,I5,1X,E12.5,1X,F10.5)') LMXRAD,ERRAD,RNNFAC
        WRITE(99,'(1X,I5,1X,E12.5,1X,F10.5)') LMXANG,ERRANG,ANNFAC
        WRITE(99,*)'------------------------------------------------'
        WRITE(99,*)'* LSYWT'
        WRITE(99,*)'* MODPART (1 -> RHO*POT/R**2, 2 -> BECKE), NDECAY'
        WRITE(99,*)'* ERRCOR'
        WRITE(99,*)'* LMXRAD, ERRAD,  RNNFAC'
        WRITE(99,*)'* LMXANG, ERRANG, ANNFAC'
        CLOSE(99)
C
C FOR EACH IDENTITY MEMBER, DETERMINE:
C * MULTIPLICITY
C * COORDINATE
C * INDEX OF IDENTITY MEMBER
C
        NNUC=0
        DO IID=1,NIDENT
         CALL GASITES(1,RIDT(1,IID),NMULT(IID),REQV,MSITES)
         NTOT=NNUC+NMULT(IID)
         IF (NTOT .GT. MX_CNT) THEN
          write(6,*)'DVPMESH: MX_CNT MUST BE AT LEAST: ',NTOT
          CALL STOPIT
         END IF
         DO INUC=1,NMULT(IID)
          RNUC(1,NNUC+INUC)=REQV(1,INUC)
          RNUC(2,NNUC+INUC)=REQV(2,INUC)
          RNUC(3,NNUC+INUC)=REQV(3,INUC)
          IDINDX(NNUC+INUC)=IID
         END DO
         NNUC=NNUC+NMULT(IID)
        END DO
C
C WRITE XYZ STRUCTURE FILE WHICH IS XMOL READABLE
C
        OPEN(99,FILE='XMOL.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(99)
        WRITE(99,'(I5)') NNUC
        WRITE(99,'(A)') 'Cluster output'
        NUM_ATMS=NNUC
        ALLOCATE(XMOL_LIST(NNUC),STAT=IERR)
        IF(IERR.NE.0)WRITE(6,*)'DVPMESH:ERROR ALLOCATING XMOL_LIST'
        DO INUC=1,NNUC 
         WRITE(99,'(1X,I4,3(1X,F10.5))') 
     &    NINT(ZNUC(IFUIDT(IDINDX(INUC)))),(AU2ANG*RNUC(I,INUC), I=1,3)
          XMOL_LIST(INUC)%ANUM=NINT(ZNUC(IFUIDT(IDINDX(INUC))))
          XMOL_LIST(INUC)%RX=REQV(1,INUC)
          XMOL_LIST(INUC)%Ry=REQV(2,INUC)
          XMOL_LIST(INUC)%Rz=REQV(3,INUC)
        END DO
        CLOSE(99)
C
C WRITE XYZ STRUCTURE FILE WHICH IS JMOL AND MOLEKEL READABLE
C
        OPEN(99,FILE='XMOL.xyz',FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(99)
        WRITE(99,'(I5)') NNUC
        WRITE(99,'(A)') 'Cluster output'
        DO INUC=1,NNUC 
          CALL GET_LETTER(NINT(ZNUC(IFUIDT(IDINDX(INUC)))),LETTER)
          WRITE(99,'(1X,A,3(1X,F10.5))')
     &      LETTER,(AU2ANG*RNUC(I,INUC), I=1,3)
        END DO
        CLOSE(99)
C
C SETUP RCDTAB
C
        DO INUC=1,NNUC
         RCDTAB(INUC,INUC)= 0.0D0
         DO JNUC=INUC+1,NNUC
          RDIF= (RNUC(1,INUC)-RNUC(1,JNUC))**2
     &         +(RNUC(2,INUC)-RNUC(2,JNUC))**2
     &         +(RNUC(3,INUC)-RNUC(3,JNUC))**2
          IF (RDIF .LT. TOLER) THEN
           write(6,*)'DVPMESH: FOUND TWO ATOMS AT IDENTICAL POSITIONS'
           CALL STOPIT
          END IF
          RDIF= 1.0D0/SQRT(RDIF)
          RCDTAB(JNUC,INUC)= RDIF
          RCDTAB(INUC,JNUC)= RDIF
         END DO
        END DO
C
C SETUP IFINDX
C
        DO INUC=1,NNUC
         IFINDX(INUC)=IFUIDT(IDINDX(INUC))
        END DO
C
C DETERMINE ONSITE OVERLAP INTEGRALS FOR RADIALLY SYMMETRIC
C FUNCTIONS R**L * EXP(-ALP*R**2) STORED IN OVLTAB
C SETUP INVERSE ONSITE OVERLAP TABLE FOR ACTUAL BASIS FUNCTIONS
C IN OVLDIA
C
        THRD=1.0D0/3.0D0
        DO IFNCT=1,NFNCT
         CALL OVLONS(IFNCT,OVLTAB(1,1,IFNCT))
         DO L=0,LSYMMAX(IFNCT)
          L1=L+1
          LS=(L1*(L1+1))/2
          FAC=1.0D0/(2*L+1)
          DO ICON=1,N_CON(L1,IFNCT)
           DO IS=1,LS
            OVLDIA(1,IS,ICON,L1,IFNCT)=FAC*OVLTAB(ICON,L1,IFNCT)
           END DO
           IF (L .EQ. 2) THEN
            DO IS=4,6
             OVLDIA(1,IS,ICON,L1,IFNCT)=OVLDIA(1,IS,ICON,L1,IFNCT)*THRD
            END DO
           END IF
           DO IS=1,LS
            OVLDIA(1,IS,ICON,L1,IFNCT)=1.0D0/OVLDIA(1,IS,ICON,L1,IFNCT)
            OVLDIA(2,IS,ICON,L1,IFNCT)=SQRT(OVLDIA(1,IS,ICON,L1,IFNCT))
           END DO
          END DO
         END DO
        END DO
C
C START ASSIGNING POINTS TO IDENTITY MEMBERS
C
        NMSH=0
        DO 800 IID=1,NIDENT
         IFNCT=IFUIDT(IID)
         PRINT '(2(A,I5))','IDENTITY MEMBER ',IID,
     &                     ', MULTIPLICITY: ',NMULT(IID)
         PRINT '(A,F11.2)','NUCLEAR CHARGE: ',ZNUC(IFNCT)
         NPTAT=0
C
C DETERMINE DISTANCE OF EACH ATOM FROM ID MEMBER, GET INATM
C
         RNNS= 1.0D30
         INATM= 0
         DO INUC= 1,NNUC
          RSQR(INUC)= (RIDT(1,IID)-RNUC(1,INUC))**2
     &               +(RIDT(2,IID)-RNUC(2,INUC))**2
     &               +(RIDT(3,IID)-RNUC(3,INUC))**2
          IF (RSQR(INUC) .GT. TOLER) THEN
           RNNS= MIN(RNNS,RSQR(INUC))
          ELSE
           INATM= INUC
          END IF
         END DO
         IF (INATM .EQ. 0) THEN
          write(6,*)'DVPMESH: RIDT <--> RNUC INCONSISTENCY'
          CALL STOPIT
         END IF
         RRSMAX= RNNFAC*RNNFAC*RNNS+TOLER
         RASMAX= ANNFAC*ANNFAC*RNNS+TOLER
C
C RADIAL MESH: SPLIT INTO TWO PARTS
C CORE REGION IS DEFINED BY RHO*POT/R**2 >= 10
C WITHIN THIS REGION, OPTIMIZATION BASED ON FREE ATOM
C
         ERROR=ERRCOR
         IF (LSYWT) ERROR=ERROR/NMULT(IID)
         RDN=  0.0D0
         RUP= 10.0D0
         DO I= 1,40
          RCOR= 0.5D0*(RUP+RDN)
          CALL RPFIT(IFNCT,RCOR,1.0D0/RCOR,RHOF,POTF)
          IF (ABS(RHOF*POTF) .LT. 10*RCOR*RCOR) THEN
           RUP= RCOR
          ELSE
           RDN= RCOR
          END IF
         END DO
         RA= 0.0D0
         CALL OPTRMSH(MAXRAD,IFNCT,RA,RCOR,OVLTAB(1,1,IFNCT),ERROR,
     &                NR1,RRAD,WRAD)
C
C NOW THE VALENCE REGION
C FIRST CALL OPTVMSH WITH (IOPT .EQ. 1) TO GET AN IDEA ABOUT HOW
C MANY POINTS MINIMUM WE'LL NEED. THEN INCREASE NUMBER OF POINTS,
C FIND SHORTEST-RANGED MESH BY CALLING OPTVMSH WITH (IOPT .EQ. 2)
C AND CHECK IF ALL DIAGONAL OVERLAPS ARE CALCULATED 
C SUFFICIENTLY ACCURATE USING AN ANGULAR GRID OF SIZE LMXRAD
C
         ERROR=ERRAD
         IF (LSYWT) ERROR=ERROR/NMULT(IID)
         NR2= 0
         CALL OPTVMSH(1,1,MAXRAD-NR1,IFNCT,RCOR,OVLTAB(1,1,IFNCT),
     &                ERROR,NR2,RRAD(NR1+1),WRAD(NR1+1),RTRAFO)
         RTRINIT= RTRAFO
         CALL LEBEDEV(MAXANG,LMXRAD,NANG,RANG,WANG)
         DO INUC= 1,NNUC
          IF (RSQR(INUC) .LT. RRSMAX) THEN
           JFNCT= IFINDX(INUC)
           DO L= 0,LSYMMAX(JFNCT)
            L1= L+1
            LS= (L1*(L1+1))/2
            DO ICON= 1,N_CON(L1,JFNCT)
             DO IS= 1,LS
              OVLSUM(1,IS,ICON,L1,INUC)= 1.0D30
             END DO
            END DO
           END DO
          END IF
         END DO
         NUP= -1
         NLW= NR2
  100    CONTINUE
          RTRAFO= RTRINIT
          CALL OPTVMSH(2,1,MAXRAD-NR1,IFNCT,RCOR,OVLTAB(1,1,IFNCT),
     &                 ERROR,NR2,RRAD(NR1+1),WRAD(NR1+1),RTRAFO)
          DO INUC= 1,NNUC
           IF (RSQR(INUC) .LT. RRSMAX) THEN
            JFNCT= IFINDX(INUC)
            DO L= 0,LSYMMAX(JFNCT)
             L1= L+1
             LS= (L1*(L1+1))/2
             DO ICON= 1,N_CON(L1,JFNCT)
              DO IS= 1,LS
               OVLSUM(2,IS,ICON,L1,INUC)= 0.0D0
              END DO
             END DO
            END DO
           END IF
          END DO
          DO IRAD= 1,NR2
           CALL SPCPART(MAXANG,NANG,RANG,WANG,RRAD(NR1+IRAD),
     &                  WRAD(NR1+IRAD),RIDT(1,IID),IFNCT,INATM,NNUC,
     &                  RNUC,IFINDX,RCDTAB,MODPART,NDECAY,STOMSH)
           DO INUC= 1,NNUC
            IF (RSQR(INUC) .LT. RRSMAX) THEN
             JFNCT= IFINDX(INUC)
             DO IANG= 0,NANG-1,NSPEED
              NPV= MIN(NSPEED,NANG-IANG)
              DO IPV= 1,NPV
               PTS(IPV,1)= STOMSH(1,IANG+IPV,2)-RNUC(1,INUC)
               PTS(IPV,2)= STOMSH(2,IANG+IPV,2)-RNUC(2,INUC)
               PTS(IPV,3)= STOMSH(3,IANG+IPV,2)-RNUC(3,INUC)
              END DO
              CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,JFNCT,GRAD)
              IF (IUPDAT) THEN
               DO L= 0,LSYMMAX(JFNCT)
                L1= L+1
                LS= (L1*(L1+1))/2
                DO ICON= 1,N_CON(L1,JFNCT)
                 IF (ICOUNT(ICON,L1)) THEN
                  DO IS= 1,LS
                   SUM= 0.0D0
                   DO IPV= 1,NPV
                    BF2= GRAD(IPV,1,IS,ICON,L1)**2
                    SUM= SUM+STOMSH(4,IANG+IPV,2)*BF2
                   END DO
                   OVLSUM(2,IS,ICON,L1,INUC)= OVLSUM(2,IS,ICON,L1,INUC)
     &                                       +SUM
                  END DO
                 END IF
                END DO
               END DO
              END IF
             END DO
            END IF
           END DO
          END DO
C
C CHECK ACCURACY
C
          ISUCC= .TRUE.
          DO INUC= 1,NNUC
           IF (RSQR(INUC) .LT. RRSMAX) THEN
            JFNCT= IFINDX(INUC)
            DO L= 0,LSYMMAX(JFNCT)
             L1= L+1
             LS= (L1*(L1+1))/2
             DO ICON= 1,N_CON(L1,JFNCT)
              DO IS= 1,LS
               DIFF= OVLSUM(2,IS,ICON,L1,INUC)-OVLSUM(1,IS,ICON,L1,INUC)
               IF (ABS(DIFF)*OVLDIA(1,IS,ICON,L1,JFNCT) .GT. ERROR) THEN
                ISUCC= .FALSE.
                GOTO 150
               END IF
              END DO
             END DO
            END DO
           END IF
          END DO
C
C IF WE DON'T HAVE AN UPPER LIMIT YET, UPDATE REFERENCE INTEGRALS
C
  150     IF (NUP .EQ. -1) THEN
           DO INUC= 1,NNUC
            IF (RSQR(INUC) .LT. RRSMAX) THEN
             JFNCT= IFINDX(INUC)
             DO L= 0,LSYMMAX(JFNCT)
              L1= L+1
              LS= (L1*(L1+1))/2
              DO ICON= 1,N_CON(L1,JFNCT)
               DO IS= 1,LS
                OVLSUM(1,IS,ICON,L1,INUC)= OVLSUM(2,IS,ICON,L1,INUC)
               END DO
              END DO
             END DO
            END IF
           END DO
          END IF
C
          IF (ISUCC) THEN
C
C THIS MESH WAS ACCURATE ENOUGH
C
           IF (NUP .EQ. -1) NUP= NR2
           IF (NUP-NLW .GE. 2) THEN
            NUP= NR2
            NR2= (NUP+NLW)/2
            GOTO 100
           END IF
          ELSE
C
C THIS MESH WAS NOT ACCURATE ENOUGH
C
           IF (NUP .EQ. -1) THEN
            NR2= NR2+10
            GOTO 100
           ELSE
            NLW= NR2
            IF (NUP-NLW .GE. 2) THEN
             NR2= (NUP+NLW)/2
             GOTO 100
            ENDIF
           ENDIF
          ENDIF
         CONTINUE
C
C THE NUMBER OF POINTS TO USE IS NUP
C
         NR2= NUP
         RTRAFO= RTRINIT
         CALL OPTVMSH(2,1,MAXRAD-NR1,IFNCT,RCOR,OVLTAB(1,1,IFNCT),
     &                ERROR,NR2,RRAD(NR1+1),WRAD(NR1+1),RTRAFO)
C
C DONE WITH RADIAL MESH
C
         NRAD= NR1+NR2
         PRINT '(A,I5)','NUMBER OF RADIAL POINTS: ',NRAD
         ERROR= ERRANG/NRAD
         IF (LSYWT) ERROR= ERROR/NMULT(IID)
C
C FOR EACH RADIAL POINT, OPTIMIZE ANGULAR MESH
C
         DO 700 IRAD=1,NRAD
C
C ZERO OVERLAP ARRAYS
C
          DO LA= 0,LSYMMAX(IFNCT)
           L1A= LA+1
           LSA= (L1A*(L1A+1))/2
           DO ICONA= 1,N_CON(L1A,IFNCT)
            DO LB= LA+1,LSYMMAX(IFNCT)
             L1B= LB+1
             LSB= (L1B*(L1B+1))/2
             DO ICONB= 1,N_CON(L1B,IFNCT)
              DO ISA= 1,LSA
               DO ISB= 1,LSB
                ONSOVL(1,ISB,ISA,ICONB,L1B,ICONA,L1A)= 1.0D30
               END DO
              END DO
             END DO
            END DO
           END DO
          END DO
          DO INUC= 1,NNUC
           IF (RSQR(INUC) .LT. RASMAX) THEN
            JFNCT= IFINDX(INUC)
            DO L= 0,LSYMMAX(JFNCT)
             L1= L+1
             LS= (L1*(L1+1))/2
             DO ICON= 1,N_CON(L1,JFNCT)
              DO IS= 1,LS
               OVLSUM(1,IS,ICON,L1,INUC)= 1.0D30
              END DO
             END DO
            END DO
           END IF
          END DO
C
C LOOP OVER ALL ANGULAR MESHES UNTIL WE FIND A GOOD ONE
C
          LCUR=0
          NSTO=0
  300     CONTINUE
           CALL LEBEDEV(MAXANG,LCUR,NANG,RANG,WANG)
           IF (LCUR .GT. LMXANG) THEN
            ISUCC= .TRUE.
            GOTO 500
           END IF
           NTOT=NMSH+NANG   
           IF (NTOT .GT. MAX_PTS) THEN
            write(6,*)'DVPMESH: MAX_PTS MUST BE AT LEAST: ',NTOT
            CALL STOPIT
           END IF
C
C SPACE PARTITIONING
C
           CALL SPCPART(MAXANG,NANG,RANG,WANG,RRAD(IRAD),WRAD(IRAD),
     &                  RIDT(1,IID),IFNCT,INATM,NNUC,RNUC,IFINDX,
     &                  RCDTAB,MODPART,NDECAY,STOMSH)
C
C CHECK IF ANGULAR RESOLUTION IS GOOD ENOUGH FOR OVERLAP INTEGRATION
C START BY DETERMINING NONDIAGONAL ONSITE ELEMENTS
C
           DO LA=0,LSYMMAX(IFNCT)
            L1A=LA+1
            LSA=(L1A*(L1A+1))/2
            DO ICONA=1,N_CON(L1A,IFNCT)
             DO LB=LA+1,LSYMMAX(IFNCT)
              L1B=LB+1
              LSB=(L1B*(L1B+1))/2
              DO ICONB=1,N_CON(L1B,IFNCT)
               DO ISA=1,LSA
                DO ISB=1,LSB
                 ONSOVL(2,ISB,ISA,ICONB,L1B,ICONA,L1A)= 0.0D0
                END DO
               END DO
              END DO
             END DO
            END DO
           END DO
C
           DO 440 IANG=0,NANG-1,NSPEED
            NPV=MIN(NSPEED,NANG-IANG)
            DO IPV=1,NPV
             PTS(IPV,1)=STOMSH(1,IANG+IPV,2)-RIDT(1,IID)
             PTS(IPV,2)=STOMSH(2,IANG+IPV,2)-RIDT(2,IID)
             PTS(IPV,3)=STOMSH(3,IANG+IPV,2)-RIDT(3,IID)
            END DO
            CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
            IF (IUPDAT) THEN
             DO LA=0,LSYMMAX(IFNCT)
              L1A=LA+1
              LSA=(L1A*(L1A+1))/2
              DO ICONA=1,N_CON(L1A,IFNCT)
               IF (ICOUNT(ICONA,L1A)) THEN
                DO LB=LA+1,LSYMMAX(IFNCT)
                 L1B=LB+1
                 LSB=(L1B*(L1B+1))/2
                 DO ICONB=1,N_CON(L1B,IFNCT)
                  IF (ICOUNT(ICONB,L1B)) THEN
                   DO ISA=1,LSA
                    DO ISB=1,LSB
                     SUM= 0.0D0
                     DO IPV=1,NPV
                      SUM=SUM+STOMSH(4,IANG+IPV,2)
     &                       *GRAD(IPV,1,ISA,ICONA,L1A)
     &                       *GRAD(IPV,1,ISB,ICONB,L1B)
                     END DO
                     ONSOVL(2,ISB,ISA,ICONB,L1B,ICONA,L1A)= 
     &                ONSOVL(2,ISB,ISA,ICONB,L1B,ICONA,L1A)+SUM
                    END DO
                   END DO
                  END IF
                 END DO
                END DO
               END IF
              END DO
             END DO
            END IF
  440      CONTINUE
C
C NEXT: ALL DIAGONAL ELEMENTS OF OVERLAP MATRIX
C
           DO 480 INUC=1,NNUC
            IF (RSQR(INUC) .LT. RASMAX) THEN
             JFNCT= IFINDX(INUC)
             DO L=0,LSYMMAX(JFNCT)
              L1=L+1
              LS=(L1*(L1+1))/2
              DO ICON=1,N_CON(L1,JFNCT)
               DO IS=1,LS
                OVLSUM(2,IS,ICON,L1,INUC)= 0.0D0 
               END DO
              END DO
             END DO
             DO IANG=0,NANG-1,NSPEED
              NPV=MIN(NSPEED,NANG-IANG)
              DO IPV=1,NPV
               PTS(IPV,1)=STOMSH(1,IANG+IPV,2)-RNUC(1,INUC)
               PTS(IPV,2)=STOMSH(2,IANG+IPV,2)-RNUC(2,INUC)
               PTS(IPV,3)=STOMSH(3,IANG+IPV,2)-RNUC(3,INUC)
              END DO
              CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,JFNCT,GRAD)
              IF (IUPDAT) THEN
               DO L=0,LSYMMAX(JFNCT)
                L1=L+1
                LS=(L1*(L1+1))/2
                DO ICON=1,N_CON(L1,JFNCT)
                 IF (ICOUNT(ICON,L1)) THEN
                  DO IS=1,LS
                   SUM= 0.0D0
                   DO IPV=1,NPV
                    BF=GRAD(IPV,1,IS,ICON,L1)
                    BF2=BF*BF
                    SUM=SUM+STOMSH(4,IANG+IPV,2)*BF2
                    BF2=BF2*OVLDIA(1,IS,ICON,L1,JFNCT)
                    STOMSH(6,IANG+IPV,2)=MAX(STOMSH(6,IANG+IPV,2),BF2)
                   END DO
                   OVLSUM(2,IS,ICON,L1,INUC)= OVLSUM(2,IS,ICON,L1,INUC)
     &                                       +SUM
                  END DO
                 END IF
                END DO
               END DO
              END IF
             END DO
            END IF
  480      CONTINUE
C
C CHECK IF CHANGE IN ONSOVL AND OVLSUM IS LESS THAN ERROR*(TOTAL INTEGRAL)
C
           ISUCC= .TRUE.
           DO LA=0,LSYMMAX(IFNCT)
            L1A=LA+1
            LSA=(L1A*(L1A+1))/2
            DO ICONA=1,N_CON(L1A,IFNCT)
             DO LB=LA+1,LSYMMAX(IFNCT)
              L1B=LB+1
              LSB=(L1B*(L1B+1))/2
              DO ICONB=1,N_CON(L1B,IFNCT)
               DO ISA=1,LSA
                DO ISB=1,LSB
                 DIFF=ONSOVL(2,ISB,ISA,ICONB,L1B,ICONA,L1A)
     &               -ONSOVL(1,ISB,ISA,ICONB,L1B,ICONA,L1A)
                 SREF=OVLDIA(2,ISA,ICONA,L1A,IFNCT)
     &               *OVLDIA(2,ISB,ICONB,L1B,IFNCT)
                 IF (ABS(DIFF)*SREF .GT. ERROR) THEN
                  ISUCC= .FALSE.
                  GOTO 500
                 END IF
                END DO
               END DO
              END DO
             END DO
            END DO
           END DO
C
           DO INUC=1,NNUC
            IF (RSQR(INUC) .LT. RASMAX) THEN
             JFNCT= IFINDX(INUC)
             DO L=0,LSYMMAX(JFNCT)
              L1=L+1
              LS=(L1*(L1+1))/2
              DO ICON=1,N_CON(L1,JFNCT)
               DO IS=1,LS
                DIFF=OVLSUM(2,IS,ICON,L1,INUC)-OVLSUM(1,IS,ICON,L1,INUC)
                DIFF=ABS(DIFF)*OVLDIA(1,IS,ICON,L1,JFNCT)
                IF (DIFF .GT. ERROR) THEN
                 ISUCC= .FALSE.
                 GOTO 500
                END IF
               END DO
              END DO
             END DO
            END IF
           END DO
C
C SUCCESS        --> USE THE LAST MESH
C LCUR >= LMXANG --> USE THE CURRENT MESH
C PURGE SYMMETRY-EQUIVALENT POINTS
C
  500      IF (ISUCC .OR. (LCUR .GE. LMXANG)) THEN
            IF (ISUCC) THEN
             NANG=NSTO
             INDX=1
            ELSE
             INDX=2
            END IF
            RRC= 1.0D0/RRAD(IRAD)
            DO IANG=1,NANG
             DO IC=1,3
              RANG(IC,IANG)=RIDT(IC,IID)
     &                     +(STOMSH(IC,IANG,INDX)-RIDT(IC,IID))*RRC
             END DO
             WANG(IANG)=STOMSH(4,IANG,INDX)
            END DO
C
C FOR ALL NON-IDENTITY SYMMETRY OPERATIONS:
C * CHECK IF SYMMETRY OPERATION LEAVES RIDT UNCHANGED
C * IF SO, CHECK IF POINTS CAN BE REDUCED
C
            DO 650 IGRP=2,NGRP
             DO I=1,3
              RTMP=0.0D0
              DO J=1,3
               RTMP=RTMP+RMAT(J,I,IGRP)*RIDT(J,IID)
              END DO
              IF (ABS(RTMP-RIDT(I,IID)) .GT. TOLER) GOTO 650
             END DO
             DO IANG=1,NANG
              IF (WANG(IANG) .GT. 0.0D0) THEN
               DO 630 JANG=IANG+1,NANG
                IF (WANG(JANG) .GT. 0.0D0) THEN 
                 DO I=1,3
                  RTMP=0.0D0
                  DO J=1,3
                   RTMP=RTMP+RMAT(J,I,IGRP)*RANG(J,JANG)
                  END DO
                  IF (ABS(RTMP-RANG(I,IANG)) .GT. TOLER) GOTO 630
                 END DO
                 WANG(IANG)=WANG(IANG)+WANG(JANG)
                 WANG(JANG)= -1.0D0
                 STOMSH(5,IANG,INDX)=STOMSH(5,IANG,INDX)
     &                              +STOMSH(5,JANG,INDX)
                 STOMSH(6,IANG,INDX)=STOMSH(6,IANG,INDX)
     &                              +STOMSH(6,JANG,INDX)
                END IF
  630          CONTINUE
              END IF
             END DO
  650       CONTINUE
            NTOT=0
            DO IANG=1,NANG
             IF (WANG(IANG) .GT. 0.0D0) THEN
              WANG(IANG)=WANG(IANG)*NMULT(IID)
              POTMX=MAX(1.0D0,STOMSH(5,IANG,INDX))
              POTMX=POTMX*STOMSH(6,IANG,INDX)
              IF (WANG(IANG)*POTMX .GT. ERROR) THEN
               NTOT=NTOT+1
               RMSH(1,NMSH+NTOT)= STOMSH(1,IANG,INDX)
               RMSH(2,NMSH+NTOT)= STOMSH(2,IANG,INDX)
               RMSH(3,NMSH+NTOT)= STOMSH(3,IANG,INDX)
               WMSH(NMSH+NTOT)= WANG(IANG)
              END IF
             END IF
            END DO
            NMSH=NMSH+NTOT
            NPTAT=NPTAT+NTOT
            GOTO 680
           END IF
C
C FAILURE: MOVE STOMSH, ONSOVL, AND OVLSUM AND INCREASE LCUR
C
           NSTO=NANG
           DO IANG=1,NSTO
            DO I=1,6
             STOMSH(I,IANG,1)=STOMSH(I,IANG,2)
            END DO
           END DO
C
           DO LA=0,LSYMMAX(IFNCT)
            L1A=LA+1
            LSA=(L1A*(L1A+1))/2
            DO ICONA=1,N_CON(L1A,IFNCT)
             DO LB=LA+1,LSYMMAX(IFNCT)
              L1B=LB+1
              LSB=(L1B*(L1B+1))/2
              DO ICONB=1,N_CON(L1B,IFNCT)
               DO ISA=1,LSA
                DO ISB=1,LSB
                 ONSOVL(1,ISB,ISA,ICONB,L1B,ICONA,L1A)=
     &            ONSOVL(2,ISB,ISA,ICONB,L1B,ICONA,L1A)
                END DO
               END DO
              END DO
             END DO
            END DO
           END DO
C
           DO INUC=1,NNUC
            IF (RSQR(INUC) .LT. RASMAX) THEN
             JFNCT= IFINDX(INUC)
             DO L=0,LSYMMAX(JFNCT)
              L1=L+1
              LS=(L1*(L1+1))/2
              DO ICON=1,N_CON(L1,JFNCT)
               DO IS=1,LS
                OVLSUM(1,IS,ICON,L1,INUC)= OVLSUM(2,IS,ICON,L1,INUC)
               END DO
              END DO
             END DO
            END IF
           END DO
C
           LCUR= LCUR+1
           GOTO 300
  680     CONTINUE
  700    CONTINUE
C
         PRINT '(A,I5,A,I10)','NUMBER OF POINTS FOR ID MEMBER ',
     &                        IID,':',NPTAT
         PRINT '(A)',' '
  800   CONTINUE
        PRINT '(A,I15)','TOTAL NUMBER OF MESH POINTS: ',NMSH
        PRINT '(A)',' '
C
C REORDER MESH POINTS
C
        PRINT '(A)','REORDERING MESH POINTS'
        CALL GTTIME(TIME2)
        CALL REORMSH
        CALL GTTIME(TIME3)
C
C TEST MESH (CHECK NUMERICAL OVERLAP MATRIX)
C
        write(6,*)' '
        PRINT '(A)','COMPARISON OF ANALYTICAL AND NUMERICAL OVERLAPS'
        PRINT '(A)','-----------------------------------------------'
        CALL TSTMSH
        write(6,*)' '
C       
C WRITE VMOLD AND LEAVE
C
        OPEN(99,FILE='OLDMESH',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(99)
        WRITE(99)NMSH,ICALC
        WRITE(99)((RMSH(J,I),J=1,3),I=1,NMSH)
        WRITE(99)(WMSH(I),I=1,NMSH)
        CLOSE(99)
  900   CALL GTTIME(TIME4)
        CALL TIMOUT('MESH REORDERING:                   ',TIME3-TIME2)
        CALL TIMOUT('MESH READING / TESTING + WRITING:  ',TIME4-TIME3)
        CALL TIMOUT('CONSTRUCTION OF VARIATIONAL MESH:  ',TIME4-TIME1)
c
c mpi: send mesh to all nodes
c
#ifdef MPI
c
        CALL SENDDATA(102)
c
#endif
c
C
C DEALLOCATE LOCAL ARRAYS
C
       DEALLOCATE(STOMSH,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING STOMSH'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING PTS'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING STOMSH'
       DEALLOCATE(OVLTAB,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING OVLTAB'
       DEALLOCATE(OVLDIA,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING OVLDIA'
       DEALLOCATE(ONSOVL,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING ONSOVL'
       DEALLOCATE(OVLSUM,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING OVLSUM'
       DEALLOCATE(RRAD,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING RRAD'
       DEALLOCATE(RANG,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING RANG'
       DEALLOCATE(WRAD,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING WRAD'
       DEALLOCATE(WANG,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING WANG'
       DEALLOCATE(RNUC,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING RNUC'
       DEALLOCATE(REQV,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING REQV'
       DEALLOCATE(RSQR,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING RSQR'
       DEALLOCATE(RCDTAB,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING RCDTAB'
       DEALLOCATE(NMULT,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING NMULT'
       DEALLOCATE(IDINDX,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING IDINDX'
       DEALLOCATE(IFINDX,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING IFINDX'
       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'DVPMESH:ERROR DEALLOCATING ICOUNT'
C

        RETURN
       END
