C UTEP Electronic Structure Lab (2020)
C
       SUBROUTINE POTGRID
C
C POTRHOGRID VERSION DIRK POREZAG JUNE 1998. 
C CALCULATION OF CHARGES WITHIN A SPHERE moved into atomsph JK 3/99
C * POTENTIAL ON A GRID -- TB 12/04
       USE XTMP2, Only : PHIG
       use common2,only : RCNT, NCNT, IGGA, ISPN, NSPN, EFIELD
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:57 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I1, I2, IBEG, ICNT, IFORM, IGRAD, ILOOP, IPTS, IS,
     & ITYPE, IUNIT, IX, IY, IZ, J, MODE, NATOM, NGRAD, NGRID, NLOOP,
     & NMSH, NSTORE
       REAL*8 :: PHIG , SYMBOL, CHR, DENTOT, DERIV, EFHERE, EXHERE,
     & FACT1, FACT4, FACTOR, PI, POT, POTIN, RBAS, TIME1, TMP, VOL, X,
     & Y, Z, RHOG
       SAVE
       PARAMETER (MAXSPH=500)
       PARAMETER (MAXRAD=1000)
       PARAMETER (MAXANG=200)
C
       LOGICAL ICOUNT,EXIST
       LOGICAL POTEN
       LOGICAL ELF
       LOGICAL LGGA
       CHARACTER*20 FNAME(6),FORMSTR
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       DIMENSION NGRID(3),DERIV(3)
       DIMENSION RBAS(3,4),TMP(3)
C      COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,KRHOG,MXSPN)
C    &  ,PHIG(MAX_PTS,2)
       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
C
C
C FOR DENSITY EVALUATIONS, USE OLD SCHEME FOR DENSITY CALCULATION
C OTHERWISE, GET DENSITY FROM COUPOT
C
       PI=4.0D0*DATAN(1.0D0)
       PRINT '(A)',' '
       PRINT '(A)','DENSITY/POTENTIAL GRIDS'
       CALL GTTIME(TIME1)
       POTEN =.FALSE.
       NSTORE=NMSH
C
C READ IN NECESSARY INPUT DATA
C
         FORMSTR= ' '

         FNAME(1)='POTGRID'
         INQUIRE(FILE=FNAME(1),EXIST=EXIST)
         IF (EXIST) POTEN=.TRUE.
         IF (.NOT.EXIST) THEN
           FORMSTR= ' --> NOTHING TO DO'
           PRINT '(2A)','POTENTIAL GRID  ',FORMSTR
         END IF
C
C DETERMINE IF THERE IS ANYTHING TO DO 
C

        IF ((.NOT.POTEN)) GOTO 900
C
C READ INPUT DATA
C
        IF (POTEN)FNAME(1)="POTGRID"
        OPEN(72,FILE=FNAME(1),FORM='FORMATTED',STATUS='OLD')
        REWIND(72)
C
C CHECK IF THIS FILE IS EMPTY
C IF YES, CREATE A DEFAULT ONE
C
   5    I=1
        READ(72,*,END=10) ITYPE,IFORM
        I=0
        REWIND(72)
   10   IF (I .EQ. 1) THEN
         DO I=1,3
          RBAS(I,1)=  1.0D30
          RBAS(I,2)= -1.0D30
         END DO
         DO ICNT=1,NCNT
          DO I=1,3
           RBAS(I,1)= MIN(RBAS(I,1),RCNT(I,ICNT))
           RBAS(I,2)= MAX(RBAS(I,2),RCNT(I,ICNT))
          END DO
         END DO
         DO I=1,3
          RBAS(I,1)= RBAS(I,1)-5.0D0
          RBAS(I,2)= RBAS(I,2)+5.0D0
          NGRID(I)= (RBAS(I,2)-RBAS(I,1))/0.5D0+2
          RBAS(I,2)= (RBAS(I,2)-RBAS(I,1))/(NGRID(I)-1)
         END DO
         REWIND(72)
         WRITE(72,1010) 1,2,'Grid mode, formatted file'
         WRITE(72,1020) (NGRID(I), I=1,3),    'Number of grid points'
         WRITE(72,1030) (RBAS(I,1), I=1,3),   'Origin' 
         WRITE(72,1030) RBAS(1,2),0.0D0,0.0D0,'Basis vector 1' 
         WRITE(72,1030) 0.0D0,RBAS(2,2),0.0D0,'Basis vector 1' 
         WRITE(72,1030) 0.0D0,0.0D0,RBAS(3,2),'Basis vector 1' 
 1010    FORMAT(2(I6,1X),8X,A)
 1020    FORMAT(3(I6,1X),1X,A)
 1030    FORMAT(3(F7.3,1X),1X,A)
         CLOSE(72)
         OPEN(72,FILE=FNAME(1),FORM='FORMATTED',STATUS='OLD')
         REWIND(72)
        END IF
C
C GRID INPUT 
C
        READ(72,*,END=880) ITYPE,IFORM
        IF((ITYPE.EQ.1).AND.(IFORM.EQ.3)) ELF=.TRUE.
        IF (ITYPE .GT. 2) ITYPE= 2
        IF (ITYPE .LT. 1) ITYPE= 1
        IF (IFORM .GT. 2) IFORM= 2
        IF (IFORM .LT. 1) IFORM= 1
        IF (ITYPE .EQ. 2) IFORM= 2
        FORMSTR='UNFORMATTED'
        IF (IFORM .EQ. 2) FORMSTR='FORMATTED'
C
C GET NUMBER OF GRID POINTS, ORIGIN, AND BASIS VECTORS
C
        READ(72,*,END=880)(NGRID(J), J=1,3)
        IF (ITYPE .EQ. 1) THEN
         DO I=1,4
          READ(72,*,END=880)(RBAS(J,I), J=1,3)
         END DO
        END IF
        DO I=1,3
         IF (NGRID(I) .LT. 1) THEN
          PRINT *,'RHOGRID: NUMBER OF GRID POINTS MUST BE >= 1'
          GOTO 890
         END IF
        END DO
        IF (NGRID(3) .GT. MAX_PTS) THEN
         PRINT *,'POTGRID: MAX_PTS MUST BEAT LEAST: ',NGRID(3)
         PRINT *,'SKIPPING GRID EVALUATION '
         GOTO 890
        END IF
        NLOOP=NGRID(1)*NGRID(2)
        IF(POTEN) THEN
          FNAME(2)='POT'
          IF (NSPN.GT.1) THEN
              FNAME(2)=FNAME(2)(1:3)//'UP'
              FNAME(3)=FNAME(2)(1:3)//'DN'
          END IF
        ENDIF
               
   
C
C OPEN OUTPUT FILES, WRITE HEADER
C

        IBEG=75
        DO IS=1,NSPN
        IUNIT=IBEG+IS
   20   OPEN(IUNIT,FILE=FNAME(IS+1),FORM=FORMSTR,STATUS='UNKNOWN')
        REWIND(IUNIT)
        IF (IFORM .EQ. 1) THEN
         WRITE(IUNIT) ITYPE,NSPN
         WRITE(IUNIT)(NGRID(J), J=1,3),NGRID(3)
         IF (ITYPE .EQ. 1) THEN
          WRITE(IUNIT)((RBAS(J,I), J=1,3), I=1,4)
         END IF
         WRITE(IUNIT) NSPN
         DO ISPN=1,NSPN
          X=ISPN
          WRITE(IUNIT) ISPN,ISPN,ISPN,ISPN,X
         END DO
        ELSE
             WRITE(IUNIT,*)'CLUSTER OUTPUT'
             IF(NSPN.EQ.1) WRITE(IUNIT,*)'SCF TOTAL POTENTIAL (ANG)'
             IF(NSPN.EQ.2) WRITE(IUNIT,*)'SCF POTENTIAL (ANG)'
         OPEN(87,FILE='XMOL.DAT')
         READ(87,*) NATOM
         READ(87,*) 
         WRITE(IUNIT,'(1X,I10,3F20.12)') NATOM,(RBAS(J,1),J=1,3)
         DO I=1,3
         WRITE(IUNIT,'(1X,I10,3F20.12)') NGRID(I),(RBAS(J,I+1),J=1,3)
         ENDDO
         DO I=1,NATOM
           READ(87,*)IZ, X, Y, Z
           CHR=REAL(IZ)
           WRITE(IUNIT,2002)IZ, CHR, X, Y, Z
         END DO
         CLOSE(87)
        ENDIF
        ENDDO
        
2001    FORMAT(A25)
2002    FORMAT(I6,4F16.10)
C
C LOOP FOR EACH PILE
C
        DENTOT=0.0D0
        DO 850 ILOOP=1,NLOOP
C
C SETUP POINTS
C
         NMSH=NGRID(3)
         IF (ITYPE .EQ. 1) THEN
          IY=MOD(ILOOP-1,NGRID(2))
          IX=(ILOOP-1)/NGRID(2)
          DO IZ=1,NMSH
           I=IZ-1
           RMSH(1,IZ)=RBAS(1,1)+IX*RBAS(1,2)+IY*RBAS(1,3)+I*RBAS(1,4)
           RMSH(2,IZ)=RBAS(2,1)+IX*RBAS(2,2)+IY*RBAS(2,3)+I*RBAS(2,4)
           RMSH(3,IZ)=RBAS(3,1)+IX*RBAS(3,2)+IY*RBAS(3,3)+I*RBAS(3,4)
           TMP(1)=RBAS(1,3)*RBAS(2,4)-RBAS(1,4)*RBAS(2,3)
           TMP(2)=RBAS(1,4)*RBAS(2,2)-RBAS(1,2)*RBAS(2,4)
           TMP(3)=RBAS(1,2)*RBAS(2,3)-RBAS(1,3)*RBAS(2,2)
           WMSH(IZ)=RBAS(3,2)*TMP(1)+RBAS(3,3)*TMP(2)+RBAS(3,4)*TMP(3)
c            WMSH(IZ)=0.0D0 
          END DO
         ELSE
          DO IZ=1,NMSH
           READ(72,*,END=870)(RMSH(I,IZ), I=1,3)
           WMSH(IZ)=0.0D0
          END DO
         END IF
C
C NOW: CALCULATE ELECTRONIC DENSITY
C DENSITY WILL BE STORED IN RHOG
C
         NGRAD=1
c
c
c MPI: need to send updated mesh data
c
c
          I1=IGGA(1)
          I2=IGGA(2)
          IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) THEN
              LGGA= .TRUE.
              NGRAD=10
          END IF

          CALL DENSOLD(VOL)
          IGGA(1)=I1
          IGGA(2)=I2

        IF (POTEN) THEN
c%ifdef MPI 
c        CALL SENDDATA(102)
c%endif
        CALL COUPOT1
        ENDIF

C
C
        FACTOR=0.5292D0**3
        FACT4=0.5292D0**4
        FACT1=0.5292D0
C
C
C UPDATE DATA IN RHOG
C
         DO IGRAD=1,NGRAD
          DO IPTS=1,NMSH
           RHOG(IPTS,IGRAD,1)=RHOG(IPTS,IGRAD,1)+RHOG(IPTS,IGRAD,NSPN)
          END DO 
          DO IPTS=1,NMSH
            DENTOT=DENTOT+RHOG(IPTS,1,1)*WMSH(IPTS)
          END DO
         END DO 
C

        IF (POTEN) THEN
          WRITE(6,*) 'WRITING POTENTIAL'
C THE FOLLOWING PART IS ONLY DONE IF (MODE .EQ. 2)
C CALCULATING KOHN-SHAM POTENTIAL POT 
C
         CALL GETVLXC(NMSH,RHOG,POT,POTIN)
         CALL AFPOT(NSPN,RHOG,POT,POTIN)
C
C ADD EFIELD POTENTIAL TO LOCAL POTENTIAL
C
         DO IPTS=1,NMSH
          CALL EXTPOT(RMSH(1,IPTS),EXHERE,DERIV)
          EFHERE=EFIELD(1)*RMSH(1,IPTS)+EFIELD(2)*RMSH(2,IPTS)
     &          +EFIELD(3)*RMSH(3,IPTS)
          POTIN(IPTS)=POTIN(IPTS)+POT(IPTS)+EFHERE+EXHERE
         END DO

C
C WRITE OUTPUT: POTIN   CONTAINS LOCAL +XC
C               COULOMB CONTAINS COULOMB 
C               POT     CONTAINS EXCHANGE-CORRELATION 
C
         IF (NSPN.EQ.1) THEN
          IF (ITYPE .EQ. 1) THEN
           IF (IFORM .EQ. 1) THEN
            WRITE(76)(POTIN(IPTS),COULOMB(IPTS),POT(IPTS),
     &                IPTS=1,NMSH)
           ELSE
            WRITE(76,9010)(POTIN(IPTS)+COULOMB(IPTS),IPTS=1,NMSH)
           END IF
          ELSE
           DO IPTS=1,NMSH
            WRITE(76,9020)(RMSH(I,IPTS), I=1,3)
            WRITE(76,9010) POTIN(IPTS),COULOMB(IPTS),POT(IPTS)
           END DO
          END IF
         ELSE                 ! SPIN-POLARIZED
          IF (ITYPE .EQ. 1) THEN
           IF (IFORM .EQ. 1) THEN
            WRITE(76)(POTIN(IPTS),COULOMB(IPTS),
     &                POT(IPTS),POT(IPTS+NMSH), IPTS=1,NMSH)
           ELSE
            WRITE(76,9010)(POTIN(IPTS)+COULOMB(IPTS),
     &          IPTS=1,NMSH)
            WRITE(77,9010)(POTIN(IPTS+NMSH)+COULOMB(IPTS),
     &          IPTS=1,NMSH)
           END IF
          ELSE    ! ITYPE 2
           DO IPTS=1,NMSH
            WRITE(76,9020)(RMSH(I,IPTS), I=1,3)
            WRITE(76,9010) POTIN(IPTS),COULOMB(IPTS),
     &                     POT(IPTS),POT(IPTS+NMSH)
           END DO
          END IF
         END IF

        END IF
C
  850   CONTINUE
        WRITE(6,*)'DENSITY IN RHOPOT INTEGRATES TO ', DENTOT
        CLOSE(74)
        IF (NSPN.EQ.2) CLOSE(75)
        CLOSE(76)
        IF (NSPN.EQ.2) CLOSE(77)
        CLOSE(81)
        IF (NSPN.EQ.2) CLOSE(82)
        GOTO 890
C
C ERROR HANDLING
C        
  870   CLOSE(76)
        IF (NSPN.EQ.2) CLOSE(77)
  880   PRINT *,'ERROR IN INPUT FILE, SKIPPING MODE ',MODE
  890   CLOSE(72)
  900  CONTINUE
       NMSH=NSTORE
       RETURN
 9010  FORMAT(3(1X,E20.12))
 9020  FORMAT(3(1X,F20.12))
       END
