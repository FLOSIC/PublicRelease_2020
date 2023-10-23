C UTEP Electronic Structure Lab (2020)
C
C ****************************************************************
C
       SUBROUTINE POTRHOGRD
C
C POTRHOGRID VERSION DIRK POREZAG JUNE 1998. 
C CALCULATION OF CHARGES WITHIN A SPHERE moved into atomsph JK 3/99
C * DENSITY AND POTENTIAL ON A GRID OF POINTS
C
       use common2,only : RCNT, NIDENT, NCNT, IGGA, ISPN, NSPN, EFIELD
       use common3,only : RMAT
       use common7,only : MODDEN, GAUSS_CUT
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:56 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MAXANG, MAXSPH, I, I1, I2, IAT, ICNT, IFORM, IGRAD,
     & ILOOP, IPTS, IS, ITYPE, IUNIT, IX, IY, IZ, J, MAXRAD, MODE,
     & NATOM, NGRAD, NGRID, NLOOP, NPTS
       REAL*8 :: SYMBOL , ANGLE, CENTER, CHR, COULOMB, DERIV, DOMEGA,
     & EFHERE, EXHERE, FACTOR, GRAD, POT, POTIN, PTS, RBAS, RDIS, RHOG,
     & RVECA, TIME1, VLOC, VOL, WTRAD, X, XRAD, Y, Z
       SAVE
       PARAMETER (MAXSPH=500)
       PARAMETER (MAXRAD=1000)
       PARAMETER (MAXANG=200)
C
       LOGICAL ICOUNT,EXIST
       CHARACTER*20 FNAME(2),FORMSTR
       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,KRHOG,MXSPN)
     &  ,RDIS(NSPEED),VLOC(NSPEED)
     &  ,XRAD(MAXRAD),WTRAD(MAXRAD),CENTER(5,MAXSPH)
     &  ,ANGLE(3,MAXANG),DOMEGA(MAXANG)
     &  ,RVECA(3,MX_GRP),PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,ICOUNT(MAX_CON,3)
       DIMENSION NGRID(3),DERIV(3)
       DIMENSION RBAS(3,4)
C
C MODE:  1 ... CALCULATE DENSITY ON A MESH OF POINTS
C        2 ... CALCULATE POTENTIAL ON A MESH OF POINTS
C
C FOR DENSITY EVALUATIONS, USE OLD SCHEME FOR DENSITY CALCULATION
C OTHERWISE, GET DENSITY FROM COUPOT
C
       PRINT '(A)',' '
       PRINT '(A)','POTENTIAL/DENSITY GRIDS'
       CALL GTTIME(TIME1)
       DO 900 MODE=1,2
C
C READ IN NECESSARY INPUT DATA
C
        IF (MODE.EQ.1) THEN
         FNAME(1)='RHOGRID'
        ELSE
         FNAME(1)='POTGRID'
        END IF
        INQUIRE(FILE=FNAME(1),EXIST=EXIST)
C
C DETERMINE IF THERE IS ANYTHING TO DO 
C
        FORMSTR= ' '
        IF (.NOT.EXIST) FORMSTR= ' --> NOTHING TO DO'
        IF (MODE .EQ. 1) PRINT '(2A)','DENSITY GRID  ',FORMSTR
        IF (MODE .EQ. 2) PRINT '(2A)','POTENTIAL GRID',FORMSTR
        IF (.NOT.EXIST) GOTO 900
C
C READ INPUT DATA
C
        OPEN(72,FILE=FNAME(1),FORM='FORMATTED',STATUS='OLD')
        REWIND(72)
C
C CHECK IF THIS FILE IS EMPTY
C IF YES, CREATE A DEFAULT ONE
C
        I=1
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
        write(6,*)(NGRID(J), J=1,3)
        IF (ITYPE .EQ. 1) THEN
         DO I=1,4
          READ(72,*,END=880)(RBAS(J,I), J=1,3)
         END DO
        END IF
        DO I=1,3
         IF (NGRID(I) .LT. 1) THEN
          write(6,*)'POTRHOGRID: NUMBER OF GRID POINTS MUST BE >= 1'
          GOTO 890
         END IF
        END DO
        IF (NGRID(3) .GT. MAX_PTS) THEN
         write(6,*)'POTRHOGRID: MAX_PTS MUST BEAT LEAST: ',NGRID(3)
         write(6,*)'SKIPPING GRID EVALUATION FOR MODE ',MODE
         GOTO 890
        END IF
        NLOOP=NGRID(1)*NGRID(2)
        IF(NSPN.EQ.1) THEN
        FNAME(1)=FNAME(1)(1:3)//'GROUT'
        ELSE
         IF (MODE.EQ.1) THEN
           FNAME(1)=FNAME(1)(1:3)//'TOT'
           FNAME(2)=FNAME(1)(1:3)//'SPN'
         ELSE
           FNAME(1)=FNAME(1)(1:3)//'MAJ'
           FNAME(2)=FNAME(1)(1:3)//'MIN'
         ENDIF
        ENDIF
C
C OPEN OUTPUT FILES, WRITE HEADER
C
        DO IS=1,NSPN
        IUNIT=73+IS
   20   OPEN(IUNIT,FILE=FNAME(IS),FORM=FORMSTR,STATUS='UNKNOWN')
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
         IF(MODE.EQ.1) THEN
             IF(NSPN.EQ.1) WRITE(IUNIT,*)'SCF TOTAL DENSITY (ANG)'
             IF(NSPN.EQ.2) WRITE(IUNIT,*)'SCF DENSITY (ANG)'
         ELSE
             WRITE(IUNIT,2001)'SCF POTENTIAL'
         END IF
         OPEN(77,FILE='XMOL.DAT')
         READ(77,*) NATOM
         READ(77,*) 
         WRITE(IUNIT,'(1X,I10,3F20.12)') NATOM,(RBAS(J,1),J=1,3)
         DO I=1,3
         WRITE(IUNIT,'(1X,I10,3F20.12)') NGRID(I),(RBAS(J,I+1),J=1,3)
         ENDDO
         DO I=1,NATOM
           READ(77,*)IZ, X, Y, Z
           CHR=REAL(IZ)
           WRITE(IUNIT,2002)IZ, CHR, X, Y, Z
         END DO
         CLOSE(77)
c         WRITE(IUNIT,'(4(1X,I10))')(NGRID(J), J=1,3),NGRID(3)
c         IF (ITYPE .EQ. 1) THEN
c          WRITE(IUNIT,'(3(1X,F20.12))')((RBAS(J,I), J=1,3), I=1,4)
c         END IF
        ENDIF
        ENDDO

2001    FORMAT(A25)
2002    FORMAT(I6,4F16.10)
C
C LOOP FOR EACH PILE
C
        DO 850 ILOOP=1,NLOOP
C
C SETUP POINTS
C
c         NMSH=NGRID(3)
         NPTS=NGRID(3)
         IF (ITYPE .EQ. 1) THEN
          IY=MOD(ILOOP-1,NGRID(2))
          IX=(ILOOP-1)/NGRID(2)
          DO IZ=1,NPTS
           I=IZ-1
           RMSH(1,IZ)=RBAS(1,1)+IX*RBAS(1,2)+IY*RBAS(1,3)+I*RBAS(1,4)
           RMSH(2,IZ)=RBAS(2,1)+IX*RBAS(2,2)+IY*RBAS(2,3)+I*RBAS(2,4)
           RMSH(3,IZ)=RBAS(3,1)+IX*RBAS(3,2)+IY*RBAS(3,3)+I*RBAS(3,4)
           WMSH(IZ)=0.0D0
          END DO
         ELSE
          DO IZ=1,NPTS
           READ(72,*,END=870)(RMSH(I,IZ), I=1,3)
           WMSH(IZ)=0.0D0
          END DO
         END IF
C
C NOW: CALCULATE ELECTRONIC COULOMB POTENTIAL AND/OR DENSITY
C DENSITY WILL BE STORED IN RHOG
C
         NGRAD=1
         MODDEN=2
         IF (MODE .EQ. 2) THEN
          MODDEN=1
          IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10
          DO IAT=1,NIDENT
           GAUSS_CUT(IAT)=1.0D30
          END DO
c
#ifdef MPI
c
c MPI: need to send updated mesh data
c
          CALL SENDDATA(102)
c
#endif
c
          CALL COUPOT1
         ELSE
          I1=IGGA(1)
          I2=IGGA(2)
          IGGA(1)=0
          IGGA(2)=0
          CALL DENSOLD(VOL)
          IGGA(1)=I1
          IGGA(2)=I2
         END IF
C
C UPDATE DATA IN RHOG
C
         DO IGRAD=1,NGRAD
          DO IPTS=1,NPTS
           RHOG(IPTS,IGRAD,1)=RHOG(IPTS,IGRAD,1)+RHOG(IPTS,IGRAD,NSPN)
          END DO 
         END DO 
C
C DENSITY GRID
C
C        DENSITY IN  ANG^-3 FOR THE GAUSSIAN CUBE FORMAT
C        
        FACTOR=0.5292D0**3
        IF (MODE .EQ. 1) THEN
        
         IF (NSPN.EQ.1) THEN
           IF (ITYPE .EQ. 1) THEN

            IF (IFORM .EQ. 1) THEN
             WRITE(74)(RHOG(IPTS,1,1), IPTS=1,NPTS)
            ELSE
             WRITE(74,9010)(RHOG(IPTS,1,1)*FACTOR, IPTS=1,NPTS)
            END IF
           ELSE

           DO IPTS=1,NPTS
             WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
             WRITE(74,9010) RHOG(IPTS,1,1)
            END DO
          END IF
         ELSE     !IF SPIN-POLARISED
          DO IS=1,NSPN
           IUNIT=73+IS
           IF (ITYPE .EQ. 1) THEN
            IF (IFORM .EQ. 1) THEN
             WRITE(IUNIT)(RHOG(IPTS,1,1),
     &                 RHOG(IPTS,1,1)-2*RHOG(IPTS,1,NSPN), IPTS=1,NPTS)
            ELSE
            IF(IS.EQ.1) WRITE(IUNIT,9010) (RHOG(IPTS,1,1)
     &               *FACTOR,IPTS=1,NPTS)
            IF(IS.EQ.2) WRITE(IUNIT,9010) ((RHOG(IPTS,1,1)-2*
     &               RHOG(IPTS,1,NSPN)) *FACTOR,IPTS=1,NPTS)
            END IF
           ELSE
            DO IPTS=1,NPTS
             WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
             WRITE(74,9010) RHOG(IPTS,1,1),
     &                 RHOG(IPTS,1,1)-2*RHOG(IPTS,1,NSPN)
            END DO
            END IF
          END DO
         END IF
         GOTO 850
         END IF
C
C THE FOLLOWING PART IS ONLY DONE IF (MODE .EQ. 2)
C CALCULATING KOHN-SHAM POTENTIAL POT 
C
         CALL GETVLXC(MAX_PTS,RHOG,POT,POTIN)
C
C ADD EFIELD POTENTIAL TO LOCAL POTENTIAL
C
         DO IPTS=1,NPTS
          CALL EXTPOT(RMSH(1,IPTS),EXHERE,DERIV)
          EFHERE=EFIELD(1)*RMSH(1,IPTS)+EFIELD(2)*RMSH(2,IPTS)
     &          +EFIELD(3)*RMSH(3,IPTS)
          POTIN(IPTS)=POTIN(IPTS)+EFHERE+EXHERE
         END DO
C
C WRITE OUTPUT: POTIN   CONTAINS LOCAL  
C               COULOMB CONTAINS COULOMB 
C               POT     CONTAINS EXCHANGE-CORRELATION 
C
         IF (NSPN.EQ.1) THEN
          IF (ITYPE .EQ. 1) THEN
           IF (IFORM .EQ. 1) THEN
            WRITE(74)(POTIN(IPTS),COULOMB(IPTS),POT(IPTS),
     &                IPTS=1,NPTS)
           ELSE
            WRITE(74,9010)(POTIN(IPTS),COULOMB(IPTS),
     &                     POT(IPTS), IPTS=1,NPTS)
           END IF
          ELSE
           DO IPTS=1,NPTS
            WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
            WRITE(74,9010) POTIN(IPTS),COULOMB(IPTS),POT(IPTS)
           END DO
          END IF
         ELSE
          IF (ITYPE .EQ. 1) THEN
           IF (IFORM .EQ. 1) THEN
            WRITE(74)(POTIN(IPTS),COULOMB(IPTS),
     &                POT(IPTS),POT(IPTS+NPTS), IPTS=1,NPTS)
           ELSE
            WRITE(74,9010)(POTIN(IPTS),COULOMB(IPTS),
     &                     POT(IPTS),POT(IPTS+NPTS), IPTS=1,NPTS)
           END IF
          ELSE
           DO IPTS=1,NPTS
            WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
            WRITE(74,9010) POTIN(IPTS),COULOMB(IPTS),
     &                     POT(IPTS),POT(IPTS+NPTS)
           END DO
          END IF
         END IF
  850   CONTINUE
        CLOSE(74)
        GOTO 890
C
C ERROR HANDLING
C        
  870   CLOSE(74)
  880   write(6,*)'ERROR IN INPUT FILE, SKIPPING MODE ',MODE
  890   CLOSE(72)
  900  CONTINUE
       RETURN
 9010  FORMAT(3(1X,E20.12))
 9020  FORMAT(3(1X,F20.12))
       END
