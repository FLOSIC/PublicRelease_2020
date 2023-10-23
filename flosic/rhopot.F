C UTEP Electronic Structure Lab (2020)
C
       SUBROUTINE RHOPOT
C
C POTRHOGRID VERSION DIRK POREZAG JUNE 1998. 
C CALCULATION OF CHARGES WITHIN A SPHERE moved into atomsph JK 3/99
C * DENSITY AND SPIN DENSITY ON A GRID OF POINTS  -- TB 04/03
C * ELF ON A GRID OF POINTS  -- TB 05/03
C * POTENTIAL ON A GRID -- TB 08/03

       use common2,only : RCNT, NCNT, IGGA, ISPN, NSPN, EFIELD
       use common3,only : RMAT
       use common7,only : MODDEN
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:58 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       REAL*8 :: SYMBOL 
       INTEGER :: I, I1, I2, IBEG, ICNT, IELF, IFORM, ILOOP, IPTS, IS,
     & ITYPE, IUNIT, IX, IY, IZ, J, JS, MODE, NATOM, NGRAD, NGRID,
     & NLOOP, NMODE, NMSH, NSTORE
       REAL*8 :: SYMBOL , CHR, COULOMB, D, D1, DENTOT, DERIV, DH,
     & EFHERE, EXHERE, FACT1, FACT4, FACTOR, PHIG, PI, POT, POTIN,
     & RBAS, RHOG, RHOSQ, TIME1, TMP, VOL, X, Y, Z
       SAVE
       PARAMETER (MAXSPH=500)
       PARAMETER (MAXRAD=1000)
       PARAMETER (MAXANG=200)
C
       LOGICAL ICOUNT,EXIST,ELF
       LOGICAL DENSTY,POTEN
       CHARACTER*20 FNAME(6),FORMSTR
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       DIMENSION NGRID(3),DERIV(3)
       DIMENSION RBAS(3,4),TMP(3)
       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,KRHOG,MXSPN)
     &  ,PHIG(MAX_PTS,2)
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
       DENSTY=.FALSE.
       POTEN =.FALSE.
       NSTORE=NMSH
C
C READ IN NECESSARY INPUT DATA
C
         FORMSTR= ' '

         FNAME(1)='RHOGRID'
         INQUIRE(FILE=FNAME(1),EXIST=EXIST)
         IF (EXIST) DENSTY=.TRUE.
         IF (.NOT.EXIST) FORMSTR= ' --> NOTHING TO DO'
         PRINT '(2A)','DENSITY GRID  ',FORMSTR

         FNAME(2)='POTGRID'
         INQUIRE(FILE=FNAME(2),EXIST=EXIST)
         IF (EXIST) POTEN=.TRUE.
         IF (.NOT.EXIST) THEN
           FORMSTR= ' --> NOTHING TO DO'
           PRINT '(2A)','POTENTIAL GRID  ',FORMSTR
         END IF
C
C DETERMINE IF THERE IS ANYTHING TO DO 
C

        IF ((.NOT.DENSTY).AND.(.NOT.POTEN)) GOTO 900
C
C READ INPUT DATA
C
        IF (.NOT.DENSTY)FNAME(1)="POTGRID"
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
          write(6,*)'RHOGRID: NUMBER OF GRID POINTS MUST BE >= 1'
          GOTO 890
         END IF
        END DO
        IF (NGRID(3) .GT. MAX_PTS) THEN
         write(6,*)'RHOGRID: MAX_PTS MUST BEAT LEAST: ',NGRID(3)
         write(6,*)'SKIPPING GRID EVALUATION '
         GOTO 890
        END IF
        NLOOP=NGRID(1)*NGRID(2)
        IF(DENSTY) THEN
          IF(NSPN.EQ.1) THEN
             FNAME(1)=FNAME(1)(1:3)//'GROUT'
          ELSE
             FNAME(1)=FNAME(1)(1:3)//'TOT'
             FNAME(2)=FNAME(1)(1:3)//'SPN'
          ENDIF
        END IF
        IF(POTEN) THEN
          FNAME(5)='POT'
          IF (NSPN.GT.1) THEN
              FNAME(5)=FNAME(5)(1:3)//'UP'
              FNAME(6)=FNAME(5)(1:3)//'DN'
          END IF
        ENDIF
        IF(ELF) THEN
              FNAME(3)='ELF'
         IF (NSPN.GT.1) THEN
              FNAME(3)=FNAME(3)(1:3)//'UP'
              FNAME(4)=FNAME(3)(1:3)//'DN'
         END IF
        END IF
               
   
C
C OPEN OUTPUT FILES, WRITE HEADER
C
        IF(POTEN) NMODE=2

        DO 700 MODE=1,NMODE
          IF(MODE.EQ.1) THEN
            IF (.NOT.DENSTY) GO TO 700 
            IBEG=73
            JS=0
          ELSE
            IBEG=75
            JS = 4
          ENDIF
        DO IS=1,NSPN
        IUNIT=IBEG+IS
   20   OPEN(IUNIT,FILE=FNAME(JS+IS),FORM=FORMSTR,STATUS='UNKNOWN')
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
         IF (MODE.EQ.1) THEN
             WRITE(IUNIT,*)'CLUSTER OUTPUT'
             IF(NSPN.EQ.1) WRITE(IUNIT,*)'SCF TOTAL DENSITY (ANG)'
             IF(NSPN.EQ.2) WRITE(IUNIT,*)'SCF DENSITY (ANG)'
         ELSE
             WRITE(IUNIT,*)'CLUSTER OUTPUT'
             IF(NSPN.EQ.1) WRITE(IUNIT,*)'SCF TOTAL POTENTIAL (ANG)'
             IF(NSPN.EQ.2) WRITE(IUNIT,*)'SCF POTENTIAL (ANG)'
         ENDIF
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
 700    CONTINUE
        
       IF (ELF) THEN
         DO IS=1,NSPN
          IELF=80+IS
          OPEN(IELF,FILE=FNAME(2+IS),FORM=FORMSTR,STATUS='UNKNOWN')
          WRITE(IELF,*)'CLUSTER OUTPUT'
          WRITE(IELF,*)'ELECTRON LOCALIZTION FUNCTION (ANG)'
           OPEN(87,FILE='XMOL.DAT')
           REWIND(87)
           READ(87,*) NATOM
           READ(87,*) 
           WRITE(IELF,'(1X,I10,3F20.12)') NATOM,(RBAS(J,1),J=1,3)
           DO I=1,3
            WRITE(IELF,'(1X,I10,3F20.12)') NGRID(I),(RBAS(J,I+1),J=1,3)
           ENDDO
           DO I=1,NATOM
             READ(87,*)IZ, X, Y, Z
             CHR=REAL(IZ)
             WRITE(IELF,2002)IZ, CHR, X, Y, Z
           END DO
           CLOSE(87)
          END DO
         END IF

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
         MODDEN=2
c
c
c MPI: need to send updated mesh data
c
c
          I1=IGGA(1)
          I2=IGGA(2)
          IF (.NOT.ELF) THEN
             IGGA(1)=0
             IGGA(2)=0
          END IF
          CALL DENSOLD(VOL)
          IGGA(1)=I1
          IGGA(2)=I2

        IF (POTEN) THEN
#ifdef MPI 
        CALL SENDDATA(102)
#endif
        CALL COUPOT1
        ENDIF

C
C  CURRENTLY, ELF is written in RHOG(IPTS,10,1)
C
        FACTOR=0.5292D0**3
        FACT4=0.5292D0**4
        FACT1=0.5292D0
         DO IS=1,NSPN
          DO IPTS=1,NMSH
           RHOG(IPTS,10,IS)=0.0D0
          END DO 
         END DO
        
       IF (ELF) THEN
        DO IPTS=1,NMSH
         DO IS=1,NSPN    
          D=RHOG(IPTS,1,IS)
          IF (D.GT.1.0D-20) THEN
           RHOSQ=RHOG(IPTS,2,IS)**2+RHOG(IPTS,3,IS)**2
     &                                  +RHOG(IPTS,4,IS)**2
           RHOSQ=RHOSQ
           D1=PHIG(IPTS,IS)-0.25D0*RHOSQ/RHOG(IPTS,1,IS)
           DH=0.6d0*(3.0d0*PI**2)**0.66666*RHOG(IPTS,1,IS)**1.6667
           RHOG(IPTS,10,IS)=1.0d0/(1.0d0+(D1/DH)**2)
          END IF
         END DO 
        END DO 
       END IF
C
C
C UPDATE DATA IN RHOG
C
c         DO IGRAD=1,NGRAD
          DO IPTS=1,NMSH
           RHOG(IPTS,1,1)=RHOG(IPTS,1,1)+RHOG(IPTS,1,NSPN)
          END DO 
          DO IPTS=1,NMSH
            DENTOT=DENTOT+RHOG(IPTS,1,1)*WMSH(IPTS)
          END DO

c         END DO 
C DENSITY GRID
C
C        DENSITY IN  ANG^-3 FOR THE GAUSSIAN CUBE FORMAT
C        
         IF (DENSTY) THEN
        
         IF (NSPN.EQ.1) THEN
           IF (ITYPE .EQ. 1) THEN

            IF (IFORM .EQ. 1) THEN
             WRITE(74)(RHOG(IPTS,1,1), IPTS=1,NMSH)
            ELSE
             WRITE(74,9010)(RHOG(IPTS,1,1)*FACTOR, IPTS=1,NMSH)
            END IF
           ELSE

           DO IPTS=1,NMSH
             WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
             WRITE(74,9010) RHOG(IPTS,1,1)
            END DO
          END IF
         ELSE     !IF SPIN-POLARISED
          DO IS=1,NSPN
           IUNIT=73+IS
           IELF =80+IS
           IF (ITYPE .EQ. 1) THEN
            IF (IFORM .EQ. 1) THEN
             WRITE(IUNIT)(RHOG(IPTS,1,1),
     &                 RHOG(IPTS,1,1)-2*RHOG(IPTS,1,NSPN), IPTS=1,NMSH)
            ELSE
            IF(IS.EQ.1) WRITE(IUNIT,9010) (RHOG(IPTS,1,1)
     &               *FACTOR,IPTS=1,NMSH)
            IF(ELF)WRITE(IELF,9010) (RHOG(IPTS,10,IS),IPTS=1,NMSH)
            IF(IS.EQ.2) WRITE(IUNIT,9010) ((RHOG(IPTS,1,1)-2*
     &               RHOG(IPTS,1,NSPN)) *FACTOR,IPTS=1,NMSH)
            END IF
           ELSE
            DO IPTS=1,NMSH
             WRITE(IUNIT,9020)(RMSH(I,IPTS), I=1,3)
             WRITE(IUNIT,9010) RHOG(IPTS,1,1),
     &                 RHOG(IPTS,1,1)-2*RHOG(IPTS,1,NSPN)
            END DO
            END IF
          END DO
         END IF
        END IF

        IF (POTEN) THEN
          WRITE(6,*) 'WRITING POTENTIAL'
C THE FOLLOWING PART IS ONLY DONE IF (MODE .EQ. 2)
C CALCULATING KOHN-SHAM POTENTIAL POT 
C
         CALL GETVLXC(MAX_PTS,RHOG,POT,POTIN)
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
  870   CLOSE(74)
        IF (NSPN.EQ.2) CLOSE(75)
  880   write(6,*)'ERROR IN INPUT FILE, SKIPPING MODE ',MODE
  890   CLOSE(72)
  900  CONTINUE
       NMSH=NSTORE
       RETURN
 9010  FORMAT(3(1X,E20.12))
 9020  FORMAT(3(1X,F20.12))
       END
