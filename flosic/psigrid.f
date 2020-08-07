C UTEP Electronic Structure Lab (2020)
C ****************************************************************
C
       SUBROUTINE PSIGRID
C
C POTRHOGRID VERSION DIRK POREZAG JUNE 1998. 
C CALCULATION OF CHARGES WITHIN A SPHERE moved into atomsph JK 3/99
C * DENSITY AND SPIN DENSITY ON A GRID OF POINTS  -- TB 04/03
C * ELF ON A GRID OF POINTS  -- TB 05/03
C
       use common2,only : RCNT, NCNT, IGGA, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI
       use common7,only : MODDEN
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:59 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I1, I2, ICNT, IELF, IFORM, ILOOP, IPTS, IS, ITYPE,
     & IUNIT, IX, IY, IZ, J, MODE, NATOM, NGRAD, NGRID, NLOOP, NMSH,
     & NSTORE
       REAL*8 :: SYMBOL , CHR, COULOMB, D, D1, DERIV, DH, FACT1, FACT4,
     & FACTOR, PHIG, PI, RBAS, RHOG, RHOSQ, TIME1, VOL, X, Y, Z
       SAVE
       PARAMETER (MAXSPH=500)
       PARAMETER (MAXRAD=1000)
       PARAMETER (MAXANG=200)
C
       LOGICAL ICOUNT,EXIST,ELF
       CHARACTER*20 FNAME(4),FORMSTR
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,KRHOG,MXSPN)
     &  ,PHIG(MAX_PTS,2)
       DIMENSION NGRID(3),DERIV(3)
       DIMENSION RBAS(3,4)
C
C
C FOR DENSITY EVALUATIONS, USE OLD SCHEME FOR DENSITY CALCULATION
C OTHERWISE, GET DENSITY FROM COUPOT
C
       PI=4.0D0*DATAN(1.0D0)
       PRINT '(A)',' '
       PRINT '(A)','PSIGRID GRIDS'
       CALL GTTIME(TIME1)
       NSTORE=NMSH
C
C READ IN NECESSARY INPUT DATA
C
         FNAME(1)='PSIGRID'
         INQUIRE(FILE=FNAME(1),EXIST=EXIST)
C
C DETERMINE IF THERE IS ANYTHING TO DO 
C
        FORMSTR= ' '
        IF (.NOT.EXIST) FORMSTR= ' --> NOTHING TO DO'
        PRINT '(2A)','DENSITY GRID  ',FORMSTR
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
         PRINT *,'PSIGRID: MAX_PTS MUST BEAT LEAST: ',NGRID(3)
         PRINT *,'SKIPPING GRID EVALUATION '
         GOTO 890
        END IF
        NLOOP=NGRID(1)*NGRID(2)
        IF(NSPN.EQ.1) THEN
        FNAME(1)=FNAME(1)(1:3)//'GROUT'
        ELSE
           FNAME(1)=FNAME(1)(1:3)//'TOT'
           FNAME(2)=FNAME(1)(1:3)//'SPN'
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
             IF(NSPN.EQ.1) WRITE(IUNIT,*)'SCF TOTAL DENSITY (ANG)'
             IF(NSPN.EQ.2) WRITE(IUNIT,*)'SCF DENSITY (ANG)'
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
         NMSH=NGRID(3)
         IF (ITYPE .EQ. 1) THEN
          IY=MOD(ILOOP-1,NGRID(2))
          IX=(ILOOP-1)/NGRID(2)
          DO IZ=1,NMSH
           I=IZ-1
           RMSH(1,IZ)=RBAS(1,1)+IX*RBAS(1,2)+IY*RBAS(1,3)+I*RBAS(1,4)
           RMSH(2,IZ)=RBAS(2,1)+IX*RBAS(2,2)+IY*RBAS(2,3)+I*RBAS(2,4)
           RMSH(3,IZ)=RBAS(3,1)+IX*RBAS(3,2)+IY*RBAS(3,3)+I*RBAS(3,4)
           WMSH(IZ)=0.0D0
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
          CALL GETPSI(VOL)
          IGGA(1)=I1
          IGGA(2)=I2

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
        
       DO IPTS=1,NMSH
        DO IS=1,NSPN    
         D=RHOG(IPTS,1,IS)
         IF (D.GT.1.0D-20) THEN
          RHOSQ=RHOG(IPTS,2,IS)**2+RHOG(IPTS,3,IS)**2+RHOG(IPTS,4,IS)**2
          RHOSQ=RHOSQ
          D1=PHIG(IPTS,IS)-0.25D0*RHOSQ/RHOG(IPTS,1,IS)
          DH=0.6d0*(3.0d0*PI**2)**0.66666*RHOG(IPTS,1,IS)**1.6667
          RHOG(IPTS,10,IS)=1.0d0/(1.0d0+(D1/DH)**2)
         END IF
        END DO 
       END DO 
C
C
C UPDATE DATA IN RHOG
C
c         DO IGRAD=1,NGRAD
          DO IPTS=1,NMSH
           RHOG(IPTS,1,1)=RHOG(IPTS,1,1)+RHOG(IPTS,1,NSPN)
          END DO 
c         END DO 
C DENSITY GRID
C
C        DENSITY IN  ANG^-3 FOR THE GAUSSIAN CUBE FORMAT
C        
c        IF (MODE .EQ. 1) THEN
        
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
            WRITE(IELF,9010) (RHOG(IPTS,10,IS),IPTS=1,NMSH)
            IF(IS.EQ.2) WRITE(IUNIT,9010) ((RHOG(IPTS,1,1)-2*
     &               RHOG(IPTS,1,NSPN)) *FACTOR,IPTS=1,NMSH)
            END IF
           ELSE
            DO IPTS=1,NMSH
             WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
             WRITE(74,9010) RHOG(IPTS,1,1),
     &                 RHOG(IPTS,1,1)-2*RHOG(IPTS,1,NSPN)
            END DO
            END IF
          END DO
         END IF
C
  850   CONTINUE
        CLOSE(74)
        IF (NSPN.EQ.2) CLOSE(75)
        CLOSE(81)
        IF (NSPN.EQ.2) CLOSE(82)
        GOTO 890
C
C ERROR HANDLING
C        
  870   CLOSE(74)
        IF (NSPN.EQ.2) CLOSE(75)
  880   PRINT *,'ERROR IN INPUT FILE, SKIPPING MODE ',MODE
  890   CLOSE(72)
  900  CONTINUE
       NMSH=NSTORE
       RETURN
 9010  FORMAT(3(1X,E20.12))
 9020  FORMAT(3(1X,F20.12))
       END
