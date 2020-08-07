C UTEP Electronic Structure Lab (2020)
       SUBROUTINE DIR 
C
C CALCULATES THE QUANTUM DIPOLE MATRICES  -- TB 05/03
C
       use mesh1,only : wmsh,rmsh,nmsh

       use common2,only : RIDT, RCNT, NCNT, N_CON, LSYMMAX, N_POS,
     &   NFNCT, NSPN, DIPOLE
       use common3,only : RMAT, NGRP
       use common5,only : PSI, NWF, EFERMI
       use common8,only : REP, N_REP, NDMREP, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:42 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NMAX, I, I_POS, IBS, ICNT, ICON, IERR, IFNCT, IGP,
     & II, ILOC, INDEX, INFOWF, IOFS, IPILE, IPTS, IPV, IRP, ISHDUM,
     & ISHDUMMY, ISHELLA, ISIZE, ISP, ITOTWF, IUNIT, IWF, IX, IY, IZ,
     & J, J_POS, JOFS, JPTS, K, L_NUC, LI, LMAX1, LPV, M_NUC, MPTS, MU,
     & MXSPEC, NATOM, NDIM, NFILE, NGRID, NOFS, NPILE, NPV, NUNSYM,
     & NWAVF
       REAL*8 :: SYMBOL , CHARGE, CHR, CSP, CSS, CST, DIM, DR, EVAL,
     & EVC, EVEC, FACT, FACTOR, FCGRP, GRAD, H, HA2EV, O, P, PHI, PS,
     & PSIG, PTS, Q, RBAS, RGRID, ROT, ROTI, RR, RVECA, RX, RY, RZ,
     & SC1, SNP, SNS, SNT, SPDIP, SPTOT, SS, SUM, SUMCHARGE, TEMP,
     & THETA, TIME1, TIME2, V, VX, VY, VZ, X, Y, Z
       SAVE
       PARAMETER (MXSPEC=10000)
       PARAMETER (NMAX=MPBLOCK)
C
       LOGICAL IUPDAT,ICOUNT,EXIST,LMKFIL
       CHARACTER*4 FILE
       CHARACTER*6 FILENAME(10)
       COMMON/TMP2/PSIG(NMAX,MAX_OCC)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       COMMON/TMP1/H(MAX_OCC,3)
     &  ,SPDIP(MXSPEC),SPTOT(MXSPEC),RVECA(3,MX_GRP)
     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,ICOUNT(MAX_CON,3),RGRID(3,MPBLOCK)
       DIMENSION O(MAX_OCC,MAX_OCC,3)
       DIMENSION EVC(10,10,3),ROT(3,3),SS(3),ROTI(3,3)
       DIMENSION ISIZE(3),ISP(10),IRP(10),IBS(10),INFOWF(4,10)
       DIMENSION P(NMAX,3),Q(NMAX,3),V(NMAX),CHARGE(MAX_OCC)
       DIMENSION RBAS(3,4),NGRID(3),EVEC(NMAX),dr(3)
       DIMENSION SUMCHARGE(3)
       DATA ISIZE/1,3,6/
       DATA HA2EV/27.2116D0/
       DATA TEMP/1.0D-4/
       DATA FACT/0.148203857088/
       write(6,*)'NUMBER OF MESH POINTS:',NMSH
C
C RETURN IF INPUT FILE DOES NOT EXIST
C
       PRINT '(A)','CALCULATING QUANTUM DIPOLE MATRICES'
       INQUIRE(FILE='DIR',EXIST=EXIST)
       IF (.NOT.EXIST) THEN
        PRINT '(2A)','QDIPOLE: FILE QDIP DOES NOT EXIST ',
     &               '--> NOTHING TO DO'
        RETURN
       END IF
C
C CREATE A STANDARD INPUT FILE IF THE CURRENT INPUT FILE IS EMPTY
C
       CALL GTTIME(TIME1)


C  GET THE IDENTITY OF HOMO AND ITS DEGENERACY
C
        INDEX=2
        NWAVF=1
C       The following two lines should be removed for a general case.
c        EFERMI(1)=-0.050516704375665364
c        EFERMI(2)=-0.050516704375665364
        CALL FINDSTATE(INDEX,NWAVF,ISP,IRP,IBS,IERR)
         DO I=1,NWAVF
           WRITE(6,*)'STATE :', ISP(I),IRP(I),IBS(I)
           DIM=NDMREP(IRP(I))
          IF ((ISP(I) .LT. 1) .OR. (ISP(I) .GT. NSPN)) THEN
           PRINT *,'WFGRID: SPIN FOR STATE ',I,' IS INVALID'
           GOTO 900
          END IF
          IF ((IRP(I).LT. 1) .OR. (IRP(I).GT. N_REP)) THEN
           PRINT *,'WFGRID: REPRESENTATION FOR STATE ',I,' IS INVALID'
           GOTO 900
          END IF
          IF ((IBS(I) .LT. 1) .OR. (IBS(I) .GT. NS_TOT(IRP(I)))) THEN
           PRINT *,'WFGRID: INDEX FOR STATE ',I,' IS INVALID'
           GOTO 900
          END IF
          NDIM=NDMREP(IRP(I))
          NUNSYM=NUNSYM+NDIM
          IF (NUNSYM .GT. MAX_OCC) THEN
           PRINT *,'WFGRID: NUMBER OF STATES IS TOO LARGE'
           PRINT *,'        INCREASE MAX_OCC TO AT LEAST: ',NUNSYM
           GOTO 900
          END IF
           INFOWF(1,I)=ISP(I)
           INFOWF(2,I)=IRP(I)
           INFOWF(3,I)=IBS(I)
           INFOWF(4,I)=NDIM
          END DO
          ITOTWF=INFOWF(4,1) 
          WRITE(6,*) ITOTWF

C ZERO H (CONTAINS THE DIPOLE MATRIX ELEMENTS)
C
       DO IX=1,3
        DO IWF=1,MAX_OCC
          H(IWF,IX)=0.0D0
        END DO
       END DO
C
C CALCULATE DIPOLE MATRIX ELEMENTS BY MESH INTEGRATION
C
       DO I=1,ITOTWF
       CHARGE(I)=0.0D0
       END DO

       NPILE=NMSH/NMAX
       FCGRP=1.0D0/NGRP
       write(6,*) 'QDIPOLE: ', NMSH, NMAX, NGRP, FCGRP
       DO 850 IPILE=0,NPILE
        NOFS=IPILE*NMAX
        MPTS=MIN(NMAX,NMSH-NOFS)
        DO IPTS=1,MPTS
         Q(IPTS,1)=RMSH(1,IPTS+NOFS)
         Q(IPTS,2)=RMSH(2,IPTS+NOFS)
         Q(IPTS,3)=RMSH(3,IPTS+NOFS)
         V(  IPTS)=WMSH(IPTS+NOFS)*FCGRP
c        PRINT 432,IPTS+NOFS,(Q(IPTS,J),J=1,3)
 432    FORMAT(' MESH TEST:',I5,3G15.6)               
        END DO
        DO 800 IGP=1,NGRP
         DO J=1,3
          DO IPTS=1,MPTS
           P(IPTS,J)=0.0D0
          END DO
          DO K=1,3
           DO IPTS=1,MPTS
            P(IPTS,J)=P(IPTS,J)+RMAT(K,J,IGP)*Q(IPTS,K)
           END DO
          END DO
         END DO
          
C
C INITIALIZE PSIG 
C
         DO IWF=1,NWF
          DO IPTS=1,MPTS
           PSIG(IPTS,IWF)=0.0D0
          END DO
         END DO  
         ISHELLA=0
         DO 86 IFNCT=1,NFNCT
          LMAX1=LSYMMAX(IFNCT)+1
          DO 84 I_POS=1,N_POS(IFNCT)
           ISHELLA=ISHELLA+1
           CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUMMY)
           DO 82 J_POS=1,M_NUC
             CALL WFRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &           RVECA,L_NUC,1,INFOWF,NWAVF,ITOTWF)
            IF(L_NUC.NE.M_NUC)THEN
             PRINT *,'QDIPOLE: PROBLEM IN UNRAVEL'
             CALL STOPIT
            END IF
            DO 80 JPTS=1,MPTS,NSPEED
             NPV=MIN(NSPEED,MPTS-JPTS+1)
             IPTS=JPTS-1
             DO LPV=1,NPV
              PTS(LPV,1)=P(IPTS+LPV,1)-RVECA(1,J_POS)
              PTS(LPV,2)=P(IPTS+LPV,2)-RVECA(2,J_POS)
              PTS(LPV,3)=P(IPTS+LPV,3)-RVECA(3,J_POS)
             END DO
             CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)

             IF (IUPDAT) THEN
              DO IWF=1,ITOTWF
              ILOC=0
              DO 78 LI=1,LMAX1
               DO MU=1,ISIZE(LI)
                DO ICON=1,N_CON(LI,IFNCT)
                 ILOC=ILOC+1
                   FACTOR=PSI(ILOC,IWF,1)
                   DO LPV=1,NPV
                    PSIG(IPTS+LPV,IWF)=PSIG(IPTS+LPV,IWF)
     &              +FACTOR*GRAD(LPV,1,MU,ICON,LI)
                   END DO  
                  END DO  
                END DO  
   78         CONTINUE
             END DO  
            END IF

   80       CONTINUE
   82      CONTINUE
   84     CONTINUE
   86    CONTINUE
C
C UPDATE CHARGE AND DIPOLE MATRICES
C
        DO IWF=1,ITOTWF
         DO IPTS=1,MPTS
************************************************************************
c       CHARGE(IWF)=CHARGE(IWF)+V(IPTS)*P(IPTS,1)*P(IPTS,3)*P(IPTS,2)
         CHARGE(IWF)=CHARGE(IWF)+V(IPTS)*PSIG(IPTS,IWF)*PSIG(IPTS,IWF)
         END DO
        END DO

         DO IX=1,3
          DO IWF=1,ITOTWF
            DO IPTS=1,MPTS
             H(IWF,IX)=H(IWF,IX)
     &       +PSIG(IPTS,IWF)*V(IPTS)*P(IPTS,IX)
            END DO
c           END DO
          END DO
         END DO
  800   CONTINUE
  850  CONTINUE
       

       DO IWF=1,ITOTWF
        WRITE(6,*) 'IWF = ',IWF,' CHARGE =', CHARGE(IWF)
       END DO
C      WRITE OUT THE MATRICES

       OPEN(UNIT=73,FILE='DIRWF')
        
       WRITE(73,*) 'HOMO IDENTIFICATION :'
       WRITE(73,*) 'SPIN = ', ISP(1), ' REPRESENTATION : ', IRP(1)
     &    , ' BASIS = ', IBS(1), ' DEGENERACY = ',INFOWF(4,1)
       WRITE(73,*)
       WRITE(73,*) 'ORIENTATION OF THE WAVEFUNCTIONS : '
       WRITE(73,*)

       DO IX=1,3
        SUM=0.0D0
         DO IWF=1,ITOTWF
           SUM=SUM+H(IWF,IX)**2
         END DO
         SUM=SQRT(SUM)
         DO IWF=1,ITOTWF
           H(IWF,IX)=H(IWF,IX)/SUM
         END DO
       END DO

       DO IWF=1,3
        WRITE(73,*) 
        WRITE(73,*) 'IWF = ', IWF
          WRITE(73,*) 
          WRITE(73,*) 'ALPHA (IWF,IX=1,3) :'
          WRITE(73,*) 
           WRITE(73,*) (H(IWF,IX),IX=1,3)
       ENDDO
C
C CALCULATE DOT PRODUCTS:
C
       WRITE(73,*)
       WRITE (73,*),'DOT PRODUCTS'
          DO I=1,3
           DO J=1,3
             SUM=0.0D0
              DO K=1,ITOTWF
                SUM=SUM+H(I,K)*H(J,K)             
              END DO
             WRITE(73,*)' I, J ', I,J, ' PRODUCT= ', SUM
            END DO
          END DO
C
        
       DO I=1,NDH
        DO J=1,NDH
         OVER(I,J)=0.0d0
         HAM(I,J)=0.0d0
        END DO
       END DO

       DO I=1,3
        EVAL(I)=0.0d0
        DO J=1, 3
         SUM=0.0D0
         DO K=1, ITOTWF
          SUM=SUM+H(K,I)*H(K,J)
         END DO
         OVER(I,J)=SUM
        END DO
         WRITE(6,*) 'OVERLAP = ', (OVER(I,J),J=1,3)
       END DO

      CALL LOWDEN(NDH,ITOTWF,OVER,HAM,EVAL,SC1)

        DO IX=1,3
            WRITE(73,*)'HAM IX = ', IX,(HAM(I,IX),I=1,3)
c           DO I=1,3
c            H(I,IX)=HAM(I,IX)
c           END DO
        END DO

C      GET THE ROTATION MATRIX
        RX=H(1,1)
        RY=H(1,2)
        RZ=H(1,3)
        RR=SQRT(RX*RX+RY*RY+RZ*RZ)
        THETA=ACOS(RZ/RR)
        PHI=ATAN(RY/RX)
        SNP=SIN(PHI)
        SNT=SIN(THETA)
        CSP=COS(PHI)
        CST=COS(THETA)
        ROT(1,1)=CSP*CST
        ROT(1,2)=SNP*CST
        ROT(1,3)=-SNT
        ROT(2,1)=-SNP
        ROT(2,2)=CSP
        ROT(2,3)=0.0
        ROT(3,1)=CSP*SNT
        ROT(3,2)=SNP*SNT
        ROT(3,3)=CST

        VX=ROT(1,1)*RX+ROT(1,2)*RY+ROT(1,3)*RZ
        VY=ROT(2,1)*RX+ROT(2,2)*RY+ROT(2,3)*RZ
        VZ=ROT(3,1)*RX+ROT(3,2)*RY+ROT(3,3)*RZ
        WRITE(73,*)'VX,VY,VZ', VX,VY,VZ
        RX=H(2,1)
        RY=H(2,2)
        RZ=H(2,3)
        VX=ROT(1,1)*RX+ROT(1,2)*RY+ROT(1,3)*RZ
        VY=ROT(2,1)*RX+ROT(2,2)*RY+ROT(2,3)*RZ
        VZ=ROT(3,1)*RX+ROT(3,2)*RY+ROT(3,3)*RZ
        WRITE(73,*)'VX,VY,VZ', VX,VY,VZ
        WRITE(73,*)'RX,RY,RZ', RX,RY,RZ
        
        do i=1,3
         write(73,*) 'rot  i', i, (ROT(i,j),j=1,3)
        end do

        PS=DATAN(VY/VX)
        WRITE(73,*)'PSI= ', PS
        CSS=COS(PS)
        SNS=SIN(PS)
        ROT(1,1)=CSS*CSP*CST-SNP*SNS
        ROT(1,2)=CSS*SNP*CST+SNS*CSP
        ROT(1,3)=-SNT*CSS
        ROT(2,1)=-SNS*CSP*CST-SNP*CSS
        ROT(2,2)=-SNS*SNP*CST+CSS*CSP
        ROT(2,3)=SNT*SNS
        ROT(3,1)=CSP*SNT
        ROT(3,2)=SNP*SNT
        ROT(3,3)=CST



        ROTI(1,1)=ROT(2,2)*ROT(3,3) -ROT(2,3)*ROT(3,2)
        ROTI(1,2)=-(ROT(1,2)*ROT(3,3)-ROT(3,2)*ROT(1,3))
        ROTI(1,3)=ROT(1,2)*ROT(2,3)-ROT(2,2)*ROT(1,3)
        ROTI(2,1)=-(ROT(2,1)*ROT(3,3)-ROT(3,1)*ROT(2,3))
        ROTI(2,2)=ROT(1,1)*ROT(3,3)-ROT(1,3)*ROT(3,1)
        ROTI(2,3)=-(ROT(1,1)*ROT(2,3)-ROT(2,1)*ROT(1,3))
        ROTI(3,1)=ROT(2,1)*ROT(3,2)-ROT(3,1)*ROT(2,2)
        ROTI(3,2)=-(ROT(1,1)*ROT(3,2)-ROT(3,1)*ROT(1,2))
        ROTI(3,3)=ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1)

        do IWF=1,3
        do i=1,3
          dr(i)=0.0d0
        end do
        do i=1,3
         do j=1,3
            dr(i)=dr(i)+H(IWF,j)*ROT(i,j)
         end do
        end do
         write(73,1000) IWF,  (DR(j),j=1,3)
        end do

        do IWF=1,3
        do i=1,3
          dr(i)=0.0d0
        end do
        do i=1,3
         do j=1,3
            dr(i)=dr(i)+H(IWF,j)*ROTI(i,j)
         end do
        end do
         write(73,1001) IWF,  (DR(j),j=1,3)
        end do
1000    format('ROT  ', I6, '=',3F16.6)
1001    format('ROTI ', I6, '=',3F16.6)
        do i=1,3
         do j=1,3
           over(i,j)=0.0d0
         end do
        end do
        do i=1,3
         do j=1,3
           do k=1,3
           over(i,j)=over(i,j)+ROTI(i,k)*ROt(k,j)
           end do
         end do
        end do

        do i=1,3
         write(73,*) 'test  i', i, (over(i,j),j=1,3)
        end do
        do i=1,3
         write(73,*) 'rot  i', i, (ROT(i,j),j=1,3)
        end do
        do i=1,3
         write(73,*) 'rotI  i', i, (ROTI(i,j),j=1,3)
        end do
C
C   NOW PLOT THE WAVEFUNCTIONS 
C

C
C  MESH FOR PLOTTING
C

       DO I=1,4
        DO J=1,3
           RBAS(I,J)=0.0D0
        END DO
       END DO

       DO I=1,3
          RBAS(I,1)=1.0D30
          RBAS(I,2)=-1.0D30
       END DO
       DO ICNT=1,NCNT
          DO I=1,3
            RBAS(I,1)=MIN(RBAS(I,1),RCNT(I,ICNT))
            RBAS(I,2)=MAX(RBAS(I,2),RCNT(I,ICNT))
          END DO
       END DO
       DO I=1,3
          RBAS(I,1)=RBAS(I,1)-5.0D0
          RBAS(I,2)=RBAS(I,2)+5.0D0
          NGRID(I)=(RBAS(I,2)-RBAS(I,1))/0.5D0+2
          RBAS(I,2)=(RBAS(I,2)-RBAS(I,1))/(NGRID(I)-1)
       END DO
       DO I=1,3
         RBAS(I,I+1)=RBAS(I,2)
       END DO
         RBAS(2,2)=0.0d0
         RBAS(3,2)=0.0d0
       
         
C
C  WRITE HEADER
C

         DO I=1,3
           DO J=1,3
             HAM(I,J)=0.0d0
           END DO
         END DO
      DO I=1,3
      SUMCHARGE(I)=0.0d0
      END DO
      NPILE=NMSH/NMAX
       FCGRP=1.0D0/NGRP
       write(6,*) 'QDIPOLE: ', NMSH, NMAX, NGRP, FCGRP
       DO 1150 IPILE=0,NPILE
        NOFS=IPILE*NMAX
        MPTS=MIN(NMAX,NMSH-NOFS)
        DO IPTS=1,MPTS
         Q(IPTS,1)=RMSH(1,IPTS+NOFS)
         Q(IPTS,2)=RMSH(2,IPTS+NOFS)
         Q(IPTS,3)=RMSH(3,IPTS+NOFS)
         V(  IPTS)=WMSH(IPTS+NOFS)*FCGRP
        END DO
        DO 1100 IGP=1,NGRP
         DO J=1,3
          DO IPTS=1,MPTS
           P(IPTS,J)=0.0D0
          END DO
          DO K=1,3
           DO IPTS=1,MPTS
            P(IPTS,J)=P(IPTS,J)+RMAT(K,J,IGP)*Q(IPTS,K)
           END DO
          END DO
         END DO

           DO IPTS=1,MPTS
               DO I=1,3
                 SS(I)=0.0D0
               END DO
               DO I=1,3
                 SS(I)=ROTI(I,1)*P(IPTS,1)+
     &             ROTI(I,2)*P(IPTS,2)+ROTI(I,3)*P(IPTS,3)
               END DO
               DO I=1,3
                 P(IPTS,I)=SS(I)
               END DO
           END DO

           DO IWF=1,MAX_OCC
            DO IPTS=1,MPTS
             PSIG(IPTS,IWF)=0.0D0
            END DO
           END DO

         ISHELLA=0
          DO IFNCT=1,NFNCT
           LMAX1=LSYMMAX(IFNCT)+1
           DO I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
            DO J_POS=1,M_NUC
             CALL WFRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &           RVECA,L_NUC,1,INFOWF,NWAVF,ITOTWF)
             IF(L_NUC .NE. M_NUC)THEN
              PRINT *,'DIR: PROBLEM IN WFRAVEL'
              CALL STOPIT
             END IF

             CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
C
C UPDATE PSIG
C

              IF (IUPDAT) THEN
              DO IWF=1,ITOTWF
               ILOC=0
               DO LI=1,LMAX1
                DO MU=1,ISIZE(LI)
                 DO ICON=1,N_CON(LI,IFNCT)
                  ILOC=ILOC+1
                    FACTOR=PSI(ILOC,IWF,1)
                    DO IPV=1,NPV
                     PSIG(JOFS+IPV,IWF)=PSIG(JOFS+IPV,IWF)
     &                                 +FACTOR*GRAD(IPV,1,MU,ICON,LI)
                    END DO
                 END DO
                END DO
               END DO
              END DO
             END IF
            END DO
           END DO
           END DO

C        NOW INTEGRATE

          DO IPTS=1,MPTS
               DO I=1,3
                 SS(I)=0.0D0
               END DO
               DO I=1,3
                 SS(I)=ROT(I,1)*P(IPTS,1)+
     &             ROT(I,2)*P(IPTS,2)+ROT(I,3)*P(IPTS,3)
               END DO
               DO I=1,3
                 P(IPTS,I)=SS(I)
               END DO
           END DO

         DO I=1,3
            DO IPTS=1,MPTS  
             DO IX=1,3
              HAM(I,IX)=HAM(I,IX)+P(IPTS,IX)*PSIG(IPTS,I)*V(IPTS)
             END DO
              SUMCHARGE(I)=SUMCHARGE(I)+PSIG(IPTS,I)**2*V(IPTS)
           END DO
         END DO

 1100   CONTINUE
 1150   CONTINUE

        WRITE(73,*) 'X*ROT(PSI)'
        DO I=1,ITOTWF
          WRITE(73,*) (HAM(I,IX),IX=1,3)
        END DO
        WRITE(73,*) 'CHARGE'
        DO I=1,ITOTWF
          WRITE(73,*) SUMCHARGE(I)
        END DO

         stop
       NFILE=0
        FILE='DIRW'
        DO II=1,ITOTWF 
           NFILE=NFILE+1
           WRITE(FILENAME(NFILE),'(A,I2.2)') FILE, II
           IUNIT=80+NFILE
           OPEN(IUNIT, FILE=FILENAME(NFILE))
           REWIND(IUNIT)
           WRITE(IUNIT,*) 'ORBITAL ROTATED'
           WRITE(IUNIT,'(A4,I1)') FILE,II
           OPEN(77,FILE='XMOL.DAT')
           REWIND(77)
           READ(77,*) NATOM
           READ(77,*)
           WRITE(IUNIT,'(1X,I10,3F20.12)') NATOM,(RBAS(J,1),J=1,3)
           DO K=1,3
            WRITE(IUNIT,'(1X,I10,3F20.12)') NGRID(K),(RBAS(J,K+1),J=1,3)
           ENDDO
           DO K=1,NATOM
            READ(77,*)IZ, X, Y, Z
            CHR=REAL(IZ)
            WRITE(IUNIT,2002)IZ, CHR, X, Y, Z
           END DO
           CLOSE(77)
         END DO
 2002  FORMAT(I6,4F16.10)

C    

       DO 795 IX=1,NGRID(1)
        DO 790 IY=1,NGRID(2)
         NMSH=NGRID(3)
         DO 780 IOFS=0,NMSH-1,MPBLOCK
          MPTS=MIN(MPBLOCK,NMSH-IOFS)
C
C SETUP GRID POINTS AND INITIALIZE PSIG
C
           DO IPTS=1,MPTS
            IZ=IOFS+IPTS
            DO I=1,3
             RGRID(I,IPTS)=RBAS(I,1)+(IX-1)*RBAS(I,2)
     &                    +(IY-1)*RBAS(I,3)+(IZ-1)*RBAS(I,4)
            END DO
           END DO

           DO IPTS=1,MPTS
               DO I=1,3
                 SS(I)=0.0D0
               END DO
               DO I=1,3
                 SS(I)=ROTI(I,1)*RGRID(1,IPTS)+
     &             ROTI(I,2)*RGRID(2,IPTS)+ROTI(I,3)*RGRID(3,IPTS)
               END DO
               DO I=1,3
                 RGRID(I,IPTS)=SS(I)
               END DO
           END DO

           DO IWF=1,MAX_OCC
            DO IPTS=1,MPTS
             PSIG(IPTS,IWF)=0.0D0
            END DO
           END DO
C
C LOOP OVER ALL FUNCTION SETS, THEIR POSITIONS, EQUIVALENT SITES
C
          ISHELLA=0
          DO IFNCT=1,NFNCT
           LMAX1=LSYMMAX(IFNCT)+1
           DO I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
            DO J_POS=1,M_NUC
             CALL WFRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &           RVECA,L_NUC,1,INFOWF,NWAVF,ITOTWF)
             IF(L_NUC .NE. M_NUC)THEN
              PRINT *,'QDIPOLE: PROBLEM IN WFRAVEL'
              CALL STOPIT
             END IF


C
C FOR ALL MESHPOINTS IN BLOCK DO A SMALLER BLOCK
C
             DO JOFS=0,MPTS-1,NSPEED
              NPV=MIN(NSPEED,MPTS-JOFS)
              DO IPV=1,NPV
               PTS(IPV,1)=RGRID(1,JOFS+IPV)-RVECA(1,J_POS)
               PTS(IPV,2)=RGRID(2,JOFS+IPV)-RVECA(2,J_POS)
               PTS(IPV,3)=RGRID(3,JOFS+IPV)-RVECA(3,J_POS)
              END DO
                   
                
c           DO IPTS=1,NPV
c               DO I=1,3
c                 SS(I)=0.0D0
c               END DO
c               DO I=1,3
c                 SS(I)=ROTI(I,1)*PTS(IPTS,1)+
c     &             ROTI(I,2)*PTS(IPTS,2)+ROTI(I,3)*PTS(IPTS,3)
c               END DO
c               DO I=1,3
c                 PTS(IPTS,I)=SS(I)
c               END DO
c           END DO

C
C GET VALUE OF BASIS FUNCTIONS
C
              CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
C
C UPDATE PSIG
C
 
              IF (IUPDAT) THEN
              DO IWF=1,ITOTWF
               ILOC=0
               DO LI=1,LMAX1
                DO MU=1,ISIZE(LI)
                 DO ICON=1,N_CON(LI,IFNCT)
                  ILOC=ILOC+1
                    FACTOR=PSI(ILOC,IWF,1)
                    DO IPV=1,NPV
                     PSIG(JOFS+IPV,IWF)=PSIG(JOFS+IPV,IWF)
     &                                 +FACTOR*GRAD(IPV,1,MU,ICON,LI)
                    END DO
                 END DO
                END DO
               END DO
              END DO
             END IF
            END DO
           END DO
          END DO
         END DO
C
C GET ORBITAL DENSITY FROM WAVEFUNCTION
C

      NFILE=0
       DO IWF=1,ITOTWF
        DO IPTS=1,MPTS
          EVEC(IPTS)=0.0D0
        END DO
        NFILE=NFILE+1
         DO IPTS=1,MPTS
c          DO IWF=1,ITOTWF
c         EVEC(IPTS)=EVEC(IPTS)+HAM(IX,IWF)*PSIG(IPTS,IWF)
c         EVEC(IPTS)=EVEC(IPTS)+H(IWF,IX)*PSIG(IPTS,IWF)
         EVEC(IPTS)=PSIG(IPTS,IWF)
c         END DO
         END DO
                  
        IUNIT=80+NFILE
        WRITE(IUNIT,'(3(1X,E20.12))')(EVEC(IPTS)*FACT,
     &        IPTS=1,MPTS)
       END DO

  780  CONTINUE
  790  CONTINUE
  795  CONTINUE

      NFILE=0
       DO IWF=1,ITOTWF
        NFILE=NFILE+1
        IUNIT=80+NFILE
        CLOSE(IUNIT)
       END DO

  900  CONTINUE
       CALL GTTIME(TIME2)
       CALL TIMOUT('WF DIR :                   ',TIME2-TIME1)
       RETURN
       END
