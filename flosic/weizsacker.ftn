C COMPUTE THE RATIO OF TAU_W AND TAU_TOTAL
C BASED ON DENSOLD ROUTINE
C MODE 1: Make RSIC file
C MODE 2: Make alpha vs. radius file (and sort in term of radius)

C DENSOLD BASED ON OLD VERSION OF APOTNL BY M. PEDERSON AND D. POREZAG
C ATTENTION: FIRST TWO ARRAYS OF COMMON BLOCK TMP1 MUST BE IDENTICAL IN 
C DENSOLD AND APOTNL SINCE THEY ARE USED TO PASS DENSITY AND COULOMB POT
C
       SUBROUTINE WEIZSACKER(MODE)
C      !use pot_dens,only : rhog
       use mesh1,only : rmsh,nmsh
       use common2,only : RIDT, N_CON, LSYMMAX, N_POS, NFNCT, IGGA,
     &     ISPN, NSPN
       use common5,only : PSI, NWF, NWFS
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:48 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       LOGICAL :: ALPHACHECK
       INTEGER,INTENT(IN) :: MODE
       INTEGER :: I_POS, ICON, IERR, IFNCT, IGR, ILOC, IPTS, ISHDUM,
     &   ISHELLA, ISIZE, IWF, J_POS, JBEG, JGR, JLOC, JPTS, JWF, KGR,
     &   KPTS, L_NUC, LI, LMAX1, LPTS, LPV, M_NUC, MPTS, MU, NDERV,
     &   NGRAD, NMAX, NPV, NSHELLA, OFFSET, COUNTER1, COUNTER2,IUNIT
       REAL*8 :: SYMBOL, TIMEGORB, APT1, FACTOR, TIME3, TIME4, PI, 
     &           RAD ,TMP
       CHARACTER*20 :: FNAME(2)
       SAVE
       PARAMETER (NMAX=MPBLOCK)
       PARAMETER (PI=3.14159265359d0)
C
C RETURN:
C RHOG(IPTS,1, 1)= rho_up   
C RHOG(IPTS,2, 1)= d rho_up/dx
C RHOG(IPTS,3, 1)= d rho_up/dy
C RHOG(IPTS,4, 1)= d rho_up/dz
C RHOG(IPTS,1, 2)= rho_dn   
C RHOG(IPTS,2, 2)= d rho_dn/dx
C RHOG(IPTS,3, 2)= d rho_dn/dy
C RHOG(IPTS,4, 2)= d rho_dn/dz
C
C       LOGICAL ICOUNT
C       COMMON/TMP1/ACOULOMB(MAX_PTS),ARHOG(MAX_PTS,KRHOG,MXSPN)
C     &  ,PHIG(MAX_PTS,2)
C       COMMON/TMP1/PHIG(MAX_PTS,2)
C
C       COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
C     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
C     &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       REAL*8,ALLOCATABLE :: PSIG(:,:,:),PTS(:,:)
     &                      ,GRAD(:,:,:,:,:),RVECA(:,:)
     &                      ,RHOG(:,:,:),TAU(:,:),TAUW(:,:)
     &                      ,TAUUNIF(:,:),ALPHA(:,:)
     &                      ,RADIUS(:),ALPHA1D(:)
       LOGICAL,ALLOCATABLE :: ICOUNT(:,:)

       LOGICAL LGGA,IUPDAT
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/
C
       TIMEGORB=0.0D0
       CALL GTTIME(APT1)
       ALLOCATE(PSIG(NMAX,10,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING PSIG'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING NSPEED'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING GRAD'
       ALLOCATE(RVECA(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING RVECA'
       ALLOCATE(ICOUNT(MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING ICOUNT'
       !LGGA= .FALSE.
       NGRAD=4

       ALLOCATE(RHOG(NMSH,NGRAD,NSPN))
       ALLOCATE( TAU(NMSH,NSPN))
       ALLOCATE(TAUW(NMSH,NSPN))
       ALLOCATE(TAUUNIF(NMSH,NSPN))
       ALLOCATE(ALPHA(NMSH,NSPN))

       RHOG=0.0d0
       TAU =0.0d0
       TAUW=0.0d0
       TAUUNIF=0.0d0
       ALPHA=0.0d0
C
C LOOP OVER ALL POINTS
C
       LPTS=0
 10    CONTINUE
        IF(LPTS+NMAX.LT.NMSH)THEN
         MPTS=NMAX
        ELSE
         MPTS=NMSH-LPTS
        END IF
C
C INITIALIZE PSIG AND RHOB
C
        DO IWF=1,NWF
         DO IGR=1,NGRAD
          DO IPTS=1,MPTS
           PSIG(IPTS,IGR,IWF)=0.0D0
          END DO
         END DO  
        END DO  
        DO ISPN=1,NSPN
         DO IGR=1,NGRAD
          DO IPTS=1,MPTS
           RHOG(LPTS+IPTS,IGR,ISPN)=0.0D0
          END DO
         END DO  
        END DO  
        ISHELLA=0
C
C FOR ALL CENTER TYPES
C
        DO 86 IFNCT=1,NFNCT
         LMAX1=LSYMMAX(IFNCT)+1
C
C FOR ALL POSITIONS OF THIS CENTER
C
         DO 84 I_POS=1,N_POS(IFNCT)
          ISHELLA=ISHELLA+1
C
C GET SYMMETRY INFO
C
          CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
          IF(NWF.GT.MAX_OCC)THEN
           write(6,*)'DENSOLD: MAX_OCC MUST BE AT LEAST:',NWF
           CALL STOPIT
          END IF
C
C FOR ALL EQUIVALENT POSITIONS OF THIS ATOM
C
          DO 82 J_POS=1,M_NUC
C
C UNSYMMETRIZE 
C
           CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                  RVECA,L_NUC,1)
           IF(L_NUC.NE.M_NUC)THEN
            write(6,*)'DENSOLD: PROBLEM IN UNRAVEL'
            CALL STOPIT
           END IF
C
C FOR ALL MESHPOINTS IN BLOCK DO A SMALLER BLOCK
C
           KPTS=0
           DO 80 JPTS=1,MPTS,NSPEED
            NPV=MIN(NSPEED,MPTS-JPTS+1)
            DO LPV=1,NPV
             KPTS=KPTS+1
             PTS(LPV,1)=RMSH(1,LPTS+KPTS)-RVECA(1,J_POS)
             PTS(LPV,2)=RMSH(2,LPTS+KPTS)-RVECA(2,J_POS)
             PTS(LPV,3)=RMSH(3,LPTS+KPTS)-RVECA(3,J_POS)
            END DO
C
C GET ORBITS AND DERIVATIVES
C
            NDERV=1
            !IF (LGGA) NDERV=2
            CALL GTTIME(TIME3)
            CALL GORBDRV(NDERV,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
            CALL GTTIME(TIME4)
            TIMEGORB=TIMEGORB+TIME4-TIME3
C
C UPDATING ARRAY PSIG
C
            IF (IUPDAT) THEN
             IPTS=JPTS-1
             ILOC=0
             DO 78 LI=1,LMAX1
              DO MU=1,ISIZE(LI)
               DO ICON=1,N_CON(LI,IFNCT)
                ILOC=ILOC+1
                IF (ICOUNT(ICON,LI)) THEN
                 DO IWF=1,NWF
                  FACTOR=PSI(ILOC,IWF,1)
                  if(abs(FACTOR).GT.1.0d-10) then 
                   DO IGR=1,NGRAD
                   DO LPV=1,NPV
                    PSIG(IPTS+LPV,IGR,IWF)=PSIG(IPTS+LPV,IGR,IWF)
     &              +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
                   END DO
                   END DO 
                  end if 
                 END DO  
                END IF
               END DO  
              END DO  
   78        CONTINUE
            END IF
   80      CONTINUE
   82     CONTINUE
   84    CONTINUE
   86   CONTINUE

C
C UPDATE DERIVATIVES IF GGA CALCULATION
C         
        DO 96 ISPN=1,NSPN
         JBEG= (ISPN-1)*NWFS(1)
         DO 94 JWF=1,NWFS(ISPN)
          JLOC=JWF+JBEG
C
C DENSITY
C
          DO IPTS=1,MPTS
           RHOG(LPTS+IPTS,1,ISPN)=RHOG(LPTS+IPTS,1,ISPN)
     &     +PSIG(IPTS,1,JLOC)**2
          END DO
C
C GRADIENT 
C
          DO IGR=2,4
           DO IPTS=1,MPTS
            RHOG(LPTS+IPTS,IGR,ISPN)=RHOG(LPTS+IPTS,IGR,ISPN)
     &      +2*PSIG(IPTS,1,JLOC)*PSIG(IPTS,IGR,JLOC)
           END DO
          END DO
C
C TAU = 1/2 * |Grad Psi|**2 
C
          DO IPTS=1,MPTS
           TAU(LPTS+IPTS,ISPN)=TAU(LPTS+IPTS,ISPN)+ 0.5d0*
     &      (PSIG(IPTS,2,JLOC)*PSIG(IPTS,2,JLOC)
     &      +PSIG(IPTS,3,JLOC)*PSIG(IPTS,3,JLOC)
     &      +PSIG(IPTS,4,JLOC)*PSIG(IPTS,4,JLOC))
          END DO

   94    CONTINUE
   96   CONTINUE

        LPTS=LPTS+MPTS
        IF (LPTS .LT. NMSH) GOTO 10

       CONTINUE

       DO ISPN=1,NSPN
        DO IPTS=1,NMSH
C
C TAU W = |Grad n|**2 / 8n (Note: TAU W < TAU)
C
           TAUW(IPTS,ISPN)=
     &      (RHOG(IPTS,2,ISPN)*RHOG(IPTS,2,ISPN)
     &      +RHOG(IPTS,3,ISPN)*RHOG(IPTS,3,ISPN)
     &      +RHOG(IPTS,4,ISPN)*RHOG(IPTS,4,ISPN))
     &      /(8.0d0*RHOG(IPTS,1,ISPN))
C
C TAU Uniform = (3/10)(3 pi**2)**(2/3)  n**(5/3)
C
           TAUUNIF(IPTS,ISPN)=
     &      (3.0d0/10.0d0)*((3.0d0*pi*pi)**(2.0d0/3.0d0))
     &      *(RHOG(IPTS,1,ISPN)**(5.0d0/3.0d0))
C
C ALPHA (Note: ALPHA is > 0 always.)
C 
           ALPHA(IPTS,ISPN)=
     &      (TAU(IPTS,ISPN) - TAUW(IPTS,ISPN))/
     &       TAUUNIF(IPTS,ISPN)
        END DO
       END DO


C
C SAVE THE INFORMATION IN FILE
C
!==== MODE 0: Save the ratio of tau on grid ===
       IF(MODE .EQ. 0) THEN
        OPEN(214,FILE='tauw_tau1.dat',STATUS='UNKNOWN')
        CLOSE(214,STATUS='DELETE')

        OPEN(215,FILE='tauw_tau2.dat',STATUS='UNKNOWN')
        CLOSE(215,STATUS='DELETE')

        OPEN(214,FILE='tauw_tau1.dat',FORM='FORMATTED',
     &       STATUS='NEW')
         REWIND(214)
         WRITE(214,*) (TAUW(IPTS,1)/TAU(IPTS,1),IPTS=1,NMSH)
        CLOSE(214)

        IF(NSPN.EQ.2) THEN
         OPEN(215,FILE='tauw_tau2.dat',FORM='FORMATTED',
     &        STATUS='NEW')
         WRITE(215,*) (TAUW(IPTS,2)/TAU(IPTS,2),IPTS=1,NMSH)
         CLOSE(215)
        END IF

        !Saving densities 
!       OPEN(214,FILE='rho1.dat',STATUS='UNKNOWN')
!       CLOSE(214,STATUS='DELETE')

!       OPEN(215,FILE='rho2.dat',STATUS='UNKNOWN')
!       CLOSE(215,STATUS='DELETE')

!       OPEN(214,FILE='rho1.dat',FORM='FORMATTED',
!    &       STATUS='NEW')
!        REWIND(214)
!        WRITE(214,*) (RHOG(IPTS,1,1),IPTS=1,NMSH)
!       CLOSE(214)

!       IF(NSPN.EQ.2) THEN
!        OPEN(215,FILE='rho2.dat',FORM='FORMATTED',
!    &        STATUS='NEW')
!        WRITE(215,*) (RHOG(IPTS,1,2),IPTS=1,NMSH)
!        CLOSE(215)
!       END IF

!==== MODE 1: Save alpha on the grid ====
       ELSE IF(MODE .EQ. 1) THEN
        OPEN(214,FILE='WEIZSACKER.dat',FORM='FORMATTED',
     &       STATUS='UNKNOWN')
        CLOSE(214,STATUS='DELETE')
        OPEN(214,FILE='WEIZSACKER.dat',FORM='FORMATTED',STATUS='NEW')
        IF(NSPN.EQ.1) THEN
         DO IPTS=1,NMSH
          RAD=sqrt(RMSH(1,IPTS)**2 + RMSH(2,IPTS)**2 + RMSH(3,IPTS)**2)
!         WRITE(214,*) RAD, TAU(IPTS,1)-TAUW(IPTS,1), ALPHA(IPTS,1)
!         WRITE(214,*) RAD, TAUW(IPTS,1)/TAU(IPTS,1),ALPHA(IPTS,1)
!         WRITE(214,*) TAU(IPTS,1),TAUW(IPTS,1),TAUUNIF(IPTS,1)
!         WRITE(214,*) RAD, TAUW(IPTS,1)/TAU(IPTS,1)
          WRITE(214,'(3F15.5,F20.6)') RMSH(1:3,IPTS),ALPHA(IPTS,1)
         END DO
        ELSE
         DO IPTS=1,NMSH
          WRITE(214,'(5F20.8,2F30.8)') RMSH(1:3,IPTS),
     &              TAUW(IPTS,1)/TAU(IPTS,1),
     &              TAUW(IPTS,2)/TAU(IPTS,2),ALPHA(IPTS,1),ALPHA(IPTS,2)
         END DO
        END IF
        CLOSE(214)
C
C 
C
!==== MODE 2: Boundary Finder based on alpha====
       ELSE IF(MODE .EQ. 2) THEN
        OPEN(215,FILE='ALPHAVSR',FORM='FORMATTED',STATUS='UNKNOWN')
        CLOSE(215,STATUS='DELETE')
        OPEN(215,FILE='ALPHAVSR',FORM='FORMATTED',STATUS='NEW')

!       OPEN(216,FILE='BOUNDARIES',FORM='FORMATTED',STATUS='UNKNOWN'
!    &  ,ACCESS='APPEND')
        OPEN(216,FILE='BOUNDARIES',FORM='FORMATTED',STATUS='UNKNOWN')
        CLOSE(216,STATUS='DELETE')
        OPEN(216,FILE='BOUNDARIES',FORM='FORMATTED',STATUS='NEW')

!       ISHELLA=0
!       ISHELLA is still in the memory
        print *,"ISHELLA is currently: ", ISHELLA
        NSHELLA=ISHELLA
! C
!           RIDT(1,ISHELLA)
! Note: this radius is centered at the coordinate origin. This needs to
! be centered at nucleii. 

        ALLOCATE(RADIUS(0:NMSH*NSPN*NSHELLA))
        RADIUS(0)=0.0D0

        DO ISHELLA=1,NSHELLA
         OFFSET=NMSH*NSPN*(ISHELLA-1)
         DO IPTS=1,NMSH
          RADIUS(IPTS+OFFSET)
     &      =sqrt((RMSH(1,IPTS)-RIDT(1,ISHELLA))**2 
     &           +(RMSH(2,IPTS)-RIDT(2,ISHELLA))**2 
     &           +(RMSH(3,IPTS)-RIDT(3,ISHELLA))**2)
         END DO
         IF(NSPN .EQ. 2) RADIUS(NMSH+1+OFFSET:2*NMSH+OFFSET)
     &                    =RADIUS(1+OFFSET:NMSH+OFFSET)
        END DO

        ALLOCATE(ALPHA1D(NSPN*NMSH*NSHELLA))
        ALPHA1D(1:NMSH)=ALPHA(:,1)
        IF(NSPN .EQ. 2) ALPHA1D(NMSH+1:2*NMSH)=ALPHA(:,2)

        COUNTER1=0
        COUNTER2=0
        DO IPTS=1,NMSH*NSPN
         DO JPTS=IPTS+1,NMSH*NSPN
          IF(RADIUS(JPTS) .LT. RADIUS(IPTS)) THEN
            TMP=RADIUS(JPTS)
            RADIUS(JPTS)=RADIUS(IPTS)
            RADIUS(IPTS)=TMP
            TMP=ALPHA1D(JPTS)
            ALPHA1D(JPTS)=ALPHA1D(IPTS)
            ALPHA1D(IPTS)=TMP
          END IF
         END DO
         IF(ABS(RADIUS(IPTS)-RADIUS(IPTS-1)) .GT. 1.0D-8) THEN
!          WRITE(215,*) RADIUS(IPTS),ALPHA1D(IPTS)
          COUNTER1=COUNTER1+1
         END IF
        END DO
        WRITE(215,*) COUNTER1
        DO IPTS=1,NMSH*NSPN
         IF(ABS(RADIUS(IPTS)-RADIUS(IPTS-1)) .GT. 1.0D-8) THEN 
          WRITE(215,*) RADIUS(IPTS),ALPHA1D(IPTS)
          IF(ALPHACHECK(ALPHA1D(IPTS))) THEN
           COUNTER2=COUNTER2+1
           RADIUS(COUNTER2)=RADIUS(IPTS)
           ALPHA1D(COUNTER2)=ALPHA1D(IPTS)
          END IF
         END IF
        END DO

        WRITE(216,*) COUNTER2
        DO IPTS=1,COUNTER2
         WRITE(216,*) RADIUS(IPTS),ALPHA1D(IPTS)
        END DO

        DEALLOCATE(ALPHA1D)
        DEALLOCATE(RADIUS)
        CLOSE(216)
        CLOSE(215)
!==== MODE 3: Save alpha for gaussian cube file  ====
       ELSE IF(MODE .EQ. 3) THEN
        FNAME(1)='ALPHAGRID'
        IF(NSPN.EQ.1) THEN
         FNAME(1)=FNAME(1)(1:3)//'UPO'
        ELSE
         FNAME(1)=FNAME(1)(1:3)//'SUP' 
         FNAME(2)=FNAME(1)(1:3)//'SDN'
        ENDIF
        DO ISPN=1,NSPN
         IUNIT=73+ISPN
         OPEN(IUNIT,FILE=FNAME(ISPN),FORM='FORMATTED',
     &              STATUS='OLD',ACCESS='APPEND')
         WRITE(IUNIT,9010) (ALPHA(IPTS,ISPN),IPTS=1,NMSH)
         CLOSE(IUNIT)
        END DO
 9010   FORMAT(3(1X,E20.12))
!==== MODE 4: Save beta to a file ====
       ELSE IF(MODE .EQ. 4) THEN
        OPEN(217,FILE='beta1.dat',STATUS='UNKNOWN')
        CLOSE(217,STATUS='DELETE')

        OPEN(217,FILE='beta2.dat',STATUS='UNKNOWN')
        CLOSE(217,STATUS='DELETE')

        OPEN(217,FILE='beta1.dat',FORM='FORMATTED',
     &       STATUS='NEW')
         REWIND(217)
         WRITE(217,*) ((TAU(IPTS,1)-TAUW(IPTS,1))/
     &                 (TAU(IPTS,1)+TAUUNIF(IPTS,1)),IPTS=1,NMSH)
        CLOSE(217)

        IF(NSPN.EQ.2) THEN
         OPEN(217,FILE='beta2.dat',FORM='FORMATTED',
     &        STATUS='NEW')
         WRITE(217,*) ((TAU(IPTS,2)-TAUW(IPTS,2))/
     &                 (TAU(IPTS,2)+TAUUNIF(IPTS,2)),IPTS=1,NMSH)
         CLOSE(217)
        END IF
!==== MODE 5: Save total density to files ====
       ELSE IF(MODE .EQ. 5) THEN
        OPEN(214,FILE='rho1.dat',STATUS='UNKNOWN')
        CLOSE(214,STATUS='DELETE')

        OPEN(215,FILE='rho2.dat',STATUS='UNKNOWN')
        CLOSE(215,STATUS='DELETE')

        OPEN(214,FILE='rho1.dat',FORM='FORMATTED',
     &       STATUS='NEW')
         REWIND(214)
         WRITE(214,*) (RHOG(IPTS,1,1),IPTS=1,NMSH)
        CLOSE(214)

        IF(NSPN.EQ.2) THEN
         OPEN(215,FILE='rho2.dat',FORM='FORMATTED',
     &        STATUS='NEW')
         WRITE(215,*) (RHOG(IPTS,1,2),IPTS=1,NMSH)
         CLOSE(215)
        END IF
       END IF

       DEALLOCATE(ALPHA)
       DEALLOCATE(TAUUNIF)
       DEALLOCATE(TAUW)
       DEALLOCATE(TAU)
       DEALLOCATE(RHOG)

       DEALLOCATE(PSIG,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING PSIG'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING NSPEED'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING GRAD'
       DEALLOCATE(RVECA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING RVECA'
       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING ICOUNT'
       RETURN
       END

!##################################################
! Function to check if you need a fine mesh 
! at a given point, i, with alpha, a(i).
!##################################################
      logical function alphacheck(x)
      real*8, intent(in) :: x
      !integer :: i

      !alphacheck=.false.
      !do i=1,10
      if(x > 0.8d0 .and. x < 1.2d0) then
        alphacheck=.true.
      else
        alphacheck=.false.
      endif
      !end do
      return
      end
