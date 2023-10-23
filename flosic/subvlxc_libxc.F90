! UTEP Electronic Structure Lab (2020)
!> @file subvlxc_libxc.F91
!> @author Yoh Yamamoto
!> @details subvlxc_libxc.F91 is a substantially edited version of 
!>  subvlxc.ftn. It is capable of calling Libxc Library. To use 
!>  libxc, look up a functional string on 
!>  https://www.tddft.org/programs/libxc/functionals/ 
!>  and use it in the CLUSTER/SYMBOL file
!
!> @note  5/12/2015 SUBVLXC modified to accomodate LIBXC library
!>       11/15/2015 Working on meta-GGA
!>       12/2015 MGGA SCAN functioanal is added
!
! ***************original description (LDA/PBE only subvlxc)****************
!
! SUBVLXC DIRK POREZAG AUGUST 1999
! CALCULATES THE LOCAL AND EXCHANGE-CORRELATION POTENTIALS
! FOR A SET OF POINTS
!
! MODE:  1: VXC ONLY
!        2: VXC AND VLOCAL
!
subroutine SUBVLXC(MODE,LPTS,MPTS,RHOV,VXCS,VLOS,EXCVEC,LDFTF,MXXD,KXXS,KXXO)

use debug1
use mesh1,only    : WMSH,RMSH
use common2,only  : RCNT, IFUCNT, NCNT, IGGA, IDFTYP, ISPN, NSPN
!use common8,only : REP, N_REP, NDMREP, NS_TOT
use global_inputs, only : libxc1   !<== Libxc switch moved to NRLMOL_INPUT.dat
use common5,only : CONVERGENCE
use XTMP2A,ONLY : TAUCHOP,ISMGGA,MIXINS,TAUW_TAU,BETA
use scaledpzsic,only : SICEXCS,SDSICON,GSICON,scaledsic
INCLUDE 'PARAMA2'
integer :: I,ICNT,IERR,IFU,IOFS,IPV,JPTS,MXDFTYP,NGRAD,NPV
real*8 :: THIRD,THIRD2,TRPI2,FRDPI,DKF2RS,EPS,DELTA,FAC, &
          EX,DEX,EXL,EXN,EC,DEC,ECL,ECN,                 &
          D,DKF,S,U,V,DGRAD,DLAP,VX,DUP,DDN,             &
          ZET,DKS,RS,G,T,UU,VV,WW,ECRS,ECZET,ALFC,       &
          REC1,REC2,REC22,REC3,DREC,                     &
          timing1,timing2!,                               &
!          sig1,sig2,sig3
SAVE
PARAMETER (MXDFTYP=15)        ! This value must match with NDFTYP in setdftyp.ftn
integer, intent(in)  :: MODE
integer, intent(in)  :: LPTS  ! Subvlxc is called with a LPTS loop from getvlxc
integer, intent(in)  :: MPTS  ! Number of points to be processed
real(8), intent(in)  :: RHOV(10*MXSPN*MPBLOCKXC)
!The exchange-correlation potentials (in phi_i phi_j terms)
real(8), intent(out) :: VXCS(MXSPN*MPBLOCKXC)        ! Only using MPTS*2 out of MPBLOCK*2
real(8), intent(out) :: VLOS(MPBLOCKXC)
real(8), intent(out) :: EXCVEC(4)
logical, intent(in)  :: LDFTF  ! T:DFT, F:SIC
integer, intent(in)  :: MXXD,KXXS,KXXO   ! needed for unravel2 in gettau

real(8) :: ISPFAC
real(8), allocatable :: RHOC(:,:),  &
                        XTMP(:,:),  &
                        DTMP(:,:),  &
                        RTMP(:),    &
                        VLOC(:)
!YY: Additional variables - dE/drho terms for exchange and correlation terms
real(8), allocatable :: vxloc(:), vclocup(:), vclocdn(:)
!YY: Inputs required for libxc calculation. Sending them to call_libxc_unp/pol
!subroutines
real(8), allocatable :: rhoin(:),      &  !Density
                        sigmain(:),    &  !Contracted gradient of density
                        sigma1(:), sigma2(:), sigma3(:), &
                        drho(:,:,:),   &  !An array to hold derivative terms
                        dense(:,:),    &
                        Eex(:),        &
                        Eco(:)
! Laplacians for meta-GGA calculation
real(8), allocatable :: laplace(:,:), laplacein(:)
! Kinetic density for meta-GGA calculation
!real(8), intent(in)  :: tauv(MPBLOCK,MXSPN)!, rhodebug(MPBLOCK,MXSPN)
real(8), allocatable :: tauin(:), tauup(:), taudn(:)
!To compute Hamiltonian matrix mixing elements, we need to gather 2*vsigma*grad
!rho and send them back to getvlxc. mixins is the array to hold that quantity.
!real(8), intent(out) :: mixins(4,MPBLOCKXC*MXSPN)
!Temporary array to hold the mixing elements from exchange term. One spin at a time.
real(8), allocatable :: mixtmp(:,:)
!Temporary arrays to hold the mixing elements: exchange and correlation terms
real(8), allocatable :: mixex(:,:), mixco(:,:)
!Temporary arrays to hold the mixing elements for correlation spin up and down terms
real(8), allocatable :: mixcorup(:,:), mixcordn(:,:)
!SCAN variables
!real(8) :: VXDD,AMUXD
!real(8) :: VCDD1,VCDD2,AMUCD1,AMUCD2
real(8) :: VCUP,VCDN,DVCUP,DVCDN
!real(8) :: vsigma1,vsigma2,vsigma3
!########################################################################
!#  Libxc switch: libxcv in NRLMOL_INPUT.DAT and libxc1 in global_input #
!########################################################################

! COMMON/PW91GAS/ IS NEEDED FOR PW91 GGA

COMMON/PW91GAS/G,EC,ECRS,ECZET
LOGICAL     ISGGA
integer ::  IPTS
real(8) ::  RHT(10,NSPN),    &
            VCOR(2),         &
            DEN(3),          &
            DG2(3),          &
            DGG(3),          &
            DLP(3)
DATA THIRD/0.3333333333333333D0/
DATA THIRD2/0.6666666666666667D0/
DATA TRPI2/29.6088132032680759D0/      ! 3*pi*pi
DATA FRDPI/1.2732395447351627D0/
DATA DKF2RS/1.9191582926775128D0/
DATA EPS/1.0D-20/

! CHECK IGGA, IDFTYP, AND MPTS
!CALL TRACER('STARTING SUBVLXC')
do I=1,2
   if(libxc1) then
      if(IGGA(I) < 0 .or. IGGA(I) > 5) then
         write(*,*)'SUBVLXC: Invalid value for IGGA (with libxc): ',IGGA(I)
         call STOPIT
      end if
   else
      if(IGGA(I) < 0 .or. IGGA(I) > 3) then
      !if (IGGA(I)*(IGGA(I)-1) /= 0) THEN
         write(6,*)'SUBVLXC: INVALID VALUE FOR IGGA: ',IGGA(I)
         CALL STOPIT
      end if
   end if
   !Skip checking IDFTYP when libxc is turned on.
   if(.not. libxc1) then
      if ((IDFTYP(I) < 0) .OR. (IDFTYP(I) > MXDFTYP)) THEN
         write(6,*)'SUBVLXC: INVALID VALUE FOR IDFTYP: ',IDFTYP(I)
         CALL STOPIT
      end if
   end if
end do
if (MPTS > MPBLOCKXC) THEN
   write(6,*)'SUBVLXC: MPTS MUST BE <= MPBLOCKXC'
   CALL STOPIT
end if

! ALLOCATE LOCAL ARRAYS

allocate(RHOC(10,MPBLOCKXC),STAT=IERR); if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR ALLOCATING RHOC'
allocate(XTMP(3,MPBLOCKXC),STAT=IERR);  if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR ALLOCATING XTMP'
allocate(DTMP(3,MPBLOCKXC),STAT=IERR);  if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR ALLOCATING DTMP'
allocate(RTMP(MPBLOCKXC),STAT=IERR);    if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR ALLOCATING RTMP'
allocate(VLOC(NSPEED),STAT=IERR);     if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR ALLOCATING VLOC'
allocate(TAUCHOP(MPBLOCKXC,MXSPN),STAT=IERR); if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR ALLOCATING TAUCHOP'

if(ISMGGA) then
!  call gttime(timing1)
   CALL GETTAU_PAR(LPTS,TAUCHOP,LDFTF,MXXD,KXXS,KXXO)
!  call GTTIME(timing2)
!  call timout('GETTAU EXECUTION:                   ',timing2-timing1)
else
  tauchop=0.0d0
end if
!CALL TRACER('ALLOCATED ARRAYS')
! LOCAL POTENTIAL
! DELTA IS USED TO PREVENT DIVISION BY ZERO FOR R=0

if (MODE > 1) THEN
   !do IPTS=1,MPTS
   !   VLOS(IPTS)= 0.0D0
   !end do
   VLOS(:)=0.0d0
   DELTA=1.0D-100
   do ICNT=1,NCNT
      IFU=IFUCNT(ICNT)
      do IPTS=0,MPTS-1,NSPEED
         NPV=MIN(NSPEED,MPTS-IPTS)
         IOFS=LPTS+IPTS
         do IPV=1,NPV
           RTMP(IPV)=(RMSH(1,IOFS+IPV)-RCNT(1,ICNT))**2    &
                    +(RMSH(2,IOFS+IPV)-RCNT(2,ICNT))**2    &
                    +(RMSH(3,IOFS+IPV)-RCNT(3,ICNT))**2
                     RTMP(IPV)=MAX(RTMP(IPV),DELTA)
           RTMP(IPV)=SQRT(RTMP(IPV))
         end do
         call VLOCAL(1,NPV,IFU,RTMP,VLOC)
!VLOS is the sum of
         do IPV=1,NPV
           VLOS(IPTS+IPV)=VLOS(IPTS+IPV)+VLOC(IPV)
         end do
      end do
   end do
end if

if(libxc1 .or. ismgga) mixins(:,:)=0.0d0 !uncommented
if((.not.LDFTF).and.scaledsic) then
    SICEXCS(:)=0.0d0
end if

! INITIALIZATION OF DATA NEEDED FOR EXCHANGE-CORRELATION

EXL= 0.0D0
EXN= 0.0D0
ECL= 0.0D0
ECN= 0.0D0
ISPFAC= 2.0d0/NSPN
! YY. ISGGA=.true. for GGA and MGGA
!ISGGA= ((IGGA(1) == 1) .OR. (IGGA(2) == 1))
ISGGA= ((IGGA(1) >= 1) .or. (IGGA(2) >= 1))
NGRAD=1
IF (ISGGA) NGRAD=10

! GET CORE DENSITY
CALL GTRHOCR(ISGGA,MPTS,RMSH(1,LPTS+1),RHOC,XTMP,DTMP,RTMP)
!CALL TRACER('GOT CORE DENSITY')
!###############################################
!# Case I: LIBXC off (Use builtin functionals) #
!###############################################
if(.NOT. libxc1 .OR. MODE==1) then

! LOOP OVER POINTS
   do 200 IPTS=1,MPTS
      VXCS(IPTS)= 0.0D0
      if (NSPN == 2) VXCS(IPTS+MPTS)= 0.0D0

! INITIALIZE DG2, DGG, DLP

      !do I=1,3
      !   DG2(I)= 0.0D0
      !   DGG(I)= 0.0D0
      !   DLP(I)= 0.0D0
      !end do
      DG2(:)=0.0d0
      DGG(:)=0.0d0
      DLP(:)=0.0d0

! MOVE DATA TO RHT

      do ISPN=1,NSPN
         FAC= 1.0D0
         if (ISPN == 2) FAC= 0.5D0
         IOFS=NGRAD*((ISPN-1)+(IPTS-1)*NSPN)
         do I=1,NGRAD
            RHT(I,ISPN)= RHOV(IOFS+I)+FAC*RHOC(I,IPTS)
         end do
         !rhodebug(IPTS,ISPN)= rhodebug(IPTS,ISPN)+FAC*RHOC(1,IPTS)
         if (RHT(1,ISPN) < 0.0D0)    RHT(1,ISPN)= 0.0D0
         if (RHT(1,ISPN) > RHT(1,1)) RHT(1,ISPN)= RHT(1,1)
      end do
!
! DENSITIES D, DUP, DDN
!
      DEN(1)= RHT(1,1)                    !Total density
      DEN(3)= DEN(1)*0.5D0
      if (NSPN == 2) DEN(3)= RHT(1,NSPN)  !Spin down density
      DEN(2)= DEN(1)-DEN(3)               !Spin up density = Total density-spin down

! ABS(GRAD(D))**2
      if (ISGGA) THEN
         DG2(1)= RHT(2,1)**2+RHT(3,1)**2+RHT(4,1)**2
         DG2(2)= DG2(1)*0.250D0
         DG2(3)= DG2(2)

! GRAD(D)xGRAD(ABS(GRAD(D)))*ABS(GRAD(D))

         DGG(1)=  RHT(5,1)*RHT(2,1)**2          &
                 +RHT(6,1)*RHT(3,1)**2          &
                 +RHT(7,1)*RHT(4,1)**2          &
                 +2*(RHT(8,1)*RHT(2,1)*RHT(3,1) &
                 +RHT(9,1)*RHT(2,1)*RHT(4,1)    &
                 +RHT(10,1)*RHT(3,1)*RHT(4,1))
         DGG(2)= DGG(1)*0.125D0
         DGG(3)= DGG(2)

! LAPLACE(D)

         DLP(1)= RHT(5,1)+RHT(6,1)+RHT(7,1)
         DLP(2)= DLP(1)*0.500D0
         DLP(3)= DLP(2)

! EQUIVALENT SPIN-POLARIZED TERMS

         if (NSPN == 2) THEN
            do I=1,10
              RHT(I,1)= RHT(I,1)-RHT(I,NSPN)
            end do
            do ISPN=1,NSPN
               DG2(ISPN+1)= RHT(2,ISPN)**2+RHT(3,ISPN)**2+RHT(4,ISPN)**2
               DGG(ISPN+1)= RHT(5,ISPN)*RHT(2,ISPN)**2   &
                           +RHT(6,ISPN)*RHT(3,ISPN)**2   &
                           +RHT(7,ISPN)*RHT(4,ISPN)**2   &
                        +2*(RHT(8,ISPN)*RHT(2,ISPN)*RHT(3,ISPN)   &
                           +RHT(9,ISPN)*RHT(2,ISPN)*RHT(4,ISPN)   &
                           +RHT(10,ISPN)*RHT(3,ISPN)*RHT(4,ISPN))
               DLP(ISPN+1)= RHT(5,ISPN)+RHT(6,ISPN)+RHT(7,ISPN)
            end do
            do I=1,10
              RHT(I,1)= RHT(I,1)+RHT(I,NSPN)
            end do
         end if
      end if

!
! EXCHANGE
!
      D= DEN(1)
      if ((IDFTYP(1) == 0) .OR. (D < EPS)) GOTO 100

      do ISPN=1,NSPN
! Computing E[2Dup/2Ddn] below. S(2Dup), U(2Dup), and V(2Dup) are calculated
! For spin unpolarized system, E[D] is calculated and factors in S ,U ,V (2,4,8 etc..) will cancel out
         D= 2.0d0*DEN(ISPN+1)
         if (D > EPS) THEN
            DKF= (TRPI2*D)**THIRD
            if (IGGA(1) == 0) THEN
               S= 0.0D0
               U= 0.0D0
               V= 0.0D0
            else
               DGRAD= 2.0d0*SQRT(DG2(ISPN+1))
               DLAP= 2.0d0*DLP(ISPN+1)
               REC2= 0.5D0/(DKF*D)
               REC1= REC2*D
               REC22= REC2*REC2
               S= DGRAD*REC2
               V= DLAP*REC1*REC2
               if (DGRAD > EPS) THEN
                  U= 8.0d0*DGG(ISPN+1)*REC1*REC22/DGRAD
               else
                  U= 0.0D0
               end if
            end if

           !if (IDFTYP(1) .LE. 6) THEN
           !   CALL PW91EX(DKF,S,U,V,EX,DEX,VX)
           !else if (IDFTYP(1) == 7) THEN
           !   CALL PBEEX(D,S,U,V,1,1,EX,DEX,VX)
           !end if
           !LB: ADDED REVPBE AND RPBE

! YY: This is where GGA functional subroutines are called
            SELECT case(IDFTYP(1))
            case(:6)
               CALL PW91EX(DKF,S,U,V,EX,DEX,VX)    !PW91
            case(7)
               CALL PBEEX(1,D,S,U,V,1,1,EX,DEX,VX) !PBE
            case(8)
               CALL PBEEX(2,D,S,U,V,1,1,EX,DEX,VX) !REVPBE
            case(9)
               CALL PBEEX(3,D,S,U,V,1,1,EX,DEX,VX) !RPBE
            case(10)
               CALL XBECKE(D,S,U,V,EX,DEX,VX)      !BECKE 88
            case(11:15)
               CALL call_mgga_ex(IPTS,MPTS,ISPN,D,DGRAD, &
                2.0D0*tauchop(IPTS,ISPN),RHT(1:4,1),RHT(1:4,NSPN),EX,VX,LDFTF)
               DEX = 0.0d0
            END SELECT
! LB: END OF NEW CODE


            JPTS= IPTS+(ISPN-1)*MPTS
            if((.not.LDFTF).and.scaledsic) then
             if(SDSICON)then
              !SDSIC Energy is scaled in apotnl
              VXCS(JPTS)= VXCS(JPTS)+VX
              EXL=EXL+ISPFAC*0.5D0*EX *D*WMSH(LPTS+IPTS)
              EXN=EXN+ISPFAC*0.5D0*DEX*D*WMSH(LPTS+IPTS)
             else if(GSICON)then
              VXCS(JPTS)= VXCS(JPTS)+VX
              EXL=EXL+ISPFAC*0.5D0*EX *D*WMSH(LPTS+IPTS)*BETA(LPTS+IPTS)
              EXN=EXN+ISPFAC*0.5D0*DEX*D*WMSH(LPTS+IPTS)*BETA(LPTS+IPTS)
             else
              !Implementation of VXCS with scaling factor
              VXCS(JPTS)= VXCS(JPTS)+VX*TAUW_TAU(LPTS+IPTS)
              EXL=EXL+ISPFAC*0.5D0*EX *D*WMSH(LPTS+IPTS)*TAUW_TAU(LPTS+IPTS)
              EXN=EXN+ISPFAC*0.5D0*DEX*D*WMSH(LPTS+IPTS)*TAUW_TAU(LPTS+IPTS)
             endif
             SICEXCS(IPTS)= SICEXCS(IPTS)+ ISPFAC*0.5D0*EX*D
            else
             VXCS(JPTS)= VXCS(JPTS)+VX  ! add returned Vx
             EXL=EXL+ISPFAC*0.5D0*EX *D*WMSH(LPTS+IPTS)   ! Add returned local energy Ex
             EXN=EXN+ISPFAC*0.5D0*DEX*D*WMSH(LPTS+IPTS)   ! Add returned non local Energy
            end if
         end if
      end do
!
! CORRELATION
!
      100 continue

      D=   DEN(1)
      DUP= DEN(2)
      DDN= DEN(3)
      if ((IDFTYP(2) == 0) .OR. (D < EPS)) GOTO 200
      DREC= 1.0D0/D
      ZET= (DUP-DDN)*DREC  !Relative spin polarization (rhoup-rhodn)/rho
      DKF= (TRPI2*D)**THIRD
      DKS= SQRT(FRDPI*DKF)
      RS= DKF2RS/DKF
      G= 0.5D0*((1.0D0+ZET)**THIRD2+(1.0D0-ZET)**THIRD2)

      !A fix for NaN issue on Cori 
      if(ZET .EQ. -1.0D0) G= 0.5D0*(2.0D0**THIRD2)
      if(ZET .EQ.  1.0D0) G= 0.5D0*(2.0D0**THIRD2)

      if (IGGA(2) == 0) THEN
         T= 0.0D0
         UU= 0.0D0
         VV= 0.0D0
         WW= 0.0D0
      else
         DGRAD= SQRT(DG2(1))
         DLAP= DLP(1)
         REC3= DREC
         REC1= 0.5D0/(DKS*G)
         REC2= REC1*REC3
         REC22= REC2*REC2
         T= DGRAD*REC2
         VV= DLAP*REC1*REC2

! NOTE THAT GRAD(D)xGRAD(ZETA)*(D**2) IS EQUAL TO
! D*(GRAD(DUP)**2-GRAD(DDN)**2)-(DUP-DDN)*(GRAD(D)**2)

         WW= D*(DG2(2)-DG2(3))
         WW= WW-(DUP-DDN)*DG2(1)
         WW= WW*REC3*REC22
         if (DGRAD > EPS) then
            UU= DGG(1)*REC1*REC22/DGRAD
         else
           UU= 0.0D0
         end if
      end if

      if (IDFTYP(2) <= 5) then
! YY:LDA subroutine
         CALL LDACOR(D,ZET,IDFTYP(2),EC,VCOR)

!YY. for LDA, set DEC, DVCUP, DVCDN to zeroes (nonlocal)
         VCUP= VCOR(1)
         VCDN= VCOR(2)
         DEC=   0.0D0
         DVCUP= 0.0D0
         DVCDN= 0.0D0
      else if (IDFTYP(2) == 6) THEN
         CALL PW91LC(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
         CALL PW91NC(RS,ZET,T,UU,VV,WW,DEC,DVCUP,DVCDN)
      else if (IDFTYP(2) == 7) THEN
         CALL PBECOR(RS,ZET,T,UU,VV,WW,1,1,EC,VCUP,VCDN,DEC,DVCUP,DVCDN)

      else if (IDFTYP(2) >= 11) THEN
         CALL call_mgga_cor(IPTS,MPTS,DUP,DDN,DG2(1:3),RHT(1:4,1),RHT(1:4,NSPN),&
                 tauchop(IPTS,1),tauchop(IPTS,NSPN),EC,VCUP,VCDN,LDFTF)
         DEC=0.0d0
         DVCUP=0.0d0
         DVCDN=0.0d0

      end if

      if((.not.LDFTF).and.(scaledsic.and.(.not.GSICON))) then
       if(SDSICON) then 
        !Scaling is done in apotnl_sic in SDSIC
        ECL= ECL+ EC*D*WMSH(LPTS+IPTS)
        ECN= ECN+DEC*D*WMSH(LPTS+IPTS)
        VXCS(IPTS)= VXCS(IPTS)+(VCUP+DVCUP)
       else
        ECL= ECL+ EC*D*WMSH(LPTS+IPTS)*TAUW_TAU(LPTS+IPTS)
        ECN= ECN+DEC*D*WMSH(LPTS+IPTS)*TAUW_TAU(LPTS+IPTS)
        !Implementation of scaled FLOSIC
        VXCS(IPTS)= VXCS(IPTS)+(VCUP+DVCUP)*TAUW_TAU(LPTS+IPTS)
        ! This will multiply both exchange and correlation by tauw/tau here.
        if(ISMGGA) then
         mixins(1,IPTS)=mixins(1,IPTS)*TAUW_TAU(LPTS+IPTS)
         mixins(2,IPTS)=mixins(2,IPTS)*TAUW_TAU(LPTS+IPTS)
         mixins(3,IPTS)=mixins(3,IPTS)*TAUW_TAU(LPTS+IPTS)
         mixins(4,IPTS)=mixins(4,IPTS)*TAUW_TAU(LPTS+IPTS)
        end if
       end if
       SICEXCS(IPTS)=SICEXCS(IPTS) + EC*D

       if (NSPN == 2) then
        if(SDSICON) then
         VXCS(IPTS+MPTS)= VXCS(IPTS+MPTS) + VCDN+DVCDN
        else
         VXCS(IPTS+MPTS)= VXCS(IPTS+MPTS) +(VCDN+DVCDN)*TAUW_TAU(LPTS+IPTS)
         if(ISMGGA) then
          mixins(1,IPTS+MPTS)=mixins(1,IPTS+MPTS)*TAUW_TAU(LPTS+IPTS)
          mixins(2,IPTS+MPTS)=mixins(2,IPTS+MPTS)*TAUW_TAU(LPTS+IPTS)
          mixins(3,IPTS+MPTS)=mixins(3,IPTS+MPTS)*TAUW_TAU(LPTS+IPTS)
          mixins(4,IPTS+MPTS)=mixins(4,IPTS+MPTS)*TAUW_TAU(LPTS+IPTS)
         end if
        end if
       end if
      else
       ECL= ECL+ EC*D*WMSH(LPTS+IPTS)    !YY. add returned Ec here
       ECN= ECN+DEC*D*WMSH(LPTS+IPTS)
       VXCS(IPTS)= VXCS(IPTS)+VCUP+DVCUP !YY. add returned Vc up here
       if (NSPN == 2) VXCS(IPTS+MPTS)= VXCS(IPTS+MPTS)+VCDN+DVCDN
      end if

   200  CONTINUE

else if(libxc1) then
!######################################################
! Case 2 - Turn LIBXC on.                             #
!######################################################
   print *, "The code is not compiled with LIBXC. Terminating."
   call stopit
else
   write(*,*) "SUBVLXC:Error Unsupported LibxcON value."
end if
!500 continue

EXCVEC(1)=EXL
EXCVEC(2)=EXN
EXCVEC(3)=ECL
EXCVEC(4)=ECN

!DEALLOCATE LOCAL ARRAYS

deallocate(RHOC,STAT=IERR); if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR DEALLOCATING RHOC'
deallocate(XTMP,STAT=IERR); if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR DEALLOCATING XTMP'
deallocate(DTMP,STAT=IERR); if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR DEALLOCATING DTMP'
deallocate(RTMP,STAT=IERR); if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR DEALLOCATING RTMP'
deallocate(VLOC,STAT=IERR); if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR DEALLOCATING VLOC'
deallocate(TAUCHOP,STAT=IERR); if(IERR/=0)WRITE(6,*)'SUBVLXC:ERROR DEALLOCATING TAUCHOP'

return
end
