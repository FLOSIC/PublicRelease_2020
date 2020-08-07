! UTEP Electronic Structure Lab (2020)
! 11/6/2015
!YY. SUBVLXC called from ATOMSCF for initializing the potential 
!using LDA. This subroutine is used along side with subvlxc2.f90

! This subroutine will initialize the Vxc using LDACOR(for LDA)
! or PBECOR(for GGA) when libxc is in use.
!
! Don't forget to call check_input before calling this subroutine
!
! This subroutine used and modified from:
! 
! SUBVLXC DIRK POREZAG AUGUST 1999
! CALCULATES THE LOCAL AND EXCHANGE-CORRELATION POTENTIALS 
! FOR A SET OF POINTS
!
!LPTS=0, MPTS=1
subroutine ATOMSCFVLXC(RHOV,VXCS,EXCVEC,isGGA)
!SUBROUTINE SUBVLXC(MODE,LPTS,MPTS,RHOV,VXCS,VLOS,EXCVEC)

use mesh1,only : WMSH,RMSH
use common2,only : IDFTYP, ISPN, NSPN !,IGGA
!NSPN is 1 in atomscf call
use global_inputs, only : libxc1
!I can eliminate IGGA since atomscf is LDA
INCLUDE 'PARAMA2'
integer :: MXDFTYP,I,IERR,IOFS,IPTS,ISPFAC,JPTS,NGRAD
real*8 :: ALFC,AX,D,DDN,DEC,DEN,DKF,DKF2RS,DREC,DUP,DVCDN,DVCUP,&
          EC,ECL,ECN,ECRS,ECZET,EPS,EX,EXL,EXN,FAC,G,RHT,RS,    &
          THIRD,TRPI2,VCDN,VCOR,VCUP,VLOS,VX,ZET
SAVE
PARAMETER (MXDFTYP=10)
real*8, intent(in) :: RHOV(10*MXSPN*MPBLOCK)
real*8, intent(out) :: VXCS(MXSPN*MPBLOCK)
dimension VLOS(MPBLOCK)
real*8, intent(out) :: EXCVEC(4)
real*8, allocatable :: RHOC(:,:),XTMP(:,:),DTMP(:,:),RTMP(:)
integer :: LPTS=0
integer :: MPTS=1
integer, intent(in) :: isGGA(2)
!
! COMMON/PW91GAS/ IS NEEDED FOR PW91 GGA
!
COMMON/PW91GAS/G,EC,ECRS,ECZET
dimension RHT(10,MXSPN),VCOR(2),DEN(3)
data THIRD/0.3333333333333333D0/
!data THIRD2/0.6666666666666667D0/
data TRPI2/29.6088132032680759D0/
!data FRDPI/1.2732395447351627D0/ 
data DKF2RS/1.9191582926775128D0/
data EPS/1.0D-20/
data Ax/-0.738558766382022405884230032680836D0/
!
! allocate LOCAL ARRAYS
!
allocate(RHOC(10,MPBLOCK),STAT=ierr);  if(ierr/=0)write(6,*)'SUBVLXC:ERROR ALLOCATING RHOC'
allocate(XTMP(3,MPBLOCK),STAT=ierr);   if(ierr/=0)write(6,*)'SUBVLXC:ERROR ALLOCATING XTMP'
allocate(DTMP(3,MPBLOCK),STAT=ierr);   if(ierr/=0)write(6,*)'SUBVLXC:ERROR ALLOCATING DTMP'
allocate(RTMP(MPBLOCK),STAT=ierr);     if(ierr/=0)write(6,*)'SUBVLXC:ERROR ALLOCATING RTMP'
!
! INITIALIZATION OF DATA NEEDED FOR EXCHANGE-CORRELATION
!
EXL= 0.0D0
EXN= 0.0D0
ECL= 0.0D0
ECN= 0.0D0
ISPFAC= 2/NSPN
NGRAD=1
!
! GET CORE DENSITY
!
call GTRHOCR(.false.,MPTS,RMSH(1,1),RHOC,XTMP,DTMP,RTMP)
!CALL GTRHOCR(ISGGA,MPTS,RMSH(1,LPTS+1),RHOC,XTMP,DTMP,RTMP)
!
! LOOP OVER POINTS (although MPTS=1 here)
!
do 200 IPTS=1,MPTS
   VXCS(IPTS)= 0.0D0
   if (NSPN .EQ. 2) VXCS(IPTS+MPTS)= 0.0D0
!
! MOVE DATA TO RHT
!
   do ISPN=1,NSPN
      FAC= 1.0D0
      if (ISPN .EQ. 2) FAC= 0.5D0
      IOFS=NGRAD*((ISPN-1)+(IPTS-1)*NSPN)
      do I=1,NGRAD
         RHT(I,ISPN)= RHOV(IOFS+I)+FAC*RHOC(I,IPTS)
      end do
      if (RHT(1,ISPN) .LT. 0.0D0)    RHT(1,ISPN)= 0.0D0
      if (RHT(1,ISPN) .GT. RHT(1,1)) RHT(1,ISPN)= RHT(1,1)
   end do
!
! DENSITIES D, DUP, DDN
!
   DEN(1)= RHT(1,1)
   DEN(3)= DEN(1)*0.5D0
   if (NSPN .EQ. 2) DEN(3)= RHT(1,NSPN)
   DEN(2)= DEN(1)-DEN(3)
!
! EXCHANGE
!
   D= DEN(1)
   !Skip if the functional is set to GGA-NONE
   if ((IDFTYP(1) .EQ. 0) .OR. (D .LT. EPS)) GOTO 100
   do ISPN=1,NSPN
      D= 2*DEN(ISPN+1)
      if (D .GT. EPS) THEN

! Use LDA exchange energy and potential
         EX=Ax*D**THIRD
         !DEX=0.0d0
         VX=EX*1.33333333333333333d0

         !Ex and Vx above are equivalent of following call
         !call pbeex(1,D,0.0d0,0.0d0,0.0d0,1,1,EX,DEX,VX)

         JPTS= IPTS+(ISPN-1)*MPTS
         VXCS(JPTS)= VXCS(JPTS)+VX
         EXL=EXL+ISPFAC*0.5D0*EX *D*WMSH(LPTS+IPTS)
         !EXN=EXN+ISPFAC*0.5D0*DEX*D*WMSH(LPTS+IPTS)
      end if
   end do
!
! CORRELATION
!
   100 continue
   D=   DEN(1)
   DUP= DEN(2)
   DDN= DEN(3)
   !Skip if the functional is set to GGA-NONE
   if ((IDFTYP(2) .EQ. 0) .OR. (D .LT. EPS)) GOTO 200
   DREC= 1.0D0/D
   ZET= (DUP-DDN)*DREC
   DKF= (TRPI2*D)**THIRD
   !DKS= SQRT(FRDPI*DKF)
   RS= DKF2RS/DKF
   !G= 0.5D0*((1.0D0+ZET)**THIRD2+(1.0D0-ZET)**THIRD2)

! Check libxc and functional family
! Then initialize potential accordingly
   if(libxc1 .AND. isGGA(2) >= 1) then
      call PBECOR(RS,ZET,0.0d0,0.0d0,0.0d0,0.0d0,0,1,EC,VCUP,VCDN,DEC,DVCUP,DVCDN)
! Case for built-in meta-GGA functionals
   else if (isGGA(2) >= 3) then
      call PBECOR(RS,ZET,0.0d0,0.0d0,0.0d0,0.0d0,0,1,EC,VCUP,VCDN,DEC,DVCUP,DVCDN)
   else if(libxc1 .AND. isGGA(2) == 0) then
! Initialize LDA-WIGNER if functional is LIBXC-LDA
      call LDACOR(D,ZET,4,EC,VCOR)
      VCUP= VCOR(1)
      VCDN= VCOR(2)
! For builtin functionals, do following:
   else if (IDFTYP(2) .LE. 5) THEN
      call LDACOR(D,ZET,IDFTYP(2),EC,VCOR)
      VCUP= VCOR(1)
      VCDN= VCOR(2)
   else if (IDFTYP(2) .EQ. 6) THEN
      call PW91LC(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
      !call PW91NC(RS,ZET,T,UU,VV,WW,DEC,DVCUP,DVCDN)
   else if (IDFTYP(2) .EQ. 7) THEN
      !YY. set LGGA flag to zero
      call PBECOR(RS,ZET,0.0d0,0.0d0,0.0d0,0.0d0,0,1,EC,VCUP,VCDN,DEC,DVCUP,DVCDN)
      !call PBECOR(RS,ZET,T,UU,VV,WW,1,1,EC,VCUP,VCDN,DEC,DVCUP,DVCDN)
   end if
   !write(*,*) 'atomvlxc', EC,VCUP,libxc1

   ECL= ECL+ EC*D*WMSH(LPTS+IPTS)
   !ECN= ECN+DEC*D*WMSH(LPTS+IPTS)
   VXCS(IPTS)= VXCS(IPTS)+VCUP+DVCUP
   if (NSPN .EQ. 2) VXCS(IPTS+MPTS)= VXCS(IPTS+MPTS)+VCDN+DVCDN
200  continue
EXCVEC(1)=EXL 
EXCVEC(2)=EXN 
EXCVEC(3)=ECL 
EXCVEC(4)=ECN
 
!
! deallocate LOCAL ARRAYS
!
deallocate(RHOC,STAT=ierr);   if(ierr/=0)write(6,*)'SUBVLXC:ERROR DEALLOCATING RHOC'
deallocate(XTMP,STAT=ierr);   if(ierr/=0)write(6,*)'SUBVLXC:ERROR DEALLOCATING XTMP'
deallocate(DTMP,STAT=ierr);   if(ierr/=0)write(6,*)'SUBVLXC:ERROR DEALLOCATING DTMP'
deallocate(RTMP,STAT=ierr);   if(ierr/=0)write(6,*)'SUBVLXC:ERROR DEALLOCATING RTMP'

RETURN
END
