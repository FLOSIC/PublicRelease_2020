!> @file call_mgga.f90
!> @author Yoh Yamamoto
!> @brief subvlxc interface for calling meta-GGA functionals
!> @note 2021-7-20 UTEP ESLab   
!> @details This subroutine call a meta-GGA functional and returns
!>          Ex, Vx, and matrix elemets in MIXINS array
!>          MIXINS array is zeroed at subvlxc
subroutine call_mgga_ex(ipts,mpts,ISPN,D,DGRAD,TAU,RHTUP,RHTDN,EX,VX,LDFTF)

use common2,only  : IDFTYP, NSPN
!use xcmod,only   : MIXINS
use XTMP2A,only   : MIXINS
IMPLICIT NONE
logical, intent(in) :: LDFTF
!integer, intent(in) :: functional_id
integer, intent(in) :: ipts
integer, intent(in) :: mpts
integer, intent(in) :: ISPN
real(8), intent(in) :: D
real(8), intent(in) :: DGRAD
real(8), intent(in) :: TAU
real(8), intent(in) :: RHTUP(4) !copy of RHT(1:10,1) upto first deriv. terms.
real(8), intent(in) :: RHTDN(4) !copy of RHT(1:10,NSPN) upto first deriv. terms.

real(8), intent(out) :: EX
real(8), intent(out) :: VX

integer :: contraction ! =1 contracted grad used for libxc style functional
                       ! =0 grad used for original SCAN style functional
real(8) :: VXDD,AMUXD,DIV
real(8) :: RHT(4,NSPN)

RHT(:,1)=RHTUP(:)
RHT(:,NSPN)=RHTDN(:)

!### exchange ###

  contraction = 0
  !DEX = 0.0d0 !Do this outside

  select case(IDFTYP(1))
  case(11) !SCAN
    call VSCANxUNP(D,DGRAD,TAU,EX,VX,VXDD,AMUXD)
  case(12) !modSCAN used for GSIC
    if(LDFTF) then
     call VSCANxUNP(D,DGRAD,TAU,EX,VX,VXDD,AMUXD)
    else !(sic)
     call MODVSCANxUNP(D,DGRAD,TAU,EX,VX,VXDD,AMUXD)
    endif
  case(13) !Regularized SCAN
    call RSCANxUNP(D,DGRAD,TAU,EX,VX,VXDD,AMUXD)
  case(14) !r2SCAN 
    call xscan_r1(D,DGRAD,TAU,VX,VXDD,AMUXD,EX,2,1,1)
    contraction=1
  case(15) !pbe integration by parts
    call exchpbe(D,DGRAD,EX,VX,VXDD)
    !call stopit
    AMUXD=0.0d0
  end select

  !DIV here is the denominator term
  if(contraction==1) then
   DIV=1.0d0
  else
   DIV=DGRAD
  end if
    !r2SCAN uses contracted gradient like libxc and doesn't 
    !need to divide the elements by DGRAD as needed for dE/|dgrad|
 
  !common operations
  !prepare hamiltonian ingridients
  if(ABS(DIV) < 1.0D-100) DIV = 1.0d-100
  if(NSPN == 1) then
  !2's cancel out: DGRAD is twice the grad rho_up. 
  !RHT(2:4) is also twice the RHT_up. NSPN is 1. 
    mixins(1,IPTS)=VXDD*RHT(2,1)/DIV
    mixins(2,IPTS)=VXDD*RHT(3,1)/DIV
    mixins(3,IPTS)=VXDD*RHT(4,1)/DIV
    mixins(4,IPTS)=AMUXD*0.5d0

  else if(NSPN == 2) then
  !2's cancel out: DGRAD is twice the grad rho_up. NSPN is 2.
    if(ISPN == 1) then
      mixins(1,IPTS)=2.0d0*VXDD*(RHT(2,1)-RHT(2,2))/DIV
      mixins(2,IPTS)=2.0d0*VXDD*(RHT(3,1)-RHT(3,2))/DIV
      mixins(3,IPTS)=2.0d0*VXDD*(RHT(4,1)-RHT(4,2))/DIV
      mixins(4,IPTS)=AMUXD*0.5d0
    else if(ISPN == 2) then
      mixins(1,IPTS+MPTS)=2.0d0*VXDD*RHT(2,2)/DIV
      mixins(2,IPTS+MPTS)=2.0d0*VXDD*RHT(3,2)/DIV
      mixins(3,IPTS+MPTS)=2.0d0*VXDD*RHT(4,2)/DIV
      mixins(4,IPTS+MPTS)=AMUXD*0.5d0
    end if
  end if

return
end subroutine call_mgga_ex

!#################################################################################

!call call_mgga_cor(DUP,DDN,DG2(1:3),RHT(1:4,1:2),tauchop(IPTS,1),tauchop(IPTS,NSPN)
subroutine call_mgga_cor(ipts,mpts,DUP,DDN,DG2,RHTUP,RHTDN,tauup,taudn,EC,VCUP,VCDN,LDFTF)

use common2,only  : IDFTYP, NSPN
!use xcmod,only   : MIXINS
use XTMP2A,only   : MIXINS
IMPLICIT NONE
!integer,intent(in) :: functional_id
logical,intent(in) :: LDFTF
integer,intent(in) :: ipts,mpts
real(8),intent(in) :: DUP,DDN,DG2(3),RHTUP(4),RHTDN(4),tauup,taudn
real(8),intent(out) :: EC,VCUP,VCDN
integer :: contraction ! =1 contracted grad used for libxc style functional
                       ! =0 grad used for original SCAN style functional
                       ! =-1 do nothing
real(8) :: VCDD1,VCDD2, &
           AMUCD1,AMUCD2, &
           sig1,sig2,sig3, &
           vsigma1,vsigma2,vsigma3, &
           DIV(3)
real(8) :: RHT(4,NSPN)

RHT(:,1)=RHTUP(:)
RHT(:,NSPN)=RHTDN(:)

  contraction=0

  !in DUP,DDN,DG2(1:3),RHT(1:4,1:2), tauchop(IPTS,1:NSPN)
  !out EC,VCUP,VCDN
  select case(IDFTYP(2))
  case(11) ! SCAN
    call VSCANc(DUP,DDN,sqrt(DG2(2)),sqrt(DG2(3)),sqrt(DG2(1)),tauup,taudn,  &
                        EC,VCUP,VCDD1,VCDN,VCDD2,AMUCD1,AMUCD2)
  case(12) ! modSCAN
    if(LDFTF) then !Do the regular SCANc
     call VSCANc(DUP,DDN,sqrt(DG2(2)),sqrt(DG2(3)),sqrt(DG2(1)),tauup,taudn,  &
                         EC,VCUP,VCDD1,VCDN,VCDD2,AMUCD1,AMUCD2)
    else !(sic) Do nothing
     contraction=-1
     EC=0.0d0
     VCUP=0.0d0
     VCDN=0.0d0
    end if
  case(13) ! Regularized SCAN
    call RSCANc(DUP,DDN,sqrt(DG2(2)),sqrt(DG2(3)),sqrt(DG2(1)),tauup,taudn,  &
                        EC,VCUP,VCDD1,VCDN,VCDD2,AMUCD1,AMUCD2)
  case(14) ! r2SCAN
    !r2SCAN uses contracted grad like libxc. We prepare contracted grad here.
    contraction=1
    if(NSPN.EQ.1) then
     sig1 = RHT(2,1)**2 +RHT(3,1)**2 +RHT(4,1)**2
     sig1 = sig1/4.0d0
     sig2 = sig1
     sig3 = sig1
    else !NSPN.EQ.2
     sig1 = (RHT(2,1)-RHT(2,2))**2 &
          + (RHT(3,1)-RHT(3,2))**2 &
          + (RHT(4,1)-RHT(4,2))**2

     sig2 = (RHT(2,1)-RHT(2,2))*RHT(2,2) &
          + (RHT(3,1)-RHT(3,2))*RHT(3,2) &
          + (RHT(4,1)-RHT(4,2))*RHT(4,2)

     sig3 = RHT(2,2)**2 &
          + RHT(3,2)**2 &
          + RHT(4,2)**2
    end if
    CALL cscan_u1(DUP,DDN,sig1,sig2,sig3,tauup,taudn,  &
                   VCUP,VCDN,vsigma1,vsigma2,vsigma3,AMUCD1,AMUCD2,EC, &
                   2,1,1)
  case(15) !pbe integratin by parts
    call ecorpbe(DUP+DDN, sqrt(DG2(1)), (DUP-DDN)/(DUP+DDN), EC, VCUP, VCDN, VCDD1, NSPN)
    !subroutine ecorpbe(rho,agrad,zet,ectot,decup,decdn,decdg,nspin)
    VCDD2=VCDD1
    AMUCD1=0.0d0
    AMUCD2=0.0d0
  end select

  !prepare hamiltonian ingridients
  select case(contraction)
  case(0)
    !Case: grad
    DIV=DG2
    if(ABS(DIV(2)) < 1.0d-100) DIV(2) = 1.0d-100
    if(ABS(DIV(3)) < 1.0d-100) DIV(3) = 1.0d-100
    if(ABS(DIV(1)) < 1.0d-100) DIV(1) = 1.0d-100
    if(NSPN == 1) then
       !RHT(:,1)'s here are up and down combined. 2 in RHT and 4 in DG2(1) will cancel out
       mixins(1,IPTS)=mixins(1,IPTS)+VCDD1*RHT(2,1)/sqrt(DIV(1))
       mixins(2,IPTS)=mixins(2,IPTS)+VCDD1*RHT(3,1)/sqrt(DIV(1))
       mixins(3,IPTS)=mixins(3,IPTS)+VCDD1*RHT(4,1)/sqrt(DIV(1))
       mixins(4,IPTS)=mixins(4,IPTS)+AMUCD1*0.5d0
    else if(NSPN == 2) then
       mixins(1,IPTS)=mixins(1,IPTS)+VCDD1*(RHT(2,1)-RHT(2,2))/sqrt(DIV(2))
       mixins(2,IPTS)=mixins(2,IPTS)+VCDD1*(RHT(3,1)-RHT(3,2))/sqrt(DIV(2))
       mixins(3,IPTS)=mixins(3,IPTS)+VCDD1*(RHT(4,1)-RHT(4,2))/sqrt(DIV(2))
       mixins(4,IPTS)=mixins(4,IPTS)+AMUCD1*0.5d0

       mixins(1,IPTS+MPTS)=mixins(1,IPTS+MPTS)+VCDD2*RHT(2,2)/sqrt(DIV(3))
       mixins(2,IPTS+MPTS)=mixins(2,IPTS+MPTS)+VCDD2*RHT(3,2)/sqrt(DIV(3))
       mixins(3,IPTS+MPTS)=mixins(3,IPTS+MPTS)+VCDD2*RHT(4,2)/sqrt(DIV(3))
       mixins(4,IPTS+MPTS)=mixins(4,IPTS+MPTS)+AMUCD2*0.5d0
    end if
  case(1)
    !Case: contracted grad
    if(NSPN.EQ.1) then
      mixins(1,IPTS)=mixins(1,IPTS)+vsigma1*RHT(2,1)+vsigma2*RHT(2,1)*0.5d0
      mixins(2,IPTS)=mixins(2,IPTS)+vsigma1*RHT(3,1)+vsigma2*RHT(3,1)*0.5d0
      mixins(3,IPTS)=mixins(3,IPTS)+vsigma1*RHT(4,1)+vsigma2*RHT(4,1)*0.5d0
      mixins(4,IPTS)=mixins(4,IPTS)+AMUCD1*0.5d0
    else !NSPN.EQ.2
      mixins(1,IPTS)=mixins(1,IPTS)+2.0d0*vsigma1*(RHT(2,1)-RHT(2,2))+vsigma2*RHT(2,2)
      mixins(2,IPTS)=mixins(2,IPTS)+2.0d0*vsigma1*(RHT(3,1)-RHT(3,2))+vsigma2*RHT(3,2)
      mixins(3,IPTS)=mixins(3,IPTS)+2.0d0*vsigma1*(RHT(4,1)-RHT(4,2))+vsigma2*RHT(4,2)
      mixins(4,IPTS)=mixins(4,IPTS)+AMUCD1*0.5d0

      mixins(1,IPTS+MPTS)=mixins(1,IPTS+MPTS)+2.0d0*vsigma3*RHT(2,2)+vsigma2*(RHT(2,1)-RHT(2,2))
      mixins(2,IPTS+MPTS)=mixins(2,IPTS+MPTS)+2.0d0*vsigma3*RHT(3,2)+vsigma2*(RHT(3,1)-RHT(3,2))
      mixins(3,IPTS+MPTS)=mixins(3,IPTS+MPTS)+2.0d0*vsigma3*RHT(4,2)+vsigma2*(RHT(4,1)-RHT(4,2))
      mixins(4,IPTS+MPTS)=mixins(4,IPTS+MPTS)+AMUCD2*0.5d0
    end if
  case(-1)
    continue ! do nothing
  end select


  !Do this outside 
  !DEC=0.0d0
  !DVCUP=0.0d0
  !DVCDN=0.0d0

end subroutine call_mgga_cor
