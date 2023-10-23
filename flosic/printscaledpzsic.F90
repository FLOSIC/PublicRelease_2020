subroutine printscaledpzsic

use scaledpzsic

! This subroutine prints out the scaled down PZSIC configs.    
! The configs. of interest are as follows.
!
! logical ::  LSICON   = .TRUE.  !LSIC
! logical ::  ORBSICON = .FALSE. !OrbSIC
! logical ::  GSICON   = .FALSE. !GSIC
! logical ::  AVGSICON = .FALSE. !Average SIC potential
! logical ::  SDSICON  = .FALSE. !sdSIC

if(scaledsic) then
  print *,"Scaled PZSIC configs."
  if(LSICON) then
    print *,"Local SIC (LSIC) is used"
  else if(SDSICON) then
    print *,"sdSIC is used"
  else if(ORBSICON) then
    print *,"ORBSIC is used"
  else if(GSICON) then
    print *,"GSIC is used"
  else
    print *,"No scaling down applied"
  end if
endif 
    
if(oneshot) then
    print *,"Oneshot calculation"
else if(AVGSICON) then
    print *,"Averaged SIC potential"
else
    print *,"SCF calculation"
end if

#ifdef GROUP
  if(scaledsic.or.LSICON.or.SDSICON.or.ORBSICON.or.GSICON) then
    print *,"Scaled SIC is not supported for GROUP calculation"
    call stopit
  end if
#endif

end subroutine printscaledpzsic

