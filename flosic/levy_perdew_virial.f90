! UTEP Electronic Structure Lab (2020)
!  This routine computes the exchange energy using 
!  Levy-Perdew virial relation as in the following 
!  paper.
!  M. Levy and J. P. Perdew, “Hellmann–Feynman, virial, and scaling requisites
!  for the exact universal density functionals. shape of the correlation potential and
!  diamagnetic susceptibility for atoms”, Phys. Rev. A 32, 2010 (1985).
!  Rajendra Zope, April 25, 2018.
!  Call this routine in apotnl after exchange potential is calculated
!  through a call to getvlxc.

subroutine  ex_levy_perdew
use mixpot1,only : POT=>POTOUT
use mesh1,only : WMSH,RMSH,NMSH
use pot_dens,only : RHOG
use common2,only : NSPN
implicit none
real(8) :: aux,beriz1,beriz2,den
integer :: i
! Raja - April 2018 - Compute exchange energy using Levy-Perdew virial theorm

beriz1 = 0.d0
beriz2 = 0.d0
do i=1,nmsh
  beriz1 =  beriz1 +  wMSH(i)*rhog(i,1,1)
  beriz2 =  beriz2 +  wMSH(i)*rhog(i,1,2)
enddo
write(6,*) 'Integration of density 1', beriz1
write(6,*) 'Integration of density 2', beriz2

beriz1 = 0.d0
beriz2 = 0.d0
write(6,*) 'NSPN: is  ', nspn
if (nspn==1) then
  do i=1,nmsh 
    aux = 0.d0
    aux =       RMSH(1,i)*rhog(i,2,1)
    aux = aux + rmsh(2,i)*rhog(i,3,1)
    aux = aux + rmsh(3,i)*rhog(i,4,1)
    beriz1 = beriz1 + (pot(i)*(3*rhog(i,1,1) + aux)*wmsh(i))
  enddo 
  write(6,*)'Virial-Exnrgy- Levy-Perdew ', beriz1
else if (nspn==2) then
  do i=1,nmsh 
! First work on up spin

    aux = 0.d0
    aux =       RMSH(1,i)*(rhog(i,2,1)-rhog(i,2,2))
    aux = aux + rmsh(2,i)*(rhog(i,3,1)-rhog(i,3,2))
    aux = aux + rmsh(3,i)*(rhog(i,4,1)-rhog(i,4,2))
! Following is needed for spin-polarization case 
! For some reason when nspn=2, rhog(i,1,1) contains total density
! but  rhog(i,1,2) contains spin down density
    den = rhog(i,1,1) - rhog(i,1,2)
    beriz1 = beriz1 + (pot(i)*(3*den + aux)*wmsh(i))
! Now work on down spin
    aux = 0.d0
    aux =       RMSH(1,i)*rhog(i,2,2)
    aux = aux + rmsh(2,i)*rhog(i,3,2)
    aux = aux + rmsh(3,i)*rhog(i,4,2)
    beriz2 = beriz2+(pot(nmsh+i)*(3*rhog(i,1,2)+aux)*wmsh(i))
  enddo 
  write(6,*)'Virial-Exenergy- Levy-Perdew for upsspin', beriz1
  write(6,*)'Virial-Exenergy- Levy-Perdew for down spin', beriz2
  write(6,*)'Virial-Exnergy- Levy-Perdew total', beriz1+beriz2
!  write(6,*)'M: Virial-Exnergy- Levy-Perdew total',beriz1/2-beriz2
endif 

end subroutine ex_levy_perdew
