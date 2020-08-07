! UTEP Electronic Structure Lab (2020)
subroutine allocate_diag1(mode)
 use for_diag1
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:00 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
 logical, intent(in) :: mode

 if(mode) then
  allocate(over(NDH,NDH))
  allocate(ham(NDH,NDH))
  allocate(filo(NDH,NDH))
  allocate(eval(NDH))
  allocate(sc1(NDH))
  allocate(sc2(NDH))
 else
  if(allocated(over)) deallocate(over)
  if(allocated(ham))  deallocate(ham)
  if(allocated(filo)) deallocate(filo)
  if(allocated(eval)) deallocate(eval)
  if(allocated(sc1))  deallocate(sc1)
  if(allocated(sc2))  deallocate(sc2)
 endif
 return
end
