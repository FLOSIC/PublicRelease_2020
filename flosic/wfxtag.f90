! UTEP Electronic Structure Lab (2020)
!###############################################################################
subroutine wfxtag(unitv,mode,string)
logical mode
integer, intent(in):: unitv
character*(*) string

if(mode) then
  write(unitv,"('<',a,'>')")string
else
  write(unitv,"('</',a,'>')")string
end if

return
end subroutine wfxtag
