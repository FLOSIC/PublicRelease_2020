! UTEP Electronic Structure Lab (2020)
subroutine read_fomesh(foidx)
 use mesh1,only : NMSH,RMSH,WMSH
 implicit none
 integer,intent(in) :: foidx
 integer :: i,j,dummy
 character*10 :: fotxt
 character*20 :: filename
 logical :: exist1

 if(foidx==-1) then
   open(90,file='VMOLD',form='unformatted',status='unknown')
   read(90) nmsh,dummy
   read(90) ((rmsh(j,i),j=1,3),i=1,nmsh)
   read(90)(wmsh(i),i=1,nmsh)
   close(90)
   write(6,*)'READ_FOMESH read origianl mesh',nmsh
 else
   write(fotxt,'(I0.3)') foidx
   filename='FOMSH'//trim(fotxt)
   inquire(file=filename,exist=exist1)
   if(exist1) then
     open(90,file=filename,form='unformatted',status='unknown')
     read(90) nmsh
     read(90) ((rmsh(j,i),j=1,3),i=1,nmsh)
     read(90)(wmsh(i),i=1,nmsh)
     close(90)
   else
     write(6,*)'ERROR: File does not exist:',filename
     call stopit
   endif
 endif
end subroutine read_fomesh
