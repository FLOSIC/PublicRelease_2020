! UTEP Electronic Structure Lab (2020)
subroutine write_fomesh(foidx)
 use debug1
 use mesh1,only : NMSH,RMSH,WMSH,FO_MESH
 implicit none
 integer,intent(in) :: foidx
 character*10 :: fotxt
 character*20 :: filename
 real*8,allocatable :: tmp_rmsh(:,:),tmp_wmsh(:)
 integer :: i,j,tally1,tally2

 write(fotxt,'(I0.3)') foidx
 filename='FOMSH'//trim(fotxt)
 tally2=0
 call tracer('NMSH',nmsh)
 do i=1,nmsh
!   write(6,*)'fo_mesh',i,fo_mesh(i)
   if(fo_mesh(i)) then
     tally2=tally2+1
   endif
 enddo
! call tracer('tally2',tally2)
 allocate(tmp_rmsh(3,tally2),stat=i)
 if(i/=0) write(6,*)'ERROR ALLOCATING TMP_RMSH'
 allocate(tmp_wmsh(tally2),stat=i)
 if(i/=0) write(6,*)'ERROR ALLOCATING TMP_WMSH'
 tally1=0
 do i=1,nmsh
   if(fo_mesh(i)) then
     tally1=tally1+1
     tmp_rmsh(:,tally1)=rmsh(:,tally1)
     tmp_wmsh(tally1)=wmsh(tally1)
   endif   
 enddo
 if(tally1/=tally2) write(6,*)'ERROR: TALLY1,TALLY2',tally1,tally2
 open(90,file=filename,form='unformatted',status='unknown')
 write(90) tally2
 write(6,*) 'Total reduced points',tally2,foidx
 write(90)((rmsh(j,i),j=1,3),i=1,tally1)
 write(90)(wmsh(i),i=1,tally1)
 close(90)
 deallocate(tmp_rmsh,stat=i)
 if(i/=0) write(6,*)'ERROR DEALLOCATING TMP_RMSH'
 deallocate(tmp_wmsh,stat=i)
 if(i/=0) write(6,*)'ERROR DEALLOCATING TMP_WMSH'
 return
end subroutine write_fomesh
