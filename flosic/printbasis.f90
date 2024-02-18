! UTEP Electronic Structure Lab (2020)
!###############################################################################
subroutine printbasis(unitv,npg,npc)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
parameter(maxpg=14000,maxpgc=10000)
integer, intent(in) :: unitv
real*8, allocatable :: expg(:),conf(:),expgc(:),confc(:)
integer, allocatable :: icnt(:),ityp(:),icmo(:)
character*1 al
character*10 tmp
real :: dummyv
integer,intent(out) :: npg  !Number of primitive Cartesian functions
integer,intent(out) :: npc  !Number of contracted Cartesian functiions

npg = 0
npc = 0

allocate( expg(maxpg) )
allocate( conf(maxpg) )
allocate( expgc(maxpgc) )
allocate( confc(maxpgc) )
allocate( icnt(maxpg) )
allocate( ityp(maxpg) )
allocate( icmo(maxpg) )

igto = 9991
imol = 9992
iconf = 9993

open(unit=igto,file='GTO.TMP',form='FORMATTED',status='old')
rewind(igto)
open(unit=imol,file='MOL.TMP',form='FORMATTED',status='old')
rewind(imol)

ipt1=0
ipt2=0

read(imol,*)tmp,occ
read(imol,*)tmp,eng
100   read(igto,*,end=500)ic1
200   read(igto,*)al,np

if(al.eq.'E'.and.np.eq.0)goto 100
do i=1,np
 read(igto,*)expgc(i),confc(i)
 !print *,expgc(i),confc(i)
end do
select case(al)
 case('S')
  nfun=1
  it=0
 case('P')
  nfun=3
  it=1
 case('D')
  nfun=6
  it=4
 ! case('F')
 !  nfun=10
 !  it=10
 ! case('G')
 !  nfun=15
 !  it=20
end select
npg=npg+np*nfun
npc=npc+nfun
print *,npg
do i=1,nfun
 ipt1=ipt2+1
 ipt2=ipt1+np-1
 it=it+1
 expg(ipt1:ipt2)=expgc(1:np)
 conf(ipt1:ipt2)=confc(1:np)
 icnt(ipt1:ipt2)=ic1
 do j=ipt1,ipt2
!  if(it.ge.14.and.it.le.19)then
!! f: xyy, xxy, xxz, xzz, yzz, yyz (MOLDEN or Gaussian-out) -->
!!    xxy, xxz, yyz, xyy, xzz, yzz (Gaussian-WFN or GAMESS-WFN)
!   ityp(j)=iorder(it)
!  else
   ityp(j)=it
!  end if
 end do
!--- read MO coeff. index
 read(imol,*)igc
 icmo(ipt1:ipt2)=igc
end do

goto 200

500   continue

open(unit=iconf,file='CONF.TMP',form='FORMATTED',status='UNKNOWN')
rewind(iconf)
write(iconf,*) ipt2
write(iconf,*) conf(1:ipt2)
write(iconf,*) icmo(1:ipt2)
close(iconf)

! Number of Primitives
 call wfxtag(unitv,.true.,"Number of Primitives")
 write(unitv,"(i8)")npg
 call wfxtag(unitv,.false.,"Number of Primitives")

! writing to wfx
call wfxtag(unitv,.true.,"Primitive Centers")
write(unitv,"(5i20)")(icnt(i),i=1,npg)
call wfxtag(unitv,.false.,"Primitive Centers")
call wfxtag(unitv,.true.,"Primitive Types")
write(unitv,"(5i20)")(ityp(i),i=1,npg)
call wfxtag(unitv,.false.,"Primitive Types")
call wfxtag(unitv,.true.,"Primitive Exponents")
write(unitv,"(5e21.12e3)")(expg(i),i=1,npg)
call wfxtag(unitv,.false.,"Primitive Exponents")

close(imol)
close(igto)

deallocate( expg )
deallocate( conf )
deallocate( expgc )
deallocate( confc )
deallocate( icnt )
deallocate( ityp )
deallocate( icmo )

return
end subroutine printbasis
