! UTEP Electronic Structure Lab (2020)
!3/30/2017
!This set of subroutines will create CLUSTER.wfx files at post-SCF processing
!Many parts of codes are referenced from MOLDEN2AIM source code and modified to
!work with our code.
!Reference, Molden2AIM: https://github.com/zorkzou/Molden2AIM
subroutine wfxdrv

use xmol, only : NUM_ATMS,XMOL_LIST,GET_LETTER,AU2ANG
use common2, only : ETOTAL, EKINONL, ENONLO, NSPN, CHNET
use common5, only : NWF
integer :: unitv
integer :: INUC
logical :: prtspn
character*2 :: LETTER
character*10 :: strnucname,strno
character*100 :: ctmp
real*8, allocatable :: cfmo(:),scalmo(:),conf(:)
integer, allocatable :: icmo(:)
!real*8 :: AU2ANG = 0.529177249D0
unitv=1732
igto=991
imol=992
ispn = 9994
iconf=9993

print *,"CREATING CLUSTER.wfx FILE"

!creates GTO.TMP file
call wfxmos
!create MOL.TMP file
call printmo(nmo,ntote)

open(unitv,file="CLUSTER.wfx")
rewind(unitv)

! Title
call wfxtag(unitv,.true.,"Title")
write(unitv,"(' NRLMOL 2017 ')")
call wfxtag(unitv,.false.,"Title")

! Keywords
call wfxtag(unitv,.true.,"Keywords")
write(unitv,"(' GTO')")
call wfxtag(unitv,.false.,"Keywords")

! Number of Nuclei
call wfxtag(unitv,.true.,"Number of Nuclei")
write(unitv,"(i8)") NUM_ATMS !nat
call wfxtag(unitv,.false.,"Number of Nuclei")

! Number of Occupied Molecular Orbitals
call wfxtag(unitv,.true.,"Number of Occupied Molecular Orbitals")
write(unitv,"(i8)")nmo
call wfxtag(unitv,.false.,"Number of Occupied Molecular Orbitals")

! Number of Perturbations
call wfxtag(unitv,.true.,"Number of Perturbations")
write(unitv,"(i8)")0
call wfxtag(unitv,.false.,"Number of Perturbations")

! Net Charge
call wfxtag(unitv,.true.,"Net Charge")
write(unitv,"(i8)") nint(CHNET)
! write(unitv,"(i8)")nint(chanet)
call wfxtag(unitv,.false.,"Net Charge")

! Number of Electrons (Core Electrons by ECP are excluded)
call wfxtag(unitv,.true.,"Number of Electrons")
write(unitv,"(i8)")ntote
call wfxtag(unitv,.false.,"Number of Electrons")

!prtspn=.false.
! if(ifspin .eq. 1 .and. ifbeta .eq. 1) then
elea = 0.d0
eleb = 0.d0
open(ispn,FILE='SPN.TMP',FORM='FORMATTED',status='OLD')
rewind(ispn)
do i=1,nmo
  read(ispn,"(a42)")ctmp
  read(ctmp(1:2),*)ix
  read(ctmp(3:22),*)x
  if(ix .eq. 1) then
    elea = elea + x !*dble(ifc4)
  else if(ix .eq. 2) then
    eleb = eleb + x !*dble(ifc4)
  end if
end do
if(NSPN .eq. 1) then
  elea = elea/2.0d0
  eleb = elea
end if
!  if(abs(elea+eleb-dble(ntote)) .lt. 0.01d0 .and. &
!     abs(elea-dble(nint(elea))) .lt. 0.01d0 .and. &
!     abs(eleb-dble(nint(eleb))) .lt. 0.01d0) prtspn = .true.
! end if
close(ispn)

! Number of Alpha Electrons
call wfxtag(unitv,.true.,"Number of Alpha Electrons")
! if(prtspn) then
write(unitv,"(i8)")nint(elea)
! else
!  write(unitv,"(' UNKNOWN')")
! end if
call wfxtag(unitv,.false.,"Number of Alpha Electrons")

! Number of Beta Electrons
call wfxtag(unitv,.true.,"Number of Beta Electrons")
! if(prtspn) then
write(unitv,"(i8)")nint(eleb)
! else
!  write(unitv,"(' UNKNOWN')")
! end if
call wfxtag(unitv,.false.,"Number of Beta Electrons")

! Electronic Spin Multiplicity (optional)
call wfxtag(unitv,.true.,"Electronic Spin Multiplicity")
! if(prtspn) then
ix = nint(abs(elea-eleb)) + 1
write(unitv,"(i8)") ix
! else
!  write(unitv,"(' UNKNOWN')")
! end if
call wfxtag(unitv,.false.,"Electronic Spin Multiplicity")

! Number of Core Electrons. Required if ECP is used. If so, nuclear charges need to be subtracted by core electron.
call wfxtag(unitv,.true.,"Number of Core Electrons")
iecp=0
 write(unitv,"(i8)")iecp       !variable not assigned(?)
! write(unitv,"(' UNKNOWN')")
call wfxtag(unitv,.false.,"Number of Core Electrons")

! Nuclear Names
call wfxtag(unitv,.true.,"Nuclear Names")
do INUC=1,NUM_ATMS
  call GET_LETTER(XMOL_LIST(INUC)%ANUM,LETTER)
  write(strno,'(I4)') INUC
  strnucname=adjustl(trim(LETTER))//adjustl(trim(strno))
  write(unitv,'(A10)') strnucname
end do
call wfxtag(unitv,.false.,"Nuclear Names")

! Atomic Numbers (Z)
call wfxtag(unitv,.true.,"Atomic Numbers")
do INUC=1,NUM_ATMS
  write(unitv,"(i8)") XMOL_LIST(INUC)%ANUM !ix
end do
call wfxtag(unitv,.false.,"Atomic Numbers")

! Nuclear Charges (Z-#core)
call wfxtag(unitv,.true.,"Nuclear Charges")
! rewind(iatm)
! rewind(icor)
! read(iatm,*)
do INUC=1,NUM_ATMS !nat
!  read(iatm,*)ctmp,j,ix
!  if(iecp .gt. 0) read(icor,*) iz, ncore
  write(unitv,"(e21.12e3)")dble(XMOL_LIST(INUC)%ANUM)
!  write(unitv,"(e21.12e3)")dble(ix-ncore)
end do
call wfxtag(unitv,.false.,"Nuclear Charges")

! Nuclear Cartesian Coordinates (in a.u.)
call wfxtag(unitv,.true.,"Nuclear Cartesian Coordinates")
do INUC=1,NUM_ATMS
  write(unitv,"(3e21.12e3)") XMOL_LIST(INUC)%RX, &
                             XMOL_LIST(INUC)%RY, &
                             XMOL_LIST(INUC)%RZ
end do
call wfxtag(unitv,.false.,"Nuclear Cartesian Coordinates")

! Number of Primitives
! call wfxtag(unitv,.true.,"Number of Primitives")
! write(unitv,"(i8)")ncar
! call wfxtag(unitv,.false.,"Number of Primitives")

! print basis functions
! file unit, wfx mode, gto file unit

call printbasis(unitv,ncar1,ncar2)

! <<<<<<<<<<<<<<<<<< ECP >>>>>>>>>>>>>>>>>> Do we need it?
! if(iecp .gt. 0)then
!   call wfxtag(unitv,.true.,"Additional Electron Density Function (EDF)")
! ! Number of EDF Primitives
!   call wfxtag(unitv,.true.,"Number of EDF Primitives")
!   write(unitv,"(i8)")nedf
!   call wfxtag(unitv,.false.,"Number of EDF Primitives")
! ! EDF Primitive Centers
!   call wfxtag(unitv,.true.,"EDF Primitive Centers")
! !  rewind(iedf)
!   do i=1,nedf
! !   read(iedf,"(i5)") ix
!    write(unitv,"(i8)",advance='no') ix
!    if(mod(i,5) .eq. 0) write(unitv,*)
!   end do
!   if(mod(nedf,5) .ne. 0) write(unitv,*)
!   call wfxtag(unitv,.false.,"EDF Primitive Centers")
! ! EDF Primitive Types
!   call wfxtag(unitv,.true.,"EDF Primitive Types")
!   write(unitv,"(5i8)") (1, i=1, nedf)
!   call wfxtag(unitv,.false.,"EDF Primitive Types")
! ! EDF Primitive Exponents
!   call wfxtag(unitv,.true.,"EDF Primitive Exponents")
! !  rewind(iedf)
!   do i=1,nedf
! !   read(iedf,"(5x,e21.12e3)") x
!    write(unitv,"(e21.12e3)",advance='no') x
!    if(mod(i,5) .eq. 0) write(unitv,*)
!   end do
!   if(mod(nedf,5) .ne. 0) write(unitv,*)
!   call wfxtag(unitv,.false.,"EDF Primitive Exponents")
! ! EDF Primitive Coefficients
!   call wfxtag(unitv,.true.,"EDF Primitive Coefficients")
!   rewind(iedf)
!   do i=1,nedf
!    read(iedf,"(26x,e21.12e3)") x
!    write(unitv,"(e21.12e3)",advance='no') x
!    if(mod(i,5) .eq. 0) write(unitv,*)
!   end do
!   if(mod(nedf,5) .ne. 0) write(unitv,*)
!   call wfxtag(unitv,.false.,"EDF Primitive Coefficients")
!   call wfxtag(unitv,.false.,"Additional Electron Density Function (EDF)")
! end if

! Molecular Orbital Occupation Numbers
call wfxtag(unitv,.true.,"Molecular Orbital Occupation Numbers")
ispn=44
open(ispn,File='SPN.TMP',FORM='FORMATTED',status='OLD')
rewind(ispn)
do i=1,nmo
  read(ispn,"(a42)")ctmp
  read(ctmp(3:22),*)x
  write(unitv,"(e21.12e3)")x !*dble(ifc4)
end do
call wfxtag(unitv,.false.,"Molecular Orbital Occupation Numbers")

! Molecular Orbital Energies
call wfxtag(unitv,.true.,"Molecular Orbital Energies")
!
rewind(ispn)
do i=1,nmo
  read(ispn,"(a42)")ctmp
  read(ctmp(23:42),*)x
  write(unitv,"(e21.12e3)")x
end do
call wfxtag(unitv,.false.,"Molecular Orbital Energies")
call wfxtag(unitv,.false.,"Molecular Orbital Energies")

! Molecular Orbital Spin Types
call wfxtag(unitv,.true.,"Molecular Orbital Spin Types")
rewind(ispn)
do i=1,nmo
!  if(prtspn) then
  if( NSPN .eq. 2) then
    read(ispn,"(a42)")ctmp
    read(ctmp(1:2),*)ix
    if(ix .eq. 1) then
      write(unitv,"(' Alpha')")
    else if(ix .eq. 2) then
      write(unitv,"(' Beta')")
    end if
  else
    write(unitv,"(' Alpha and Beta')")
  end if
end do
call wfxtag(unitv,.false.,"Molecular Orbital Spin Types")
!close(ispn) !YY moved down
! normalization factor
!  do i=1,ncar
!   cn(i)=fnorm_lmn(expg(i),ityp(i))
!  end do

! MO
call wfxtag(unitv,.true.,"Molecular Orbital Primitive Coefficients")
open(imol,FILE='MOL.TMP',FORM='FORMATTED',status='OLD')
rewind(imol)
open(iconf,FILE='CONF.TMP',FORM='FORMATTED',status='OLD')
rewind(iconf)
allocate(conf(ncar1))
allocate(icmo(ncar1))
read(iconf,*) ipt2
read(iconf,*) conf(1:ipt2)
read(iconf,*)  icmo(1:ipt2)

print *,nmo, NWF, ncar2, ipt2, ncar1
if(ipt2 .ne. ncar1) print *, "ipt2 does not match number &
                                of primitive functions"
allocate(cfmo(ncar2))
ix=0
do i=1,NWF !nmo !nmotot
  read(imol,*)ctmp,occ   !read OCCUP
  read(imol,*)           !read ENE
  do j=1,ncar2 !ncarc
    read(imol,*)ingc,cfmo(j)
!    cfmo(j)=cfmo(j)*scalmo(j)
  end do
  tolocc = 5.0d-5 !YY just because
  if(abs(occ) .ge. tolocc)then
    ix = ix + 1
! MO Number
    call wfxtag(unitv,.true.,"MO Number")
    write(unitv,"(i8)")ix
    call wfxtag(unitv,.false.,"MO Number")
!---      N*cgto*cmo
    !write(unitv,"(5e21.12e3)") (cn(j)*conf(j)*cfmo(icmo(j)),j=1,ncar)
    write(unitv,"(5e21.12e3)") (conf(j)*cfmo(icmo(j)),j=1,ncar1)
  end if
end do

deallocate(conf)
deallocate(cfmo)
close(ispn,status='delete')
close(iconf,status='delete')
close(imol,status='delete')
!close(ispn)
!close(iconf)
!close(imol)
call wfxtag(unitv,.false.,"Molecular Orbital Primitive Coefficients")

! Energy = T + Vne + Vee + Vnn
call wfxtag(unitv,.true.,"Energy = T + Vne + Vee + Vnn")
!  write(unitv,"(' UNKNOWN')")         !DFT Energy
write(unitv,"(e21.12)") ETOTAL       !DFT Energy
call wfxtag(unitv,.false.,"Energy = T + Vne + Vee + Vnn")

! Virial Ratio (-V/T)
call wfxtag(unitv,.true.,"Virial Ratio (-V/T)")
!write(unitv,"(' UNKNOWN')")
write(unitv,"(e21.12)") -1.0d0*(ETOTAL-EKINONL+ENONLO)/(EKINONL-ENONLO)
call wfxtag(unitv,.false.,"Virial Ratio (-V/T)")

close(98)
return
end subroutine wfxdrv
