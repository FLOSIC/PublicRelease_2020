! UTEP Electronic Structure Lab (2020)
subroutine init_inputs
  use global_inputs
  use hstor1,only : MAXCLUSTERSIZE
  implicit none
! Set initial values for variables checking
! This is to avoid repetitive warning messages

  cexcited=.TRUE.
  cdosjnt=.TRUE.
  cwfgrid=.TRUE.
  cdosoccu=.TRUE.
  cformfak=.TRUE.
  catomsph=.TRUE.
  cmatdipole=.TRUE.
  cnonscf=.TRUE.
  cnonscfforces=.TRUE.
  cdftd3=.TRUE.
  ccalctype=.TRUE.
  cfragment=.TRUE.
  cbasis=.TRUE.
  cpolarizability=.TRUE.
  crhogrid=.TRUE.
  csolvent=.TRUE.
  cmaxscf=.TRUE.
  cscftol=.TRUE.
  cmolden=.TRUE.
  cnbo=.TRUE.
  csymmetry=.TRUE.
  cmpi_io=.TRUE.
  cspnorb=.TRUE.
  cfixm=.TRUE.

! Set initial default values for variables
! These values change whenever NRLMOL_INPUT.DAT is processed

  idiag1=1
  idiag2=1
  idiag3=0
  calc_basis=1
  dosjnt1=.FALSE.
  wfgrid1=.FALSE.
  dosoccu1=.FALSE.
  formfak1=.FALSE.
  atomsph1=.FALSE.
  matdipole1=.FALSE.
  excited1=.FALSE.
  nonscf1=.FALSE.
  dftd31=.FALSE.
  calctype1=1
  fragment1=.FALSE.
  polarizability1=.FALSE.
  rhogrid1=.FALSE.
  solvent1=.FALSE.
  molden1=.FALSE.
  nbo1=.FALSE.
  symmetry1=.FALSE.
  mpi_io1=.FALSE.
  spnorb1=.FALSE.
  fixm1=.FALSE.
  scanmesh1=0
  population1=.FALSE.
  scaledlbfgs1=.TRUE.
  uniham1=.FALSE.
! 
  MAXCLUSTERSIZE=1
end subroutine init_inputs
