!> This is the new input reader


module inputvars
  ! this module contains routines to read a new input file for flosic 
  
  
  use global_inputs
  use mod_periodic_table
  use common2, only: natoms, chnet, spnnet, efield
!  use sicmat, only: lfrm
!  use frm, only: lfrm

!  include 'PARAMA2'


  implicit none
  ! use common2, only : natoms
!  integer :: natoms
  ! this in common2
!  real(dp) :: efield(3)
  !<LA: replace with 
!  integer :: lfrm(mxspn)
!  use common2, only chnet, spnnet
!  real(dp) :: chnet, spnnet

  !end of modules to add

!  character(300) :: cwd
!  logical :: average
!  real(dp) :: avg
!  real(dp) :: tmptre
  logical :: simplified_input
  integer :: num_occ, num_up, num_down
  character(60) :: functype
  character(5) :: symgrp
!  real(dp), allocatable :: geometry(:,:)
  real(8), allocatable :: geometry(:,:) 
  character(2), allocatable :: atoms(:) 
  ! atom integer values
  integer, allocatable :: atoms_int(:)
  character(6), allocatable :: spinpseudo(:) 

  contains

  subroutine set_initial_values()

!    use global_inputs
!    use mod_periodic_table
!    use common2

    ! this initilization is taken from check_inputs.f90
    ! it repeats some variables from init_inputs but preserved here
    ! as I am unsure if changing init_inputs will mess something up 

    idiag1 = 1
    idiag2 = 1
    idiag3 = 0
    dosjnt1 = .false.
    wfgrid1 = .false.
    dosoccu1 =.false.
    formfak1 =.false.
    atomsph1 =.false.
    matdipole1=.false.
    excited1=.false.
    nonscf1 =.false.
    dftd31  =.false.
    fragment1=.false.
    rhogrid1=.false.
    calc_basis=1
    solvent1=.false.
    maxscf=100
    scftol=1.0d-6
    molden1=.true.
!    cmolden=.false.
    nbo1=.true.
!    cnbo=.false.
    symmetry1=.true.
!    csymmetry=.false.
    mpi_io1=.false.
    spnorb1=.false.
    fixm1=.false.
    pcm1 =.false.
    wffrm1=.false.
    libxc1=.false.
    dmat1 =.false.
    mixing1=0
    nwfout=10
    spnpol1=.false.
    scanmesh1=0
    efp1=.false.
    population1=.false.
    scaledlbfgs1=.false.
    uniham1=.false.
    frozendensity1=.false.
    nonscfforces1=.false.
    fod_loop1 = .false.
    fod_opt2 = .false.
    fod_opt3 = 0

  end subroutine

  !> This subroutine reads the input file
  subroutine read_file(filename)

!    include 'PARAMA2'

!    use global_inputs
!    use mod_periodic_table
!    use common2
    
    character(len=*), intent(in) :: filename
    integer :: inputio, ios, ios2, frmorbio
    character(1) temp
    character(20) temp2
    character(300) iom
    integer :: i, j, k, l
    integer :: input_len, nrows
    character(1) :: line

    integer :: geo_start, geo_end, geom_lines
    integer :: frmorb_start, frmorb_end, xyz_end
    integer :: input_start, input_end, input_lines
    character(10) :: trash
!    real(dp), allocatable :: fodpos(:,:,:)
    real(8), allocatable :: fodpos(:,:,:) 
    integer :: lfrm2(2)
    character(50), allocatable :: keywargs(:)
    character(50), allocatable :: keywargs2(:)
    character(50), allocatable :: keywvalues(:)
    logical :: frmorb_exist, geo_exist, input_exist
    type(element_t) :: element

!    call init_inputs()
    call set_initial_values()

    ios = 0
    ios2 = 0
    ! openfile, stop if there is a problem
    open(newunit=inputio, file=trim(filename), status='old', action='read', iostat=ios, iomsg=iom)
    if (ios /= 0) then
      write(*,*) 'Fatal error!!!'
      write(*,*)  trim(iom)
      stop
    endif

    ! get the number of lines in a file w/o using a external bash call
    input_len = 0
    if (ios == 0) then
      do while (ios2 == 0)
        if (ios2 == 0) then
          read(inputio,*,iostat=ios2) temp
          input_len = input_len + 1
        endif
      enddo
    endif
    ! TODO error handling 
    input_len = input_len - 1
    rewind(inputio)

    nrows = input_len

    allocate(keywargs(nrows))

    ! read the first word on each row and rewind
    do i=1, nrows
      read(inputio,*) keywargs(i)
    enddo 
  
    rewind(inputio)

    geo_exist = .false.

    do i=1, nrows
      if (trim(keywargs(i)) == 'geometry') then
        geo_start = i
        geo_exist = .true. 
      else if (trim(keywargs(i)) == 'end') then
        if (geo_exist .eqv. .true.) then
          geo_end = i
          exit 
        endif 
      endif 
    enddo 

    if (geo_exist .eqv. .false.) then
      write(*,*) 'Geometry not present in input file!'
      write(*,*) 'Check input for typos '
      stop
    endif 

    do i=1, geo_start
      read(inputio,*) trash
    enddo 

    geom_lines = geo_end - geo_start - 1
    ! TODO ensure mistakes in input are caught
    !      have @file option to read a file instead of lines 

    ! first two lines should be symmetry, then charge and spin
    ! natoms will be the number of geometry lines minus the 2 for symm and charge
    natoms = geom_lines - 2
    write(*,*) 'Number of atoms detected:', natoms
    allocate(geometry(3,natoms))
    allocate(atoms(natoms))
    allocate(atoms_int(natoms))
    allocate(spinpseudo(natoms))
    ! read symmmetry then charge and spin 
    ! symgrp = symmetry 
    read(inputio,*) symgrp
    ! TODO incorporate alternative for pseudopotential option

    read(inputio,*) chnet, spnnet

    do i=1, natoms
      read(inputio,*) atoms(i), (geometry(j,i),j=1,3), spinpseudo(i) 
      element = get_element(atoms(i))
      atoms_int(i) = element%n_protons
    enddo 

    ! test
!    write(*,*) 'Input geometry'
!    do i=1, natoms
!      write(*,*) atoms(i), (geometry(j,i),j=1,3), spinpseudo(i)
!    enddo
!    write(*,*) 'Atom numbers'
!    do i=1, natoms
!      write(*,*) atoms_int(i)
!    enddo
    ! end test

    
    frmorb_exist = .false.
  
    do i=geo_end+1, nrows
      if (trim(keywargs(i)) == 'frmorb') then
        frmorb_exist = .true. 
        frmorb_start = i
      else if (trim(keywargs(i)) == 'end') then
        if (frmorb_exist .eqv. .true.) then
          frmorb_end = i
          exit
        endif
      endif
    enddo
    !TODO change frmorb to read top line for nspn up and down 

    if (frmorb_exist) then

      do i=geo_end, frmorb_start
        read(inputio,*) trash
      enddo

      read(inputio,*) num_up, num_down
      lfrm2 = [num_up,num_down]
!      num_occ = maxval([num_up,num_down])
      num_occ = maxval(lfrm2)
!      write(*,*) 'numocc is ', num_occ
      
      allocate(fodpos(3,num_occ,2))
!      frmorb_end - frmorb_start

      do i=1,2
        do k=1, lfrm2(i)
          read(inputio,*) (fodpos(j,k,i),j=1,3)
        enddo
      enddo

!      write(*,*) 'FOD positions'
      open(newunit=frmorbio, file='FRMORB', status='unknown', action='write')
      write(frmorbio,*) lfrm2
      do i=1,2
        do k=1, lfrm2(i)
          write(frmorbio,*) (fodpos(j,k,i),j=1,3)
        enddo
      enddo
      close(frmorbio)

      xyz_end = frmorb_end

    else

      xyz_end = geo_end

    endif 

    input_exist = .false.

    do i=xyz_end+1, nrows
      if (trim(keywargs(i)) == 'inputdata') then
        input_exist = .true.
        input_start = i
      else if (trim(keywargs(i)) == 'end') then
        if (input_exist .eqv. .true.) then
          input_end = i
          exit
        endif 
      endif
    enddo

    if (input_exist) then

      do i=xyz_end, input_start
        read(inputio,*) trash
      enddo

      input_lines = input_end - input_start - 1
      
      allocate(keywargs2(input_lines))
      allocate(keywvalues(input_lines))

      do i=1, input_lines
        read(inputio,*) keywargs2(i), trash, keywvalues(i)
        write(*,*) 'Read in: ', trim(keywargs2(i)), ' ', trim(keywvalues(i))
      enddo 

    else
      write(*,*) 'No input parameter section found, using default values'
    endif 

    rewind(inputio)

    ! check for efield 
    do i=1, nrows
      read(inputio,*) trash
      if (trim(keywargs(i)) == 'efield') then
        read(inputio,*) (efield(j),j=1,3)
        write(*,*) 'Found efield parameter'
        write(*,*) (efield(j),j=1,3)
        exit
      endif 
    enddo

    rewind(inputio)

    ! uncomment later
!    avg = 0.15_dp
!    average = .false.
!
!    ! check for avrgdat
!    do i=1, nrows
!      read(inputio,*) trash
!      if (trim(keywargs(i)) == 'avrgdat') then
!        read(inputio,*)  avg, average
!        exit
!      endif
!    enddo
!
!    rewind(inputio)
!
!    ! check for temperature
!    do i=1, nrows
!      read(inputio,*) trash
!      if (trim(keywargs(i)) == 'temperature') then
!        read(inputio,*)  tmptre
!        write(*,'((a,1x),(es16.7))') 'Found temperature:', tmptre
!        exit
!      endif
!    enddo
!
!    rewind(inputio)

    close(inputio)
    
    ! process the keyword arguments 
    ! implementation of new keywords will go here 
  
    if (input_exist) then
      do i=1, input_lines

        if (trim(keywargs2(i)) == 'scftolv' .or. trim(keywargs2(i)) == 'SCFTOLV') then
          call str2real(keywvalues(i), scftol)
        elseif (trim(keywargs2(i)) == 'FUNC' .or. trim(keywargs2(i)) == 'func') then
          functype = trim(keywvalues(i))
          write(*,*) functype

        elseif (trim(keywargs2(i)) == 'ATOMSPHV' .or. trim(keywargs2(i)) == 'atomsphv') then
          call read_logic(keywvalues(i), atomsph1)

        elseif (trim(keywargs2(i)) == 'DFTD3V' .or. trim(keywargs2(i)) == 'dftd3v') then
          call read_logic(keywvalues(i), dftd31)

        elseif (trim(keywargs2(i)) == 'JNTDOSV' .or. trim(keywargs2(i)) == 'jntdosv') then
          call read_logic(keywvalues(i), dosjnt1)

        elseif (trim(keywargs2(i)) == 'DMATV' .or. trim(keywargs2(i)) == 'dmatv') then
          call read_logic(keywvalues(i), dmat1)

        elseif (trim(keywargs2(i)) == 'DOSOCCUV' .or. trim(keywargs2(i)) == 'dosoccuv') then
          call read_logic(keywvalues(i), dosoccu1)

        elseif (trim(keywargs2(i)) == 'FIXMV' .or. trim(keywargs2(i)) == 'fixmv') then
          call read_logic(keywvalues(i), fixm1)

        elseif (trim(keywargs2(i)) == 'FOD_LOOPV' .or. trim(keywargs2(i)) == 'fod_loopv') then
          call read_logic(keywvalues(i), fod_loop1)

          ! mod here 
        elseif (trim(keywargs2(i)) == 'FOD_OPT1V' .or. trim(keywargs2(i)) == 'fod_opt1v') then
          if (trim(keywvalues(i)) == 'LBFGS' .or. trim(keywvalues(i)) == 'lbfgs') then
            fod_opt1 = .true.
          else
            fod_opt1 = .false.
          endif 

        elseif (trim(keywargs2(i)) == 'FOD_OPT2V' .or. trim(keywargs2(i)) == 'fod_opt2v') then
          call read_logic(keywvalues(i), fod_opt2)

        elseif (trim(keywargs2(i)) == 'FRAGMENTV' .or. trim(keywargs2(i)) == 'fragmentv') then
          call read_logic(keywvalues(i), fragment1)

        elseif (trim(keywargs2(i)) == 'NONSCFV' .or. trim(keywargs2(i)) == 'nonscfv') then
          call read_logic(keywvalues(i), nonscf1)

        elseif (trim(keywargs2(i))=='NONSCFFORCESV'.or.trim(keywargs2(i)) == 'nonscfforcesv') then
          call read_logic(keywvalues(i), nonscfforces1)

        elseif (trim(keywargs2(i)) == 'POPULATIONV' .or. trim(keywargs2(i)) == 'populationv') then
          call read_logic(keywvalues(i), population1)

        elseif (trim(keywargs2(i)) == 'RHOGRIDV' .or. trim(keywargs2(i)) == 'rhogridv') then
          call read_logic(keywvalues(i), rhogrid1)

        elseif (trim(keywargs2(i)) == 'SPNPOLV' .or. trim(keywargs2(i)) == 'spnpolv') then
          call read_logic(keywvalues(i), spnpol1)

        elseif (trim(keywargs2(i)) == 'SYMMETRYV' .or. trim(keywargs2(i)) == 'symmetryv') then
          call read_logic(keywvalues(i), symmetry1)

        elseif (trim(keywargs2(i)) == 'UNIHAMV' .or. trim(keywargs2(i)) == 'unihamv') then
          call read_logic(keywvalues(i), uniham1)

        elseif (trim(keywargs2(i)) == 'WFGRIDV' .or. trim(keywargs2(i)) == 'wfgridv') then
          call read_logic(keywvalues(i), wfgrid1)

        elseif (trim(keywargs2(i)) == 'WFFRMV' .or. trim(keywargs2(i)) == 'wffrmv') then
          call read_logic(keywvalues(i), wffrm1)

        ! modded here
        elseif (trim(keywargs2(i)) == 'BASISV' .or. trim(keywargs2(i)) == 'basisv') then
!          call str2int(keywvalues(i), calc_basis)
          call read_basistype(keywvalues(i),calc_basis)
        ! modded here
        elseif (trim(keywargs2(i)) == 'CALCTYPEV' .or. trim(keywargs2(i)) == 'calctypev') then
          call read_calctype(keywvalues(i), calctype1)

        elseif (trim(keywargs2(i)) == 'DIAG1V' .or. trim(keywargs2(i)) == 'diag1v') then
          call str2int(keywvalues(i), idiag1)

        elseif (trim(keywargs2(i)) == 'DIAG2V' .or. trim(keywargs2(i)) == 'diag2v') then
          call str2int(keywvalues(i), idiag2)

        elseif (trim(keywargs2(i)) == 'DIAG3V' .or. trim(keywargs2(i)) == 'diag3v') then
          call str2int(keywvalues(i), idiag3)

        elseif (trim(keywargs2(i)) == 'FOD_OPT3V' .or. trim(keywargs2(i)) == 'fod_opt3v') then
          call str2int(keywvalues(i), fod_opt3)
          if (fod_opt3 < 0 .or. fod_opt3 > 4) then
            fod_opt3 = 0 
          endif

        elseif (trim(keywargs2(i)) == 'MAXSCFV' .or. trim(keywargs2(i)) == 'maxscfv') then
          call str2int(keywvalues(i), maxscf)

        ! mod here
        elseif (trim(keywargs2(i)) == 'MIXINGV' .or. trim(keywargs2(i)) == 'mixingv') then
          if (keywvalues(i) == 'P' .or. keywvalues(i) == 'p') then 
            mixing1 = 0 
          elseif (keywvalues(i) == 'D' .or. keywvalues(i) == 'd') then 
            mixing1 = 1 
          elseif (keywvalues(i) == 'H' .or. keywvalues(i) == 'h') then 
            mixing1 = 2 
          else
            mixing1 = 0 
          endif
!          call str2int(keywvalues(i), mixing1)

        elseif (trim(keywargs2(i)) == 'NWFOUTV' .or. trim(keywargs2(i)) == 'nwfoutv') then
          call str2int(keywvalues(i), nwfout)

        elseif (trim(keywargs2(i)) == 'MESHSETTINGV'.or. trim(keywargs2(i)) == 'meshsettingv') then
          call str2int(keywvalues(i), scanmesh1)

        ! end of standard flosic keywords, new additions can go here 

        endif 
      enddo 
    endif 

    ! testing area
!    write(*,*)
!    write(*,*) 'TESTING'
!    write(*,*) 'basis found was: ', calc_basis
!    write(*,*) 'calctype is: ', calctype1

    return 

  end subroutine

  elemental subroutine str2int(str,intval)
    implicit none
    character(len=*), intent(in) :: str
    integer, intent(out) :: intval

    read(str,*) intval

  end subroutine str2int

  elemental subroutine str2real(str,realval)
    implicit none
    character(len=*), intent(in) :: str
    real(8), intent(out) :: realval
!    real(dp), intent(out) :: realval

    read(str,*) realval

  end subroutine str2real

  subroutine read_logic(str, logicval)
    implicit none
    character(len=*), intent(in) :: str
    logical, intent(inout) :: logicval

    if (trim(str) == 'yes' .or. trim(str) == 'y') then
      logicval = .true.
    elseif (trim(str) == 'T' .or. trim(str) == 't') then
      logicval = .true.
    elseif (trim(str) == 'Y') then
      logicval = .true.
    elseif (trim(str) == 'Yes' .or. trim(str) == 'YES') then
      logicval = .true.
    elseif (trim(str) == 'no' .or. trim(str) == 'n') then
      logicval = .false.
    elseif (trim(str) == 'N') then
      logicval = .false.
    elseif (trim(str) == 'F' .or. trim(str) == 'f') then
      logicval = .false.
    elseif (trim(str) == 'NO' .or. trim(str) == 'No') then
      logicval = .false.
    !  this
    else
      write(*,*) 'Input for '//trim(str)//' not recognized'
      write(*,*) 'Using default value of: ', logicval
    endif

  end subroutine read_logic

  subroutine read_calctype(ctype,calctype)
    implicit none
    character(len=*), intent(in) :: ctype
    integer, intent(inout) :: calctype

    if (trim(ctype) == 'LBFGS' .or. trim(ctype) == 'lbfgs') then
      calctype = 1
    elseif (trim(ctype) == 'SCF-ONLY' .or. trim(ctype) == 'scf-only') then
      calctype = 2
    elseif (trim(ctype) == 'VERLET' .or. trim(ctype) == 'verlet') then
      calctype = 3
    elseif (trim(ctype) == 'VIBRATIONAL' .or. trim(ctype) == 'vibrational') then
      calctype = 4
    elseif (trim(ctype) == 'REDUNDANT' .or. trim(ctype) == 'redundant') then
      calctype = 5
    else
      write(*,*) "Unrecognized value for calctypev: ", trim(ctype)
      write(*,*) "Assuming scf-only"
      calctype = 2
    endif

  end subroutine read_calctype

  subroutine read_basistype(basisv, calc_basis)
    implicit none
    character(len=*), intent(in) :: basisv
    integer, intent(inout) :: calc_basis

    if(basisv=='DEFAULT'.or.basisv=='default')then
      calc_basis=1
    elseif(basisv=='3-21G')then
       calc_basis=2
    elseif(basisv=='3-21Gs')then
       calc_basis=3
    elseif(basisv=='3-21GSP')then
       calc_basis=4
    elseif(basisv=='3-21Gs_Polarization')then
       calc_basis=5
    elseif(basisv=='3-21ppG')then
       calc_basis=6
    elseif(basisv=='3-21ppGs')then
       calc_basis=7
    elseif(basisv=='6-311G2df2pd')then
       calc_basis=8
    elseif(basisv=='6-311G')then
       calc_basis=9
! Number 10 was removed
    elseif(basisv=='6-311Gss')then
       calc_basis=11
    elseif(basisv=='6-311Gss_Polarization')then
       calc_basis=12
    elseif(basisv=='6-311pGs')then
       calc_basis=13
    elseif(basisv=='6-311ppG2d2p')then
       calc_basis=14
    elseif(basisv=='6-311ppG3df3pd')then
       calc_basis=15
    elseif(basisv=='6-311ppGss')then
       calc_basis=16
    elseif(basisv=='6-31G_3df3pd')then
       calc_basis=17
    elseif(basisv=='6-31G-Blaudeau')then
       calc_basis=18
    elseif(basisv=='6-31Gs')then
       calc_basis=19
    elseif(basisv=='6-31Gs-Blaudeau')then
       calc_basis=20
    elseif(basisv=='6-31Gs_Polarization')then
       calc_basis=21
    elseif(basisv=='6-31Gss')then
       calc_basis=22
    elseif(basisv=='6-31Gss_Polarization')then
       calc_basis=23
    elseif(basisv=='6-31pGs')then
       calc_basis=24
    elseif(basisv=='6-31ppG')then
       calc_basis=25
    elseif(basisv=='6-31ppGs')then
       calc_basis=26
    elseif(basisv=='6-31ppGss')then
       calc_basis=27
    elseif(basisv=='6-32G')then
       calc_basis=28
    elseif(basisv=='Ahlrichs_Coulomb_Fitting')then
       calc_basis=29
    elseif(basisv=='Ahlrichs_Polarization')then
       calc_basis=30
    elseif(basisv=='Ahlrichs_pVDZ')then
       calc_basis=31
    elseif(basisv=='Ahlrichs_TZV')then
       calc_basis=32
    elseif(basisv=='Ahlrichs_VDZ')then
       calc_basis=33
    elseif(basisv=='Ahlrichs_VTZ')then
       calc_basis=34
    elseif(basisv=='DZVP2')then
       calc_basis=35
    elseif(basisv=='DZVP')then
       calc_basis=36
    elseif(basisv=='GAMESS_PVTZ')then
       calc_basis=37
    elseif(basisv=='GAMESS_VTZ')then
       calc_basis=38
    elseif(basisv=='Huzinaga_polarization')then
       calc_basis=39
    elseif(basisv=='IGLO-II')then
       calc_basis=40
    elseif(basisv=='IGLO-III')then
       calc_basis=41
    elseif(basisv=='McLeanChandler_VTZ')then
       calc_basis=42
    elseif(basisv=='MIDI_Huzinaga')then
       calc_basis=43
    elseif(basisv=='MINI_Huzinaga')then
       calc_basis=44
    elseif(basisv=='Partridge_Uncontracted_1')then
       calc_basis=45
    elseif(basisv=='Partridge_Uncontracted_2')then
       calc_basis=46
    elseif(basisv=='Partridge_Uncontracted_3')then
       calc_basis=47
    elseif(basisv=='Partridge_Uncontracted_4')then
       calc_basis=48
    elseif(basisv=='Saddlej')then
       calc_basis=49
    elseif(basisv=='STO-2G')then
       calc_basis=50
    elseif(basisv=='STO-3G')then
       calc_basis=51
    elseif(basisv=='STO-3Gs')then
       calc_basis=52
    elseif(basisv=='STO-3Gs_Polarization')then
       calc_basis=53
    elseif(basisv=='STO-6G')then
       calc_basis=54
    elseif(basisv=='TZVP')then
       calc_basis=55
    elseif(basisv=='STUTTGART_RLC') then
       calc_basis=56
    elseif(basisv=='LANL2DZ') then
       calc_basis=57
    elseif(basisv=='LANL2DZDP') then
       calc_basis=58
    elseif(basisv=='UGBS') then
       calc_basis=59
    else
      write(*,*) 'Unknown value found for basisv, using default'
      calc_basis=1
    endif

  end subroutine

end module
