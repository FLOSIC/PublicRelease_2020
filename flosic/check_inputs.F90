! UTEP Electronic Structure Lab (2020)
subroutine check_inputs
  use global_inputs
  use common2,only : IGGA,IDFTYP 
  use common2,only : OPTTYP
  use inputvars, only: simplified_input
  implicit none
  logical   :: exist,lopen
  integer   :: ios
  character*80 :: line,linebuff
  character*10 :: calctypev  = 'LBFGS'
  character*30 :: basisv     = 'DEFAULT'
  character :: jntdosv       = 'F'
  character :: wfgridv       = 'F'
  character :: dosoccuv      = 'F'
  character :: excitedv      = 'F'
  character :: formfakv      = 'F'
  integer   :: diag1v        =  1
  integer   :: diag2v        =  1
  integer   :: diag3v        =  0
  character :: atomsphv      = 'F'
  character :: matdipolev    = 'F'
  character :: nonscfv       = 'F'
  character :: dftd3         = 'F'
  character :: nonscfforcesv = 'F'
  character :: dftd3v        = 'F'
  character :: rhogridv      = 'F'
  character :: fragmentv     = 'F'
  character :: solventv      = 'F'
  integer   :: maxscfv       = -1
  real*8    :: scftolv       = -1.0
  character :: moldenv       = 'F'
  character :: nbov          = 'F'
  character :: symmetryv     = 'F'
  character :: mpi_iov       = 'F'
  character :: spnorbv       = 'F'
  character :: fixmv         = 'F'
  character :: pcmv          = 'F'
  character :: wffrmv        = 'F'
  character :: libxcv        = 'F'
  character :: dmatv         = 'F'
  character :: mixingv       = 'P'
  integer   :: nwfoutv       = -1
  character :: spnpolv       = 'F'
  character :: efpv          = 'F'
  character :: populationv   = 'F'
! character :: scaledlbfgsv  = 'F'
  integer   :: meshsettingv  = -1
  character :: unihamv       = 'F'
  character :: frozendenv    = 'F'
  character :: fod_loopv     = 'F'
  character*5 :: fod_opt1v     = 'F'
  character :: fod_opt2v     = 'F'
  integer   :: fod_opt3v     = -1
  character :: symmmodulev = 'F'
  namelist/input_data/atomsphv,diag1v,diag2v,diag3v,dosoccuv,excitedv,formfakv,jntdosv,matdipolev,nonscfv,nonscfforcesv,wfgridv, &
                      dftd3v,calctypev,fragmentv,basisv,rhogridv,solventv,maxscfv,scftolv,moldenv,nbov,symmetryv,mpi_iov,spnorbv, &
                      fixmv,pcmv,wffrmv,dmatv,mixingv,nwfoutv,libxcv,spnpolv,efpv,populationv,meshsettingv,unihamv,&
                      frozendenv, &
                      fod_loopv,fod_opt1v,fod_opt2v,fod_opt3v, symmmodulev

  if (simplified_input) return

  INQUIRE(FILE='NRLMOL_INPUT.DAT',EXIST=EXIST)
  IF (.NOT.EXIST) THEN
    PRINT '(2A)','CHECK_INPUTS: FILE NRLMOL_INPUT.DAT DOES NOT EXIST ',&
    'WILL CREATE A DEFAULT ONE'
    OPEN(68,FILE='NRLMOL_INPUT.DAT')
    write(68,'(A)')'# Put Y,N or number next to the equal sign to determine execution'
    write(68,'(A)')'# Don''t forget the quotation marks for the letters'
    write(68,'(A)')'# All variables in this list end with v'
    write(68,'(A)')
    write(68,'(A)')'&input_data'
    write(68,'(A)')'ATOMSPHV      = ''N'''
    write(68,'(A)')'BASISV        = ''DEFAULT'' ! Specify basis for calculation(basis.txt)'
    IF (OPTTYP.EQ.'CARTESIAN') THEN
      write(68,'(A)')'CALCTYPEV     = ''LBFGS'''
    ELSE IF (OPTTYP.EQ.'REDUNDANT') THEN
      write(68,'(A)')'CALCTYPEV     = ''REDUNDANT'''
    END IF
    write(68,'(A)')'DFTD3V        = ''N'' ! Set to Y to do include Grimmes DFT-D3 dispersion'
    write(68,'(A)')'DIAG1V        =  1  ! diagonalization to use on regular arrays (diagge.f90)'
    write(68,'(A)')'DIAG2V        =  1  ! diagonalization to use on packed arrays (diag_dspgv.f90)'
    write(68,'(A)')'DIAG3V        =  0  ! diagonalization to use on parallel (sdiagge_n.f90)'
    write(68,'(A)')'DMATV         = ''N'' ! Create/use/mix density matrix'
    write(68,'(A)')'DOSOCCUV      = ''N'' ! Controls wether to calculate density of states'
    write(68,'(A)')'FIXMV         = ''N'' ! Fix spin moment'
    write(68,'(A)')'FOD_LOOPV     = ''N'' ! Inernal FOD loop for frozen density optimization'
    write(68,'(A)')'FOD_OPT1V     = ''CG'' ! FOD_OPT: algorithm (LBFGS/CG)'
    write(68,'(A)')'FOD_OPT2V     = ''N'' ! FOD_OPT: scaling of r and F'
    write(68,'(A)')'FOD_OPT3V     =  0  ! FOD_OPT (0)no constraint (1)fix1s (2)fullForce (3)shellOPT (4)freeFOD'
    write(68,'(A)')'FRAGMENTV     = ''N'' ! Process CLUSTER in fragments'
!   write(68,'(A)')'FROZENDENV    = ''N'' ! Frozen density mode (SIC only)'
    write(68,'(A)')'JNTDOSV       = ''N'' ! This calculates joint density of states (DFA only)'
    write(68,'(A)')'MAXSCFV       = 100 ! Maximum SCF iterations'
    write(68,'(A)')'MIXINGV       = ''H'' ! (H)amiltonian (P)otential (D)ensity matrix mixing'
    write(68,'(A)')'NONSCFV       = ''N'' ! Set to Y to do a non SCF calculation' 
    write(68,'(A)')'NONSCFFORCESV = ''N'' ! Set to Y to calculate forces in a non SCF calculation' 
    write(68,'(A)')'NWFOUTV       = 10  ! Write WFOUT file for every N-th iteration'    
    write(68,'(A)')'POPULATIONV   = ''N'' ! Population analysis'
    write(68,'(A)')'RHOGRIDV      = ''N'' ! Set to Y to execute RHOGRID'
!   write(68,'(A)')'SCALEDLBFGSV  = ''Y'' ! Set to Y to scaled LBFGS (SIC only)'
    write(68,'(A)')'SCFTOLV       = 1.0D-6 ! SCF tolerance'
    write(68,'(A)')'SPNPOLV       = ''N'' ! Run spin polarized calculation from CLUSTER'
    write(68,'(A)')'MESHSETTINGV  =  0  ! Mesh recommended for (0)LDA/PBE, (1)SCAN, (2)rSCAN'
    write(68,'(A)')'SYMMETRYV     = ''N'' ! Set to Y to detect symmetry'
    write(68,'(A)')'SYMMMODULEV   = ''N'' ! (T) Use symmtery and approx Ham. (F) use Jacobi rotation (SIC only)'
    write(68,'(A)')'UNIHAMV       = ''N'' ! Set to Y to use unified Hamiltonian formalism in SCF-SIC (SIC only)'
    write(68,'(A)')'WFGRIDV       = ''N'' ! set to Y to write orbitals in cube format (DFA only)'
    write(68,'(A)')'WFFRMV        = ''N'' ! set to Y to write Fermi orbitals in cube format (SIC only)'
    write(68,'(A)')'&end'
    close(68)
! These variables are defined in the global_inputs modules
! Here they are being initialized since no NRLMOL_INPUT.DAT file was found

    idiag1=diag1v
    idiag2=diag2v
    idiag3=diag3v
    dosjnt1=.FALSE.
    wfgrid1=.FALSE.
    dosoccu1=.FALSE.
    formfak1=.FALSE.
    atomsph1=.FALSE.
    matdipole1=.FALSE.
    excited1=.FALSE.
    nonscf1=.FALSE.
    dftd31=.FALSE.
    fragment1=.FALSE.
    rhogrid1=.FALSE.
    calc_basis=1
    solvent1=.FALSE.
    maxscf=100
    scftol=1.0d-6
    molden1=.TRUE.
    cmolden=.FALSE.
    nbo1=.TRUE.
    cnbo=.FALSE.
    symmetry1=.TRUE.
    csymmetry=.FALSE.
    mpi_io1=.FALSE.
    spnorb1=.FALSE.
    fixm1=.FALSE.
    pcm1=.FALSE.
    wffrm1=.FALSE.
    libxc1=.FALSE.
    dmat1=.FALSE.
    mixing1=0
    nwfout=10
    spnpol1=.FALSE.
    scanmesh1=0
    efp1=.FALSE.
    population1=.FALSE.
!   SCALEDLBFGS1=.FALSE.
    symmetrymodule1=.FALSE.
    UNIHAM1=.FALSE.
    FROZENDENSITY1=.FALSE.
  ELSE
    inquire(unit=68,opened=lopen)
    if(lopen) write(6,'(A)')'CHECK_INPUTS:Unit 68 is already open!!'
5   continue
    OPEN(68,FILE='NRLMOL_INPUT.DAT')
    read(68,input_data,ERR = 10,IOSTAT=ios)
    !read(68,input_data,IOSTAT=ios)
    !read(68,input_data)
    if(ios .eq. 0) goto 20
10  continue
    print *,"#################### WARNING ####################"
    print *," THE FOLLOWING NRLMOL_INPUT.DAT VARIABLES IS NOT "
    print *," AVAILABE IN THIS SOFTWARE VERSION AND IGNORED.  "
    backspace(68)
    read(68,'(A)') line
    print*,"Undefiend line: "//trim(line)
    print *,"#################################################"
!Create a temp NRLMOL_INPUT file to process remaining inputs after the problematic line
!This will update NRLMOL_INPUT.DAT with the problematic option commented out
    call system('cat NRLMOL_INPUT.DAT > NRLMOL_INPUT.tmp')
    rewind(68)
    open(69,FILE='NRLMOL_INPUT.tmp')
    do 
      read(69,'(A)', END=15) linebuff
      if(linebuff .eq. line) then
        write(68,'(A)') "!"//trim(linebuff)
      else
        write(68,'(A)') linebuff
      end if
    end do
15  continue
    close(69,status='DELETE')
    close(68)
    goto 5

20  continue
    CLOSE(68)
    idiag1=diag1v
    idiag2=diag2v
    idiag3=diag3v
! Check for DOSJNT
    if(jntdosv=='Y'.or.jntdosv=='y')then
      dosjnt1=.TRUE.
    elseif(jntdosv=='N'.or.jntdosv=='n')then
      dosjnt1=.FALSE.
    else
      if(cdosjnt) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for JNTDOSV'
        write(6,'(A)')'Assuming N'
        dosjnt1=.FALSE.
        cdosjnt=.FALSE.
      endif
    endif
! Check for WFGRID
    if(wfgridv=='Y'.or.wfgridv=='y')then
      wfgrid1=.TRUE.
    elseif(wfgridv=='N'.or.wfgridv=='n')then
      wfgrid1=.FALSE.
    else
      if(cwfgrid) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for WFGRIDV'
        write(6,'(A)')'Assuming N'
        wfgrid1=.FALSE.
        cwfgrid=.FALSE.
      endif
    endif
! Check for DOSOCCU
    if(dosoccuv=='Y'.or.dosoccuv=='y')then
      dosoccu1=.TRUE.
    elseif(dosoccuv=='N'.or.dosoccuv=='n')then
      dosoccu1=.FALSE.
    else 
      if(cdosoccu) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for DOSOCCUV'
        write(6,'(A)')'Assuming N'
        dosoccu1=.FALSE.
        cdosoccu=.FALSE.
      endif
    endif
! Check for FORMFAK
    if(formfakv=='Y'.or.formfakv=='y')then
      formfak1=.TRUE.
    elseif(formfakv=='N'.or.formfakv=='n')then
      formfak1=.FALSE.
    else
      if(cformfak) then 
        !write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for FORMFAKV'
        !write(6,'(A)')'Assuming N'
        formfak1=.FALSE.
        cformfak=.FALSE.
      endif
    endif
! Check for ATOMSPH
    if(atomsphv=='Y'.or.atomsphv=='y')then
      atomsph1=.TRUE.
    elseif(atomsphv=='N'.or.atomsphv=='n')then
      atomsph1=.FALSE.
    else
      if(catomsph) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for ATOMSPHV'
        write(6,'(A)')'Assuming N'
        atomsph1=.FALSE.
        catomsph=.FALSE.
      endif
    endif
! Check for MATDIPOLE
    if(matdipolev=='Y'.or.matdipolev=='y')then
      matdipole1=.TRUE.
    elseif(matdipolev=='N'.or.matdipolev=='n')then
      matdipole1=.FALSE.
    else 
      if(cmatdipole) then
        !write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for MATDIPOLEV'
        !write(6,'(A)')'Assuming N'
        matdipole1=.FALSE.
        cmatdipole=.FALSE.
      endif
    endif
! Check for excited state run
    if(excitedv=='Y'.or.excitedv=='y')then
      excited1=.TRUE.
    elseif(excitedv=='N'.or.excitedv=='n')then
      excited1=.FALSE.
    else
      if(cexcited)then 
        !write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for EXCITEDV'
        !write(6,'(A)')'Assuming N'
        excited1=.FALSE.
      endif
    endif
! Check for NONSCF
    if(nonscfv=='Y'.or.nonscfv=='y')then
      nonscf1=.TRUE.
    elseif(nonscfv=='N'.or.nonscfv=='n')then
      nonscf1=.FALSE.
    else
      if(cnonscf)then 
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for NONSCFV'
        write(6,'(A)')'Assuming N'
        nonscf1=.FALSE.
        cnonscf=.FALSE.
      endif
    endif
! Check for non scf forces
    if(nonscfforcesv=='Y'.or.nonscfforcesv=='y')then
      nonscfforces1=.TRUE.
    elseif(nonscfforcesv=='N'.or.nonscfforcesv=='n')then
      nonscfforces1=.FALSE.
    else
      if(cnonscfforces) then 
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for NONSCFFORCESV'
        write(6,'(A)')'Assuming N'
        nonscfforces1=.FALSE.
        cnonscfforces=.FALSE.
      endif
    endif
! Check for dftd3
    if(dftd3v=='Y'.or.dftd3v=='y')then
      dftd31=.TRUE.
    elseif(dftd3v=='N'.or.dftd3v=='n')then
      dftd31=.FALSE.
    else
      if(cdftd3) then 
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for DFTD3V'
        write(6,'(A)')'Assuming N'
        dftd31=.FALSE.
        cdftd3=.FALSE.
      endif
    endif
! Check for calc_type
    if(calctypev=='LBFGS'.or.calctypev=='lbfgs')then
      calctype1 = 1
    elseif(calctypev=='SCF-ONLY'.or.calctypev=='scf-only')then
      calctype1 = 2
    elseif(calctypev=='VERLET'.or.calctypev=='verlet')then
      calctype1 = 3
    elseif(calctypev=='VIBRATIONAL'.or.calctypev=='vibrational')then
      calctype1 = 4
    elseif(calctypev=='REDUNDANT'.or.calctypev=='redundant')then
      calctype1 = 5
    else 
      if(ccalctype) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for CALCTYPEV'
        write(6,'(A)')'Assuming LBFGS'
        calctype1 = 1
        ccalctype = .FALSE.
      endif
    endif
! Check for fragments
    if(fragmentv=='Y'.or.fragmentv=='y')then
      fragment1=.TRUE.
    elseif(fragmentv=='N'.or.fragmentv=='n')then
      fragment1=.FALSE.
    else
      if(cfragment) then 
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for FRAGMENTV'
        write(6,'(A)')'Assuming N'
        fragment1=.FALSE.
        cfragment=.FALSE.
      endif
    endif
! Check for basis set (this has to be the last input to be processed)
! Newer additions to NRLMOL_INPUT.DAT must be place above this
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
      !grab basis set file name
      calc_basis=100
      basis_filename=basisv
      !if(cbasis)then
      !  write(6,'(A)')'CHECK_INPUTS:WRONG value for BASISV'
      !  write(6,'(A)')'Assuming DEFAULT'
      !  calc_basis=1
      !  cbasis=.FALSE.
      !endif
    endif
! Check for RHOGRID
    if(rhogridv=='Y'.or.rhogridv=='y')then
      rhogrid1=.TRUE.
    elseif(rhogridv=='N'.or.rhogridv=='n')then
      rhogrid1=.FALSE.
    else
      if(crhogrid) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for RHOGRIDV'
        write(6,'(A)')'Assuming N'
        rhogrid1=.FALSE.
        crhogrid=.FALSE.
      endif
    endif
! Check for SOLVENT
    if(solventv=='Y'.or.solventv=='y')then
      solvent1=.TRUE.
    elseif(solventv=='N'.or.solventv=='n')then
      solvent1=.FALSE.
    else
      if(csolvent) then
        !write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for SOLVENTV'
        !write(6,'(A)')'Assuming N'
        solvent1=.FALSE.
        csolvent=.FALSE.
      endif
    endif
! Check for MAXSCF
    if(maxscfv==-1)then
      if(cmaxscf) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for MAXSCFV'
        write(6,'(A)')'Assuming 100'
        maxscf=100
        cmaxscf=.FALSE.
      endif
    else
        maxscf=maxscfv
    endif
! Check for SCFTOL
    if(scftolv==-1.0)then
      if(cscftol) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for SCFTOLV'
        write(6,'(A)')'Assuming 1.0D-6'
        scftol=1.0D-6
        cscftol=.FALSE.
      endif
    else
      scftol=scftolv
    endif
! Check for MOLDENV
    if(moldenv=='Y'.or.moldenv=='y')then
      molden1=.TRUE.
    elseif(moldenv=='N'.or.moldenv=='n')then
      molden1=.FALSE.
    else
      if(cmolden) then
        !write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for MOLDENV'
        !write(6,'(A)')'Assuming N'
        molden1=.FALSE.
        cmolden=.FALSE.
      endif
    endif
! Check for NBOV
    if(nbov=='Y'.or.nbov=='y')then
      nbo1=.TRUE.
    elseif(nbov=='N'.or.nbov=='n')then
      nbo1=.FALSE.
    else
      if(cnbo) then
        !write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for NBOV'
        !write(6,'(A)')'Assuming N'
        nbo1=.FALSE.
        cnbo=.FALSE.
      endif
    endif
! Check for SYMMETRYV
    if(symmetryv=='Y'.or.symmetryv=='y')then
      symmetry1=.FALSE.
    elseif(symmetryv=='N'.or.symmetryv=='n')then
      symmetry1=.TRUE.
    else
      if(csymmetry) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for SYMMETRYV'
        write(6,'(A)')'Assuming Y'
        symmetry1=.TRUE.
        csymmetry=.FALSE.
      endif
    endif
! Check for SYMMETRYMODULEV
    if(symmmodulev=='Y'.or.symmmodulev=='y')then
      symmetrymodule1=.TRUE.
    elseif(symmmodulev=='N'.or.symmmodulev=='n')then
      symmetrymodule1=.FALSE.
    else
      symmetrymodule1=.FALSE.
    endif
! Check for MPI_IOV
    if(mpi_iov=='Y'.or.mpi_iov=='y')then
      mpi_io1=.TRUE.
    elseif(mpi_iov=='N'.or.mpi_iov=='n')then
      mpi_io1=.FALSE.
    else
      if(cmpi_io) then
        !write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for MPI_IOV'
        !write(6,'(A)')'Assuming N'
        mpi_io1=.FALSE.
        cmpi_io=.FALSE.
      endif
    endif
! Check for SPNORBV
    if(spnorbv=='Y'.or.spnorbv=='y')then
      spnorb1=.TRUE.
    elseif(spnorbv=='N'.or.spnorbv=='n')then
      spnorb1=.FALSE.
    else
      if(cspnorb) then
        !write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for SPNORBV'
        !write(6,'(A)')'Assuming N'
        spnorb1=.FALSE.
        cspnorb=.FALSE.
      endif
    endif
! Check for FIXMV
    if(fixmv=='Y'.or.fixmv=='y')then
      fixm1=.TRUE.
    elseif(fixmv=='N'.or.fixmv=='n')then
      fixm1=.FALSE.
    else
      if(cfixm) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for FIXMV'
        write(6,'(A)')'Assuming N'
        fixm1=.FALSE.
        cfixm=.FALSE.
      endif
    endif
! Check for WFFRMV
    if(wffrmv=='Y'.or.wffrmv=='y')then
      wffrm1=.TRUE.
    elseif(wffrmv=='N'.or.wffrmv=='n')then
      wffrm1=.FALSE.
    else 
      wffrm1=.FALSE.
    endif
! YY. Check for LIBXC  
!    if(libxcv=='Y'.or.libxcv=='y')then
!      libxc1=.TRUE.
!    elseif(libxcv=='N'.or.libxcv=='n')then
!      libxc1=.FALSE.
!    else
!     if(clibxc) then
!      write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for LIBXCV'
!      write(6,'(A)')'Assuming N'
!      libxc1=.FALSE.
!      clibxc=.FALSE.
!     endif
!    endif
! Check for DMATV
    if(dmatv=='Y'.or.dmatv=='y')then
      dmat1=.TRUE.
    elseif(dmatv=='N'.or.dmatv=='n')then
      dmat1=.FALSE.
    else
      if(cdmat) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for DMATV'
        write(6,'(A)')'Assuming N'
        dmat1=.FALSE.
        cdmat=.FALSE.
      endif
    endif
! Check for MIXINGV
    if(mixingv=='P'.or.mixingv=='p')then
      mixing1=0
    elseif(mixingv=='D'.or.mixingv=='d')then
      mixing1=1
    elseif(mixingv=='H'.or.mixingv=='h')then
      mixing1=2
    else
      if(cmixing) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for MIXINGV'
        write(6,'(A)')'Assuming P'
        mixing1=0
        dmat1=.FALSE.
        cdmat=.FALSE.
      endif
    endif
! If Meta-GGA or Libxc is used, automatically turn on Ham. mixing.
    if(mixing1.ne.2)then
     if(IGGA(1).eq.3 .or. IGGA(2).eq.3)then
      !write(6,'(A)')'>>>> Meta-GGA detected. Hamiltonian mixing will be used <<<<'
      mixing1=2
     elseif(libxcv=='Y'.or.libxcv=='y')then
      !write(6,'(A)')'>>>> Libxc detected. Hamiltonian mixing will be used <<<<'
      mixing1=2
     endif
    endif
! Check for NWFOUTV
    if(nwfoutv < 0)then
        !write(6,'(A)')'Assuming 10'
        nwfout=10
    else
        nwfout=nwfoutv
    endif
! Check for spnpolv
    if(spnpolv=='Y'.or.spnpolv=='y')then
      spnpol1=.TRUE.
    elseif(spnpolv=='N'.or.spnpolv=='n')then
      spnpol1=.FALSE.
    else
      spnpol1=.FALSE.
    endif
! Check for MESHSETTINGV (SCANMESHV)                                               
    if(meshsettingv .ge. 0 .and. meshsettingv .le. 2)then                   
      scanmesh1=meshsettingv                                              
    else                                                                 
      scanmesh1=0
    endif
    if(scanmesh1 .ne. 1) then
      !if(IGGA(1).eq.3 .or. IGGA(2).eq.3)then mgga->scan
      if(IDFTYP(1).eq.11 .or. IDFTYP(2).eq.11) then
        !write(6,'(A)')'>>>> SCAN functional detected. SCAN mesh will be used <<<<'
        scanmesh1=1
      endif
    endif
! Check for EFPV
    if(efpv=='Y'.or.efpv=='y')then
      efp1=.TRUE.
    elseif(efpv=='N'.or.efpv=='n')then
      efp1=.FALSE.
    else
      if(cefp) then
        write(6,'(A)')'CHECK_INPUTS:Wrong or missing value for EFPV'
        write(6,'(A)')'Assuming N'
        efp1=.FALSE.
        cefp=.FALSE.
      endif
    endif
! Check for POPULATION
    if(populationv=='Y'.or.populationv=='y')then
      population1=.TRUE.
    elseif(populationv=='N'.or.populationv=='n')then
      population1=.FALSE.
    else
      population1=.FALSE.
    endif
! Check for SCALEDLBFGS
!   if(scaledlbfgsv=='Y'.or.scaledlbfgsv=='y')then
!     SCALEDLBFGS1=.TRUE.
!   else
!     SCALEDLBFGS1=.FALSE.
!   end if
! Check for UNIHAM
    if(unihamv=='Y'.or.unihamv=='y')then
      UNIHAM1=.TRUE.
    else
      UNIHAM1=.FALSE.
    end if
! Check for FROZENDENSITY
    if(frozendenv=='Y'.or.frozendenv=='y')then
      FROZENDENSITY1=.TRUE.
    else
      FROZENDENSITY1=.FALSE.
    end if
    if(FOD_LOOPV=='Y'.or.FOD_LOOPV=='y')then
      FOD_LOOP1=.TRUE.
    else
      FOD_LOOP1=.FALSE.
    end if
    if(FOD_OPT1V=='LBFGS'.or.FOD_OPT1V=='lbfgs')then
      FOD_OPT1=.TRUE.
    else
      FOD_OPT1=.FALSE.
    end if
    if(FOD_OPT2V=='Y'.or.FOD_OPT2V=='y')then
      FOD_OPT2=.TRUE.
    else
      FOD_OPT2=.FALSE.
    end if
    !FOD_OPT3V default is 0
    FOD_OPT3=FOD_OPT3V
    if(FOD_OPT3V.LT.0 .OR. FOD_OPT3V.GT.4)then
      FOD_OPT3=0
    end if
   ENDIF
!Print out name_list for testing and debugging this subroutine
!  open(5,file='NRLMOL_OUTPUT.DAT')
!  write(5,input_data)
!  close(5)
!  call stopit
  RETURN
end subroutine check_inputs
