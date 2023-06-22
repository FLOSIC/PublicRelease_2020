! UTEP Electronic Structure Lab (2020)
module debug1

  USE MPIDAT1,only : IRANK
! DEBUG:  determines if debug info will be printed
  LOGICAL :: DEBUG
  CONTAINS
!  SUBROUTINE TRACER(MSG,VALUE1)
  SUBROUTINE TRACER(MSG,VALUE1,VALUE2)
  CHARACTER(LEN=*) :: MSG
  INTEGER,OPTIONAL :: VALUE1
  REAL*8,OPTIONAL  :: VALUE2
  IF(PRESENT(VALUE2)) THEN
    WRITE(6+IRANK,*) 'TRACER:',MSG,VALUE1,VALUE2
  ELSE
    IF(PRESENT(VALUE1)) THEN
      WRITE(6+IRANK,*) 'TRACER:',MSG,VALUE1
    ELSE
      WRITE(6+IRANK,*) 'TRACER:',MSG
    ENDIF
  ENDIF
  CALL FLUSH(6+IRANK)
  END SUBROUTINE TRACER

  subroutine tick(t)
  integer, intent(OUT) :: t

  call system_clock(t)
  end subroutine tick

! returns time in seconds from now to time described by t
  real function tock(t)
  integer, intent(in) :: t
  integer :: now, clock_rate

  call system_clock(now,clock_rate)

  tock = real(now - t)/real(clock_rate)
  end function tock

  pure real(8) function mat_trace(A)
    implicit none
    real(8),intent(in) :: A(:,:)
    integer :: i,n

    n=size(A,1)
    mat_trace=0.0d0
    do i=1,n
      mat_trace=mat_trace+A(i,i)
    end do
    return
  end function

  pure real(8) function mat_trace_cube(A,idx)
    implicit none
    real(8),intent(in) :: A(:,:,:)
    integer,intent(in) :: idx
    integer :: i,n

    n=size(A,1)
    mat_trace_cube=0.0d0
    do i=1,n
      mat_trace_cube=mat_trace_cube+A(i,i,idx)
    end do
    return
  end function

end module debug1

module common1

! ALPCOR: largest exponent of NLCC core density
! PSPANG: lmax of different radial zones used for nonlocal PSPs
! PSPINF: general PSP info 
! PSPNLO: tabulated nonlocal PSPs
! BHSPSP: data for local potential of type BHS
! TABPSP: data for local potential of type TAB (including NLCC)

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! ALPCOR:
  REAL*8  :: ALPCOR(MAX_FUSET)
! PSPANG:
  REAL*8,DIMENSION(4)  :: PSRZONE=(/0.2D0, 0.4D0, 0.8D0, 1.6D0/)
  INTEGER,DIMENSION(5) :: LMXPSRZ=(/7, 9, 11, 13, 15/)
! PSPINF: 
  CHARACTER*7 :: PSPSYM(MAX_FUSET)
  INTEGER :: ISITPSP,ISNLCC
! PSPNLO:
  REAL*8  :: RPSNLO(MXRPSP,MAX_FUSET),WPSNLO(MXRPSP,MAX_FUSET),&
             VPSNLO(MXLPSP+1,MXRPSP,MAX_FUSET)
  INTEGER :: LMAXNLO(MAX_FUSET),NRPSP(MAX_FUSET)
! BHSPSP:
  REAL*8  :: BHSALP(2,MAX_FUSET),BHSCOF(2,MAX_FUSET)
! TABPSP:
  REAL*8  :: RRADTAB(MXPTAB,MAX_FUSET),VLRTAB(2,MXPTAB,MAX_FUSET),&
             RHOCOR(3,MXPTAB,MAX_FUSET)
  INTEGER :: NRADTAB(MAX_FUSET),NLCC(MAX_FUSET)

end module common1

module common2

!
! NUCLEI:   information about atom locations and types
! BASET:    basis set / nuclear charge data for atom types
! FERMIONS: number of electrons, wavefunction output file
! DFTYP:    information about density functional
! SPIN:     spin info
! DIPOLE:   dipole moment and external homegeneous electric field
! ENERG:    stores a lot of different energies
! FORCES:   force information
! NUMFRC:   used for calculation Pulay corrections to the forces
! OPTGEOM:  data needed in geometry optimizations
! SYMBNAME: contains the symbols used for the identity members
!JUR
! IHOMOP:   TOTAL NUMBER OF MO IS 2*IHOMOP+1, HOMO STATE IS IN THE MIDDLE
!JUR
!
! define maximum number of symbols
!
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
  integer MXSYMBS 
  PARAMETER (MXSYMBS=MAX_IDENT+30)

!
! NUCLEI:   information about atom locations and types
  REAL*8  :: RIDT(3,MAX_IDENT),RCNT(3,MX_CNT)
!
!CJUR
! ATOMNO:    Atomic numbers.
! ATOMLAB:   Atomic labels.
! ATOMBAS:   Atomic basis flag.
! ATOMMAP:   Pointer for Cartesian -> Z-Matrix order.
! ATOMSPIN:  Atom sping definition, i.e. SUP, SDN and UPO
! CHNET:     Net atomic charge.
! COORD:     Atomic coordinates.
! CONMAT:    Connectivity matrix.
! DE:        Energy gradient of optimization step.
! DR:        Geometry change of optimization step.
! EGMIN:     Minium gap for printing MO orbitals.
! EGMAX:     Maximum gap for printing MO orbitals.
! EPREDICT:  Predicted energy change.
! ESTEPS:    Energy of optimization steps.
! HESSIAN:   Start Hessian for geometry optimization.
! HFOLLOW:   Hessian mode to follow for geometry optimization.
! INDICES:   Indices labeling MO's in Molden output.
! INPTYP :   Flag for input type.
! MATOM:     Atomic mass.
! MAXDR:     Maximum step size for optimization.
! MOLDEN:    Flag for Molden output.
! NA,NB,NC:  Atomic connectivity in Z-Matrix
! NATOM:     Number of atoms.
! NATOMS:    Total number of atoms
! NDUMMY:    Number of dummy atoms.
! NACT:      Number of active optimization variables. 
! NDEG:      Number of degrees of freedom, including dummy atoms. 
! NNEG:      Number of negative eigenvalues in Hessian.
! NPRI:      Number of primitive internal coordinates.
! NVAR:      Number of primitive variables.
! NVARRST:   Number of processed variables.
! NXYZ:      Number of Cartesian (real atom) degrees of freedom.
! NFXYZ:     Number of frozen Cartesian coordinates.
! NCODE:     Flag for angle/dihedral angle substitution
! NBO:       Natural bonding orbital analysis.
! OPTCYC:    Number of optimization cycle.
! OPTFLAG:   Z-Matrix optimization flags for internal coordinates
! OPTTYP:    Geometry optimization coordinate type.    
! OPTTOL:    Tolerance for geometry optimization convergence.
! OPTVAR:    Z-Matrix optimization variable and constant names
! OPTXYZ:    Cartesian freezing string for coordinates
! PRILAB:    Primitive internal coordinate label
! PRIPTR:    Pointer for the primitive entry in a non-redundant
!            delocalized coordinate.
! PRTOPT:    If true, print optimization information.
! PRTPRI:    If true, print primitive coordinate information. 
! PRTMOS:    If true, print primitive Molecular orbitals.
! PRTPMAT:   If true, print density matrix.
! PRTSMAT:   If true, print overlap matrix.
! SSCALED:   If true, internal displacement is scaled.
! SLENGTH:   Step length.
! SPNNET:    Net spin magnetic moment.
! UNITS:     Flag for BOHR or ANGSTROM unit system.
! VALENCE:   Atomic valence matrix
! ZATOM:     Nuclear charges
! CPRIMS:    Connector list of primitive coordinates. 
! NPRIMS:    Number of primitive length and angle connectors
! QPRIMS:    Primitive internal coordinates.
! OPTIMIZED: Logical Flag: If true, optimization converged.
! OPTRST:    If true, restart geometry optimization.
! RTRUST:    Trust radius.
! STEPTYP:   Step type for geometry optimization (DEFAULT).
!
  INTEGER, PARAMETER :: MAXPRI = 960
  INTEGER, PARAMETER :: MAXKEY = 18
  INTEGER, PARAMETER :: MAXVAR = 15
  INTEGER, PARAMETER :: MAXNUM = 15
  INTEGER, PARAMETER :: MAXOPT = 31
  INTEGER, PARAMETER :: INP =    32

  LOGICAL :: NBO,OPTIMIZED,OPTRST,PRTOPT,PRTPRI,PRTMOS,PRTPMAT,PRTSMAT
  LOGICAL :: SSCALED
!
  CHARACTER*2  :: ATOMLAB(MAX_IDENT)
  CHARACTER*3  :: ATOMBAS(MAX_IDENT),ATOMSPIN(MAX_IDENT)
  CHARACTER*3  :: NCODE(MAX_IDENT),OPTXYZ(MAX_IDENT)
  CHARACTER*10 :: PRILAB(MAXPRI),OPTVAR(MAX_IDENT,3)
  CHARACTER*10 :: INPTYP,HESSIAN,MOLDEN,OPTTYP,STEPTYP
  CHARACTER*15 :: UNITS
!
  INTEGER :: HFOLLOW,NATOM,NATOMS,NDEG,NDUMMY,NNEG,NPRI,NACT,NVAR,&
             NXYZ,NFXYZ,OPTCYC,NVARRST
  INTEGER :: NA(MAX_IDENT),NB(MAX_IDENT),NC(MAX_IDENT),NPRIMS(-3:3),&
             ATOMNO(MAX_IDENT),ZATOM(MAX_IDENT),ATOMMAP(MAX_IDENT,4),&
             OPTFLAG(MAX_IDENT,3),CPRIMS(MAXPRI,6),INDICES(4,MAX_OCC),&
             CONMAT(MAX_IDENT,MAX_IDENT),PRIPTR(MAXPRI)
  REAL*8  :: CHNET,EGMIN,EGMAX,EPREDICT,MAXDR,OPTTOL,RTRUST,SLENGTH
  REAL*8  :: SPNNET,ESTEPS(2),MATOM(MAX_IDENT)
  REAL*8  :: VALENCE(MAX_IDENT,MAX_IDENT),COORD(3,MAX_IDENT,5)
  REAL*8  :: QPRIMS(MAXPRI,2),DE(3*MAX_IDENT,2),DR(3*MAX_IDENT,2)
!
  TYPE NRLMOL_KEYWORD
    SEQUENCE
    CHARACTER*32 :: FULL_KEYWORD
    INTEGER      :: NUMBER_OF_OPTIONS
    LOGICAL      :: WRITE_KEYWORD_IN_NEW_INPUT
  END TYPE NRLMOL_KEYWORD
!
  TYPE (NRLMOL_KEYWORD) :: KEYWORDS(MAXKEY)
!
!CJUR
! 
  INTEGER :: IFUIDT(MAX_IDENT),IFUCNT(MX_CNT),NIDENT,NCNT
!
! BASET:    basis set / nuclear charge data for atom types
  REAL*8  :: ZELC(MAX_FUSET),ZNUC(MAX_FUSET),&
             BFCON(MAX_BARE,MAX_CON,LDIM,MAX_FUSET),&
             BFALP(MAX_BARE,MAX_FUSET)
  INTEGER :: N_BARE(MAX_FUSET),N_CON(LDIM,MAX_FUSET),&
             LSYMMAX(MAX_FUSET),N_POS(MAX_FUSET),NFNCT
! FERMIONS: number of electrons, wavefunction output file
  REAL*8  :: E_UP,E_DN
  CHARACTER*40 :: WFFILE
! DFTYP:    information about density functional
  INTEGER :: IGGA(2),IDFTYP(2)
! SPIN:     spin info
  INTEGER :: ISPN,NSPN
! DIPOLE:   dipole moment and external homegeneous electric field
  REAL*8,DIMENSION(3)  :: DIPOLE=(/0.0,0.0,0.0/)
  REAL*8,DIMENSION(3)  :: EFIELD=(/0.0,0.0,0.0/)
! ENERG:    stores a lot of different energies
  REAL*8  :: ETOTAL,ENNUC,ELOCAL,ECOUL,ERGFLD,EKINONL,ENONLO,&
             ERGXL,ERGXN,ERGCL,ERGCN,EDISP,ESOLTOT,ESOLC,ENUCSOL,&
             EPCM,EKINONL2
! FORCES:   force information
  REAL*8  :: FHELLF(3,MAX_IDENT),FNONL(3,MAX_IDENT),&
             FTOT(3,MAX_IDENT)
! NUMFRC:   used for calculation Pulay corrections to the forces
  REAL*8  :: FRC1(3,MAX_IDENT),FRC2(3,MAX_IDENT)
! Force 0.0
  REAL*8  :: over1(max_occ,max_occ),ek(max_occ,max_occ)
  REAL*8  :: dftV(max_occ,max_occ),allH(max_occ,max_occ)
! OPTGEOM:  data needed in geometry optimizations
  REAL*8  :: ATMSYM(MAX_IDENT),RIDSYM(3,MAX_IDENT),GTOL,MOVING(3,MAX_IDENT)
  CHARACTER*80 :: SYMBS(2,MXSYMBS)
  INTEGER :: NSCHAR(2,MXSYMBS),NIDSYM,NSYMBS
! SYMBNAME: contains the symbols used for the identity members
  CHARACTER :: SYMATM(10,MXSYMBS)
  REAL*8 :: GSUM,GMAX

end module common2

module common3

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
!
! GROUP: point group operations
!
  REAL*8 :: RMAT(3,3,MX_GRP)
  INTEGER :: NGRP,MULTAB(MX_GRP,MX_GRP)

end module common3

module common4

! STPOT:   data for spin dependence of start potential
! RHOPFIT: data for fitted atomic potentials and density

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! STPOT:   data for spin dependence of start potential
  INTEGER :: ISPIDT(MAX_IDENT)
! RHOPFIT: data for fitted atomic potentials and density
  REAL*8 :: RPFALP(MAX_FUSET),RPFCMX(MAX_FUSET),RPFCOF(2,MAXLSQF,MAX_FUSET)
  INTEGER :: NRPFIT(MAX_FUSET),LDIVR(MAX_FUSET)

end module common4

module common5
! SCFDAT:   SCF and startup data
! CBLK4:    contains wavefunctions and their occupation numbers
! UNSYM:    contains unsymmetrized wavefunctions
! CBLK7:    used to store the hamiltonian and overlap matrices (GONE: now in
!            module hstor1)
! CBLK30:   hamiltonian submatrix and temporary array needed during H setup
! EFRMI:    Fermi energy for each spin system
! ELEVELS:  contains occupied eigenvalues
! FOR_DIAG: used during matrix diagonalization (GONE: now in module for_diag1)
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! SCFDAT:   SCF and startup data
  LOGICAL :: CONVERGENCE,HAVEHAM,FOD_LOOP
  INTEGER :: ISTSCF,IHIPOL
! CBLK4:    contains wavefunctions and their occupation numbers
  REAL*8  :: PSI_COEF(NDH,MAX_VIRT_PER_SYM,MAX_REP,MXSPN),&
             OCCUPANCY(MAX_VIRT_PER_SYM*MAX_REP*MXSPN)
  INTEGER :: N_OCC(MAX_REP,MXSPN)
! UNSYM:    contains unsymmetrized wavefunctions
  REAL*8  :: PSI(MAXUNSYM,MAX_OCC,2)
  INTEGER :: NWF,NWFS(MXSPN)
! CBLK30:   hamiltonian submatrix and temporary array needed during H setup
  REAL*8  :: HOLD(MAXUNSYM,MAXUNSYM),HTEMP(MTEMP)
! EFRMI:    Fermi energy for each spin system
  REAL*8  :: EFERMI(2)
! ELEVELS:  contains occupied eigenvalues
  REAL*8  :: EVLOCC(MAX_OCC)

end module common5

module common6
! MESH:   mesh points and weights
! LB: Common block MESH now used trough module mesh1
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! SPH:    data for spheres that contains subsets of mesh points
  REAL*8  :: TSPH(4,MX_SPH)
  INTEGER :: LIMSPH(2,MX_SPH),NSPHERES

end module common6

module common7
! FORDEN:   determines if the density is created in COUPOT or DENSOLD
! CUTOFF:   defines cutoff for exponents in COUPOT
! WASTED:   time wasted in UNRAVEL

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! FORDEN:   determines if the density is created in COUPOT or DENSOLD
  INTEGER :: MODDEN
! CUTOFF:   defines cutoff for exponents in COUPOT
  REAL*8  :: GAUSS_CUT(MAX_IDENT)
! WASTED:   time wasted in UNRAVEL
  REAL*8  :: T1UNRV=0.0D0
  REAL*8  :: T2UNRV=0.0D0

end module common7

module common8
!
! symmetry stuff 
!
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! BASE_REP
  REAL*8  :: S_REP(MX_GRP),P_REP(3,3,MX_GRP),D_REP(6,6,MX_GRP)
! CBLK11
  REAL*8  :: REP(5,5,MX_GRP,MAX_REP)
  INTEGER :: N_REP,NDMREP(MAX_REP)
! CBKL12
  REAL*8  :: U_MAT(LOCMAX,ISMAX,MAXSYMSALC,3,2)
  INTEGER :: N_SALC(MAXSYMSALC,3,MAX_IDENT)
! CBLK18
  INTEGER :: IGEN(MX_GRP)
! INFOR
  REAL*8  :: RDENT(3,MAX_IDENT)
  INTEGER :: NUMSITES(MAX_IDENT),IGGEN(MX_GRP,MAX_IDENT),N_IDNT
! INTS_INDEX
  INTEGER :: INDBEG(MAX_IDENT,MAX_REP),NS_TOT(MAX_REP)
! REDREP
  INTEGER  :: LDMREP(MAX_REP)

end module common8

module dosjnt_mod

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
 REAL*8,ALLOCATABLE :: H(:,:,:),PSIG(:,:),SPTOT(:),SPDIP(:),SOS_FREQ(:), &
                       TEMP_RCV(:),TEMP_RCV2(:,:,:)
 REAL*8  :: P(MPBLOCK,3),Q(MPBLOCK,3),V(MPBLOCK)
 REAL*8,ALLOCATABLE :: RVECA(:,:),PTS(:,:),GRAD(:,:,:,:,:)
 LOGICAL,ALLOCATABLE :: ICOUNT(:,:)
 REAL*8  :: ESTEP,EALP,SOS_POL,VFAC
 REAL*8  :: ENJD,EXJD,TEMP,FCGRP,CHARGE
 REAL*8,PARAMETER :: HA2EV=27.2116D0
 INTEGER :: NSPEC
 INTEGER,PARAMETER,DIMENSION(3) :: ISIZE=(/1,3,6/)

end module dosjnt_mod

module fragment

 INTEGER :: NFRAGMENT
 INTEGER,ALLOCATABLE :: FRAGMENTS(:)

end module fragment

module solvent

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! REAL*8 :: POTSOL(MAX_PTS)
! REAL*8 :: PCMSOL(MAX_PTS)
! Currently not used and changed into allocatable arrays. Allocate those as needed.
 REAL*8,ALLOCATABLE :: POTSOL(:)
 REAL*8,ALLOCATABLE :: PCMSOL(:)

end module solvent

module unravl

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
 LOGICAL :: OCCL(MAX_OCC),OCCV(MAX_OCC)
 
end module unravl


