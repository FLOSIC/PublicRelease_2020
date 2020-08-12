C> @file electronic_geometry.f
C> @note YY. This subroutine is a draft. Not meant for production.
C> FOD force calculation and optimization 
C> SIC
C> create an empty file, FODLBF, to use LBFGS
C> @param[out] energy
C> @param xvec : FOD positions
C> @param gvec : FOD forces
C> @param msite(2) : This array holds numbers of FODs for spin up/down. 

C Modifications KT
C August 2nd, 2019
C Constrained FOD optimization
C Use force on one FOD per structural element (no projections) on all FODs of that element -> automatically treat rotations and stretching!
C Use rotation matrix
C symmetrize FODs after optimization
C SEPARATE UP AND DN CHANNEL !!
C Combine all approached: Unconstrained, fix1s and fully constrained (fullForce). Make all available by choice of the user. Using CG or LBFGS and using scaled or unscaled forces.
C
C October 11th, 2019
C Use internal coordinates -> easier to optimize, no need to carry all
C forces/positions
C October 11th, 2019
C Start introducing contraints for bond FODs 
C
C October 21th, 2019
C force component along bond axis for bond FODs (in heteronuclear
C       bonding situations)!
C -> move them along the axis as well as the breathing mode
C
C October 28th, 2019
C Implemented scaling of coordinates/forces (by Kushantha
C Withanage, see scaled_lbfgs.f as well)
C
C October 30th, 2019
C Spin-restricted case. If UP and DN are very close after optimization
C -> put them at the same position
C Reset LBFGS: If LBFGS is stuck, reset it.
C
C Think about aborting a calculation when it is converged (energy and
C forces) ---- TBD
C
C November 19th, 2019
C Introduce two internal coordinates for structural motifs with four or
C more points -> include a proper rotation of the motif.
C
C December 10th, 2019
C Only use one internal coordinate for tetrahedra. Only breathing mode
C This overrules the last change
C
C December 17th, 2019
C Internal coordinates for all structural motif.
C TBD: average forces for any structural motif. -> Justlike
C for the tetrahedra. More consistent!
C
C January 14th, 2020
C Use proper internal coordinates. If the internal coordinates and
C forces of a given structural motif are verysimilar to another
C structural motif -> use only one set of coordinates and forces TBD
C (average them). 
C What we do so far: Set certain forces equal to similar forces which
C have already been found. NOW: average them!
C
C January 15th, 2020 
C Idea Sebastian: 'fixCore' mode. fix 1s and breath tetrehedra. leave rest free
C Maybe fix all tetrahedra. Leave the rest of the valence free 
C (e.g. Kr 3s3p3d)
C TBD !!!!!!!
C
C January 20th, 2020
C identify X-H bond FODs as bonds. Special case for H-bonds.
C     January 22nd: Almost there
C     January 23rd: Closer
C     January 24th: Should be done
C
C January 31th, 2020
C     TO BE DONE :
C   - Include Kushantha Withanages scaled scheme in this subroutine (no need
C     for scaledlbfgs.f).  -- DONE
C   - Read atomic coordinates from file -> avoid the need for common
C     blocks               -- DONE 
C   - Reduce the number of internal coordinates to be used in the
C     optimizer to the bare minimum. Just have a look atthe internal
C     coordinates right before the optimizer is called. If
C     coordinates/forces are the same -> only use one of them. -- DONE
C February 7th, 2020
C   take correct orientations of the internal forces/gradients. DONE
C
C February 13th, 2020
C    - Introduce ShellOPT:
C      optimize different structural object independently -> should give
C      smoother optimization behaviour, and no scaling is needed.
C      Thus, only feed one internal coordinate to the optimizer and
C      optimize it. Use the one with largest force!
C      If the forces on a given object drop below 10^-4 -> take next
C      internal coordinate
C      Continue until all are converged.
C        Change file names to FDIAG1, FDIAG2,same for FSEARCH and FSTEP
C  February 17th: If all objects have forces below 10**-4: Optimize them
C                 all together. With smaller LBFGS_stepsize. Serves as
C                 final convergence. Let's see whether it works
C  February 20th: Introduce shellOPT as an option.
C 
C February 21th, 2020
C    ensure that LBFGS always get the right number of coordinates. If
C    the number of coordinates changes -> reset LBFGS
C
C    Idea: get array called free_FOD. unconstrain the FODs who's index
C    is in this array. Partially constrain, partially do not. Could be
C    read from input file, and could be post-processed like the fix1s
C    method.   

         subroutine fod_opt(energy,fod_converge)
         use global_inputs,only : fod_opt1,fod_opt2,fod_opt3
         use xmol,only : AU2ANG,NUM_ATMS,XMOL_LIST,GET_LETTER
         implicit none
         logical :: exist, exist1
         ! Read this from input file. For now, call it FOD_OPT
         logical :: do_lbfgs, scaled, reset         ! decide what kind of optimization to perform. lbfgs_cg decides which optimizer to use   (former LFODLBF)
         logical :: fod_converge
         logical :: first
                                                    ! scaled decides whether to use the Kushanthas scaling of the forces/positions
                                                    ! reset determines whether LBFGS is stuck and needs to be reset
         integer :: constraint                      ! constraint decides what constraints shall be used (0 .. unconstrained, 1 .. fix1s, 2 .. fullForce constrained)
         character(len=20) :: junk1, junk2          ! to read the file correctly, and to write nicer screen output
         !
         ! Read atomic coordinates and species from file, e.g. SYMBOL.
         ! Do not use common blocks
         ! Use a subroutine to do this, it will be neater
         ! 
         integer             :: n_atoms             ! number of atoms
         real(8),allocatable :: r_atoms(:,:)        ! atomic coodinates
         character(len=2),allocatable :: species(:) ! atomic species


         real(8), intent(in) :: energy
         integer, allocatable :: n_motifs(:)        ! number of structural motifs per atom -> determines number of degrees of freedom to optimize!
C KT
         integer :: i,j,k,l, ind1, ind2, fod_UP, fod_DN    ! loop variables, atomic index, number of UP and DN FODs
         integer :: counter, count_hess, counter2          ! count the number of FODs with similar distance to an atom, use to write fande.out.Counter for entries in hess_D array. counter for tetrahedral stuff
         !
         ! For shellwise OPT
         !
         integer :: object, object_file             ! decide which object shall be optimized, and which one is usedto write agiven file
         character(len=15) :: fdiag, fsearch, fstep ! names of files to distinguish between objects
         real(8) :: r_all_tmp(1), f_all_tmp(1)      ! internal coordinates to feed into LBFGS. One at a time
         !
         ! For free_FOD, partially unconstrained optimization
         !
         logical, allocatable :: freeFOD(:)         ! Store indices of FODs that shall not be constrained
         integer              :: N_free             ! number of unconstrained FODs
         !
         !
         !
         real(8), allocatable :: r_FOD(:,:)         ! FOD positions                                      (former r)
         real(8), allocatable :: r_org_int_UP(:,:)  ! internal FOD positions at used to reconstruct structural motifs
         real(8), allocatable :: r_org_int_DN(:,:)  ! internal FOD positions at used to reconstruct structural motifs
         real(8), allocatable :: f_org_int_UP(:,:)  ! internal FOD forces at used to reconstruct structural motifs
         real(8), allocatable :: f_org_int_DN(:,:)  ! internal FOD forces at used to reconstruct structural motifs
         real(8), allocatable :: r_tmp(:,:)         ! original FOD positions (all) at used to reconstruct structural motifs
         real(8), allocatable :: f_FOD(:,:)         ! FOD forces                                         (former f)
         real(8), allocatable :: r_all(:)           ! FOD positions to be fed into the LBFGS algortihm   (former xvec)
         real(8), allocatable :: f_all(:)           ! FOD forces to be fed into the LBFGS algortihm      (former gvec)
         real(8), allocatable :: r_internal(:,:,:)  ! FOD positions as internal coordinates (one per structural object), 
         real(8), allocatable :: f_internal(:,:,:)  ! FOD forces    as internal coordinates (one per structural object)
         real(8), allocatable :: r_int1(:,:,:)      ! structural element distance, per atom, per spin
         real(8), allocatable :: r_bond(:,:,:)      ! per bonded atoms, per spin - for bond FODs, there is always just one structural object!
         logical, allocatable :: is_bond(:,:)       ! array to define whether a FOD is a bond FOD or not -> preventsincorrect treatmeent later on
         real(8) :: fmax, ftest, dist1, dist2, dist_tmp    ! maximum force, dummy force magnitude, distances of FODs to nuclei
         real(8) :: tmp_vec(3), tmp_vec2(3), tmp_len! temporary vector
         real(8) :: cutoff, junk, cut_bond          ! maximal distance considered for an FOD to be close to an atom, junk variable, cutoff for bond FODs
         real(8) :: f_ave, factor                   ! average force, factor to determine whether force points towards an atom or away from it
         real(8) :: len1, len2, len3, delta         ! lengths of forces etc., difference in angle
         real(8) :: len4, len5, len6                ! evaluation of energy convergence -> reset LBFGS or not
         real(8) :: bond_center(3), fod_center(3)   ! center of two atoms, center for a given set of bond FODs
         ! For rotation matrix
         real(8) :: vec1(3)                         ! Vector to rotate (FOD position 1)
         real(8) :: vec2(3)                         ! Vector to rotate to (other FOD position)    -- overwrite these vectors, such that you don't overwrite the FOD positions/forces
         real(8) :: cross_prod(3)                   ! cross product of the vectors
         real(8) :: cos_angle, angle                ! cosine of the angle between vectors, angle
         real(8) :: full_rot(3,3)                   ! intermediate and full roation matrix
         ! For LBFGS
!         logical :: do_lbfgs                        ! logical to determine whether to do LBFGS or CG     (former LFODLBF)
         double precision :: lbfgs_step, lbfgs_epsilon, lbfgs_precision ! further variables for LBFGS    (former DGUESS, ACCSOLN, XTOL)
         integer :: additional_info(2)              ! output variable for LBFGS algorithm                (former IPRINT)
         integer :: n_corrections                   ! number of corrections in the LBFGS                 (former MUPDATE)
         integer :: n_variables                     ! number of variables in the LBFGS                   (former NPAR)
         integer :: status_error                    ! determine whether there have been any errors       (former iflag)
         real(8), allocatable :: work_space(:)      ! work space used in the LBFGS algorithm             (former work)
         real(8), allocatable :: diag(:)            ! DIAG array used in the LBFGS algorithm             (former DIAG)
         real(8), allocatable :: hess_D(:)          ! HESS_D used in the scaledLBFGS                     (former HESS_D)
         ! For CG only
         real(8) :: gtol, ftol
         integer :: mopt,istat
         real(8) :: scrv(18000)  ! What dimensions are these ??
         ! timings
         real(8) :: start, finish                   ! start and end time for the subroutine
         data first/.true./

         character(6) :: FODFILESTR
         character(2)  :: LETTER
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Get basic parameters for the optimization
         !
         ! Structure of the file FOD_OPT:
         ! Line1: XXX		Optimizer        - options: LBFGS, CG
         ! Line2: YYY		Scaled           - options: yes, no
         ! Line3: ZZZ		Constraint       - options: unconstrained -> full optimization, no constraints
         !                                                  fix1s         -> keeping the 1s FODs fixed
         !                                                  fullForce     -> using internal coordinates, feed them into the optimizer
         !                                                  shellOPT      -> like fullForce, but optimizing internal coordinates individually (beta version)
         !                                                  freeFOD N     -> unconstrain certain FODs, which need to be listed below (i.e., their
         !                                                                   indices). Also, the number of unconstrained FODs N needs to be provided right after the keyword
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!YY This section can be moved to NRLMOL_INPUT.DAT
! 1st variable - LBFGS/CG
! 2nd variable - Y/N scale x/f
! 3rd variable - fix1s/fullForce/shellOPT/freeFOD
         !default
         do_lbfgs   = .true.
         scaled     = .false.
         constraint = 0

         call check_inputs
         do_lbfgs = fod_opt1
         scaled = fod_opt2
         constraint = fod_opt3
         
!        inquire(file='FOD_OPT',exist=exist1)
!        if (exist1) then
!          open(90,file='FOD_OPT',status='old',action='read')
!          ! Read what optimizer shall be used
!          read(90,*) junk1 
!          if (junk1 == 'LBFGS' .or. junk1 == 'lbfgs') then
!            do_lbfgs = .true.
!          else if (junk1 == 'CG' .or. junk1 == 'cg') then
!            do_lbfgs = .false.
!          end if
!          ! Read whether to scale the forces/positions according to Kushanthas scheme
!          read(90,*) junk1
!          if (junk1 == 'no') then
!            scaled = .false.
!          else if (junk1 == 'yes') then
!            scaled = .true.
!          end if
!          ! Read what constraint which shall be used
!          read(90,*) junk1
!          if (junk1 == 'unconstrained') then      ! full optimization, no constraints
!            constraint = 0
!          else if (junk1 == 'fix1s') then         ! fix 1s FODs at nuclear positions
!            constraint = 1
!          else if (junk1 == 'fullForce') then     ! constrain all structural motifs
!            constraint = 2
!          else if (junk1 == 'shellOPT') then       ! fix1s + breathing of inner core. Leave the rest free
!            constraint = 3
!          else if (junk1 == 'freeFOD') then       ! fix1s + breathing of inner core. Leave the rest free
!            constraint = 4
!          end if
!        ! If the file doesn't exist -> do unscaled, unconstrained LBFGS
!        else
!          do_lbfgs   = .true.
!          scaled     = .false.
!          constraint = 0
!        end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Write to screen what is done !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (do_lbfgs) then
           if (scaled) then
             if (constraint == 0) then
               write(6,*) 'Unconstrained FOD-OPT - scaled LBFGS'
             else if (constraint == 1) then
               write(6,*) 'Fix1s constrained FOD-OPT - scaled LBFGS'
             else if (constraint == 2) then
               write(6,*) 'FullForce constrained FOD-OPT - scaled LBFGS'
             else if (constraint == 3) then
               write(6,*) 'ShellOPT constrained FOD-OPT - scaled LBFGS'
             else if (constraint == 4) then
               write(6,*) 'freeFOD constrained FOD-OPT - scaled LBFGS'
             end if
           else
             if (constraint == 0) then
               write(6,*) 'Unconstrained FOD-OPT - LBFGS'
             else if (constraint == 1) then
               write(6,*) 'Fix1s constrained FOD-OPT - LBFGS'
             else if (constraint == 2) then
               write(6,*) 'FullForce constrained FOD-OPT - LBFGS'
             else if (constraint == 3) then
               write(6,*) 'ShellOPT constrained FOD-OPT - LBFGS'
             else if (constraint == 4) then
               write(6,*) 'freeFOD constrained FOD-OPT - LBFGS'
             end if
           end if
         else
           if (scaled) then
             if (constraint == 0) then
               write(6,*) 'Unconstrained FOD-OPT - scaled CG'
             else if (constraint == 1) then
               write(6,*) 'Fix1s constrained FOD-OPT - scaled CG'
             else if (constraint == 2) then
               write(6,*) 'FullForce constrained FOD-OPT - scaled CG'
             else if (constraint == 3) then
               write(6,*) 'ShellOPT constrained FOD-OPT - scaled CG'
             else if (constraint == 4) then
               write(6,*) 'freeFOD constrained FOD-OPT - scaled CG'
             end if
           else
             if (constraint == 0) then
               write(6,*) 'Unconstrained FOD-OPT - CG'
             else if (constraint == 1) then
               write(6,*) 'Fix1s constrained FOD-OPT - CG'
             else if (constraint == 2) then
               write(6,*) 'FullForce constrained FOD-OPT - CG'
             else if (constraint == 3) then
               write(6,*) 'ShellOPT constrained FOD-OPT - CG'
             else if (constraint == 4) then
               write(6,*) 'freeFOD constrained FOD-OPT - CG'
             end if
           end if
         end if
         call flush(6)

         inquire(file='FRMIDT',exist=exist)
         if(.not.exist)then
           write(FODFILESTR,'(A)')'FRMORB'
         else
           write(FODFILESTR,'(A)')'FRMIDT'
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Get atomic information. Read XMOL.xyz file and extract
         ! number of atoms and species. Then, read SYMBOL for 
         ! coordinates (in Bohr)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!YY This block has potential problems. SYMBOL can have more than one
!geometry and > 100 lines. Edited to take geometry from module.
         !NUM_ATMS,XMOL_LIST
         allocate(r_atoms(NUM_ATMS,3))
         allocate(species(NUM_ATMS))
         !print *,"check geometry"
         !checked and it returns the same values in the same order.
         do i=1,NUM_ATMS
          CALL GET_LETTER(XMOL_LIST(i)%ANUM,LETTER)
          species(i)=LETTER
          r_atoms(i,1)=XMOL_LIST(i)%RX
          r_atoms(i,2)=XMOL_LIST(i)%Ry
          r_atoms(i,3)=XMOL_LIST(i)%Rz
          print *,species(i), r_atoms(i,1:3)
         end do

!        open(unit=90,file='XMOL.xyz',status='old',action='read')
!        read(90,*) n_atoms
!        ! Allocate arrays
!        allocate(r_atoms(n_atoms,3))
!        allocate(species(n_atoms))
!        ! Skip one line, then read information
!        read(90,*)
!        do i = 1, n_atoms
!          read(90,*) species(i)
!        end do
!        close(90)
!        ! Get atomic coordinates
!        open(unit=90,file='SYMBOL',status='old',action='read')
!        ! Counter for the atomic index
!        j = 0
!        loop1234: do i = 1, 100    
!          read(90,*) junk1
!          ! Abort if no more coordinates are going to follow
!          if (junk1 == 'ELECTRONS') then
!            exit loop1234
!          else if (junk1(1:3)=='ALL') then
!            ! Go back one line. Read coordinates
!            backspace(90)
!            j = j + 1
!            read(90,*) junk1, junk2, r_atoms(j,1:3)
!          end if
!        end do loop1234
!        close(90)

         !!!!!!!!!!!!!!!!
         ! Start timing
         !!!!!!!!!!!!!!!!
         call cpu_time(start)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Initialize some arrays which are needed later on
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         allocate(n_motifs(n_atoms))
         allocate(r_int1(n_atoms,100,2))        ! initialize size with some number for now
         allocate(r_bond(n_atoms,n_atoms,2))     ! initialize size with some number for now
         r_int1(:,:,:)   = 100.0D0
         r_bond(:,:,:)   = 100.0D0
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Initialize LBFGS or CG parameters
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (do_lbfgs) then
           lbfgs_step = 0.1D0                ! step size in the LBFGS routine.
           additional_info(1) = -1           ! no additional output will be printed
           additional_info(2) =  0           ! type of additional output
           n_corrections = 4                 ! Number of corrections used         -- is this enough ??
           lbfgs_epsilon = 1.0D-6            ! accuracy of the found solution
           lbfgs_precision = 1.0D-16         ! estimated machine precision
         else
           gtol = 0.000001d0
           ftol = 0.000001d0
           mopt = 1000 !1000 => 3000
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Read FRMORB and fforce.dat -> store FOD positions and forces
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         open(91,file=FODFILESTR,status='old',action='read')
         open(92,file='fforce.dat',status='old',action='read')
         read(91,*) fod_UP, fod_DN
         allocate(r_FOD(fod_UP+fod_DN,3))
         allocate(r_tmp(fod_UP+fod_DN,3))
         allocate(f_FOD(fod_UP+fod_DN,3))
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Read FOD positions and forces into the respective arrays
         !
         ! Idea: Read the forces in with a format, like 3F20.5
         ! -> constrain used forces, avoid numerical issues?
         ! Coordinates with 3F10.7 -> write them later like this as well
         ! 
         ! To do this, first read them unformatted. Write the files 
         ! again, with the right format. Then, read them in in this
         ! format. There's gotta be a better way, but it should work for
         ! now.
         !
         ! Thanks to Kushantha
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do i = 1, fod_UP+fod_DN
           read(91,*) r_FOD(i,1:3)
           read(92,*) f_FOD(i,1:3)
           r_tmp(i,1:3) = r_FOD(i,1:3)
         end do
         close(91)
         close(92)
         print *, 'in electronic geometry'
         print *, 'total fod', fod_UP+fod_DN
         print *, 'FODs'
         do i = 1,fod_UP+fod_DN
           print *, (r_FOD(i,j),j=1,3)
         end do
         print *, 'FOD forces'
         do i = 1,fod_UP+fod_DN
           print *, (f_FOD(i,j),j=1,3)
         end do
         if ((constraint >= 2)) then
           open(91,file=FODFILESTR,status='unknown',action='write')
           open(92,file='fforce.dat',status='unknown',action='write')
           write(91,*) fod_UP, FOD_DN
           do i = 1, fod_UP+fod_DN
             write(91,fmt='(3F12.7)') r_FOD(i,1:3)
             write(92,fmt='(3F12.5)') f_FOD(i,1:3)
           end do
           close(91)
           close(92)

           open(91,file=FODFILESTR,status='old',action='read')
           open(92,file='fforce.dat',status='old',action='read')
           read(91,*) fod_UP, FOD_DN
           do i = 1, fod_UP+fod_DN
             read(91,fmt='(3F12.7)') r_FOD(i,1:3)
             read(92,fmt='(3F12.5)') f_FOD(i,1:3)
             r_tmp(i,1:3) = r_FOD(i,1:3)
           end do
           close(91)
           close(92)
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! For partially unconstrained OPT (freeFOD)
         !
         if (constraint == 4) then
           allocate(freeFOD(fod_UP+fod_DN))
           freeFOD(:) = .false.
           open(90,file='FOD_OPT',status='old',action='read')
           ! skip the first two lines
           read(90,*)
           read(90,*)
           read(90,*) junk1, N_free
           do i = 1, N_free
             read(90,*) ind1
             freeFOD(ind1) = .true.
           end do
           close(90)
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Initialize internal coordinate arrays. Maximum size: Number
         ! of FODs (maximum number of internal coordinates = number of
         ! FODs, if all FODs are independent)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         allocate(r_internal(fod_UP+fod_DN,3,2))        ! initialize size, with 3 elements each (x,y,z), for each spin (1,2)
         allocate(f_internal(fod_UP+fod_DN,3,2))        ! initialize size, with 3 elements each (x,y,z), for each spin (1,2)
         allocate(is_bond(fod_UP+fod_DN,2))             ! initialize size, for each spin (1,2)
         r_internal(:,:,:) = 1000.0D0
         f_internal(:,:,:) = 1000.0D0
         is_bond(:,:) = .false.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! If scaled is enabled:
         ! Create Table file if it doesn't exist yet.
         ! And setting up HESS_DIAG file to scale forces and coordinates
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (scaled) then
           call write_table_constr ! included this in this subroutine here. no need for an extra subroutine
           inquire(file='HESS_DIAG',exist=exist1)
           if (.not.exist1) then
             !
             ! If HESS_DIAG doesn't exist -> create it
             !
             call write_HD_constr()
           end if
           !
           ! Allocate hess_D array -> storing the scaling factors
           !
           allocate(hess_D(3*(fod_UP+fod_DN)))
           !
           ! open HESS_DIAG file
           !
           open(95,file='HESS_DIAG',status='old',action='read')
           !
           ! Write the corresponding array
           !
           count_hess = 0
           do i = 1, fod_UP+fod_DN
             do j = 1, 3
               count_hess = count_hess + 1
               read(95,*) hess_D(count_hess)
             end do
           end do
           close(95)
         end if
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! separate things for unconstrained, fix1s and fullForce
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! If unconstrained -> just take all coordinates and forces
         !
         if (constraint == 0) then
           do i = 1, fod_UP
             r_internal(i,1:3,1) = r_FOD(i,1:3)
             f_internal(i,1:3,1) = f_FOD(i,1:3)
           end do
           do i = 1, fod_DN
             r_internal(i,1:3,2) = r_FOD(i+fod_UP,1:3)
             f_internal(i,1:3,2) = f_FOD(i+fod_UP,1:3)
           end do
         !
         ! In any other case (fix1s, fullForce)
         !
         else

!
! 1. Loop over FODs -> evaluate distance to nearest nucleus
! 2. Loop over FODs -> evaluate FODs with the same distance to that nucleus
! for all of these  -> evaluate force per structual motifs
!
! Find bond FODs -> symmetrize them as well
!
           if (n_atoms == 1) cutoff = 50.0D0      ! For atoms -> take all FODs into account
           if (n_atoms > 1)  cutoff =  2.5D0      ! For molecules -> Only take FODs close to one atom   make this larger for bigger atoms?
           if (n_atoms > 1)  cut_bond =  1.50D0   ! For molecules -> Only take FODs to be in a bond which have a distance difference to two atoms of less than this

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! 1. loop over FODs - UP channel !
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !
           ! determine which atom the FOD is closest to
           !
           ! First of all, get all bond FODs
           ! 
           do i = 1, fod_UP
             dist1 = 1000.0                  ! initial value for the distance. Going to find minimum for all atoms
             do j = 1, n_atoms                ! Loop over atoms
               dist_tmp = sqrt(sum((r_FOD(i,:)-r_atoms(j,:))**2))
               if ((dist_tmp < dist1) .and. (dist_tmp < cutoff)) then  ! if distance to atom b is smaller than dist1 AND is is smaller than a cutoff radius
                 dist1 = dist_tmp                                      ! the cutoff radius ensures that we are only taking FODs for one atom 
                 ind1 = j
               end if
             end do
             !
             ! For bond FODs -> get FOD information
             ! Count how many bond FODs there are
             ! 
             ind2 = 10000
             counter = 0
             if (n_atoms > 1 .and. r_internal(i,3,1).ne.100.0D0) then ! exclude FODs which have already been found as bond FODs
               bond_center(:) = (/ 0.0D0, 0.0D0, 0.0D0 /)
               fod_center(:) = (/ 0.0D0, 0.0D0, 0.0D0 /)
               do j = 1, n_atoms                ! Loop over atoms
                 dist_tmp = sqrt(sum((r_FOD(i,:)-r_atoms(j,:))**2))
                 !
                 ! For X-H bonds: see whether there are any H in the
                 ! bond. If so -> single bond FODs!
                 !
                 ! If either is H
                 ! atom j is NOT atom ind1
                 ! FOD i it is not already a bond
                 ! FOD is not a 1s for either of the atoms
                 !
                 if (((trim(adjustl(species(ind1)))=='H').or.
     &(trim(adjustl(species(j)))=='H'))
     &.and.(j/=ind1).and.(is_bond(i,1).eqv..false.)
     &.and.(dist1>0.01).and.(dist_tmp>0.01)) then
                   !
                   ! Check whether atom j is the closest atom to atom
                   ! ind1. If so, continue here. If not -> do not do
                   ! anything
                   !
                   counter2 = 0
                   do k = 1, n_atoms
                     if ((k.ne.j).and.(k.ne.ind1)) then
                       !
                       ! If any atom is closer to atom ind1
                       ! -> abort
                       !
                       if(sqrt(sum((r_atoms(k,:)-r_atoms(ind1,:))**2))<
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))) then
                         counter2 = 1
                       end if
                     end if
                   end do
                   !
                   ! Only if atom j is closest to atom ind1->continue 
                   !
                   if (counter2 == 0) then
                     !
                     ! Exclude evaluation where the FOD has a larger
                     ! distance to any of the two atoms than the distance
                     ! between the two atoms
                     !
                     if ((dist1 > 
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))).or.((dist_tmp > 
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))))) then
                     else
                       ! 
                       ! Condition for X-H bonds:
                       !  1. Vectors X-FOD and H-FOD are parallel to the
                       !     bond axis
                       !  2. The direction of these vectors is exactly
                       !     parallel to the bond axis (by construction)
                       !
                       ! Vector between atoms. Pointing from initial to j
                       ! Distance between atoms.
                       ! Normalize the vector
                       !
                       tmp_vec = r_atoms(ind1,:)-r_atoms(j,:)
                       len3=sqrt(sum((r_atoms(ind1,:)-r_atoms(j,:))**2))
                       tmp_vec = tmp_vec/len3
                       !
                       ! Vector between FOD and first atom
                       ! And its length
                       ! Normalize it
                       !
                       tmp_vec2 = r_FOD(i,:)-r_atoms(ind1,:)
                       len2 = sqrt(sum((r_FOD(i,:)-r_atoms(ind1,:))**2))
                       tmp_vec2 = tmp_vec2/len2
                       !
                       ! Form dot porduct to bond vector between atoms
                       ! Gives the cosine of the angle between them
                       !
                       len2 = dot_product(tmp_vec,tmp_vec2)
                       !
                       ! The cosine needs to be 1, or very close to it 
                       ! Avoid numercial noise, by setting the length to 1
                       ! if is is numerically larger
                       !
                       if (abs(len2) > 1.00D0) len2 = 1.0D0
                       !
                       ! DO the same for the second atom
                       !
                       tmp_vec2 = r_FOD(i,:)-r_atoms(j,:)
                       len3 = sqrt(sum((r_FOD(i,:)-r_atoms(j,:))**2))
                       tmp_vec2 = tmp_vec2/len3
                       !
                       len3 = dot_product(tmp_vec,tmp_vec2)
                       !
                       if (abs(len3) > 1.00D0) len3 = 1.0D0
                       !
                       ! The cosines needs to be 1, or very close to it 
                       !
                       if ((abs(len2)>0.99).and.(abs(len3)>0.99)) then
                         !
                         counter = 1
                         r_internal(i,3,1) = 100.0D0
                         is_bond(i,1) = .true.
                       end if

                       if (counter == 1) then
                         !
                         ! Use force component along the FOD-atom axis.
                         ! -> use: F_i \cdot (a_i-R_i) / |a_i-R_i|
                         ! -> Store as internal coordinates
                         !
                         len1 = dot_product(r_FOD(i,:)-r_atoms(ind1,:),
     &f_FOD(i,:))/sqrt(sum(((r_FOD(i,:)-r_atoms(ind1,:))**2)))
                         !
                         ! If the force along the bond is tiny -> numerical
                         ! 'noise' ->> ignore. Or we are in a homonuclear
                         ! bonding situation!
                         ! 
                         ! HERE: make this 1.0D-4. For X-H bonds
                         ! only
                         !
                         if (abs(len1) < 1.0D-4) then
                           f_internal(i,1,1) = 0.0D0
                         else
                           f_internal(i,1,1) = len1
                         end if
                         !
                         ! Use distance of FOD to the atom as coordinate
                         ! -> Store as internal coordinates
                         !
                         len1 = 
     &sqrt(sum(((r_FOD(i,:)-r_atoms(ind1,:))**2)))
                         r_internal(i,1,1) = len1
                         !
                         ! If a coordinate and a force similar to these ones
                         ! has been found before -> use that one instead !
                         ! Go through the other internal coordinates and
                         ! check. UP or DN channel!
                         !
                         ! UP Channel
                         ! If force is 0 -> don't do this.   TBD
                         !
                         do l = 1, fod_UP
                           if (l.ne.i) then
                             !
                             ! small subroutine comparing and averaging
                             !
                             call comp_int(r_internal(i,1,1),
     &r_internal(l,1,1),f_internal(i,1,1),f_internal(l,1,1))
                           end if
                         end do
                         !
                         ! DN Channel
                         !
                         do l = 1, fod_DN
                           call comp_int(r_internal(i,1,1),
     &r_internal(l,1,2),f_internal(i,1,1),f_internal(l,1,2))
                         end do
                         !
                         ! overwrite any existing other values for this 
                         ! internal coordinate pair. Ensure that only the
                         ! ones above are used
                         !
                         r_internal(i,2,1) = 1000.0D0
                         r_internal(i,3,1) = 1000.0D0
                         f_internal(i,2,1) = 1000.0D0
                         f_internal(i,3,1) = 1000.0D0
                       end if
                     end if
                   end if
                   !!!!!
                   ! For any other bonding situation
                   !!!!!
                 else
                   !
                   ! if distance to atom j is equal to dist1 : same distance to two atoms -> bond FOD
                   !
                   if ((abs(dist_tmp-dist1)<cut_bond).and.(j/=ind1)
     &.and.(is_bond(i,1).eqv..false.))then
                     ind2 = j
                     counter = 0
                     !
                     ! Get center between atoms
                     !
                     bond_center(:) = 
     &(r_atoms(ind1,:)+r_atoms(ind2,:))/2.0D0
                     !
                     ! Get distance of initial FOD to this bond center
                     !
                     dist2 = sqrt(sum((r_FOD(i,:)-bond_center(:))**2))
                     !
                     ! Use the center between the atoms. Analyze
                     ! distance of the FODs to that center. They should
                     ! be very similar!
                     !
                     do k = 1, fod_UP
                       !
                       ! Check distance to FOD center
                       !
                       dist_tmp=
     &sqrt(sum((r_FOD(k,:)-bond_center(:))**2))
                       if ((abs(dist_tmp-dist2) < 0.05D0)) then
                         !
                         ! Check distance to corresponding atom
                         !
                         dist_tmp=
     &sqrt(sum((r_FOD(k,:)-r_atoms(ind1,:))**2))
                         if ((abs(dist_tmp-dist1) < 0.05D0)) then
                           !
                           ! Get information about the center of all bond
                           ! FODs
                           !
                           fod_center(:) = fod_center(:)+r_FOD(k,:)
                           !
                           ! Increase the counter for the bond FODs
                           !
                           counter = counter + 1
                           !
                           ! Set internal z-coordinates for these FODs to 100.0D0
                           ! -> differentiate them later on. Do not evaluate
                           ! these FODs later again
                           !
                           r_internal(k,3,1) = 100.0D0
                           is_bond(k,1) = .true.
                         end if
                       end if
                     end do
                     !
                     ! Calculate the actual bond center
                     !
                     fod_center(:) = fod_center(:)/counter
                     !
                     ! If counter = 1: Single bond. Project force onto the
                     ! bond axis. Thus, only move the FOD on the axis
                     ! Effectively: Project FOD force on FOD coordinate !
                     !
                     ! Use only ONE coordinates as internal one,
                     ! not three
                     !
                     if (counter == 1) then
                       !
                       ! Use force component along the FOD-atom axis.
                       ! -> use: F_i \cdot (a_i-R_i) / |a_i-R_i|
                       ! -> Store as internal coordinates
                       !
                       len1 = dot_product(r_FOD(i,:)-r_atoms(ind1,:),
     &f_FOD(i,:))/sqrt(sum(((r_FOD(i,:)-r_atoms(ind1,:))**2)))
                       !
                       ! If the force along the bond is tiny -> numerical
                       ! 'noise' ->> ignore. Or we are in a homonuclear
                       ! bonding situation!
                       !
                       if (abs(len1) < 2.5D-4.or.
     &(species(ind1)==species(ind2))) then
                         f_internal(i,1,1) = 0.0D0
                       else
                         f_internal(i,1,1) = len1
                       end if                     
                       !
                       ! Use distance of FOD to the atom as coordinate
                       ! -> Store as internal coordinates
                       !
                       len1=sqrt(sum(((r_FOD(i,:)-r_atoms(ind1,:))**2)))
                       r_internal(i,1,1) = len1
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       !
                       ! UP Channel
                       ! If force is 0 -> don't do this.   TBD
                       !
                       do l = 1, fod_UP
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,1),
     &r_internal(l,1,1),f_internal(i,1,1),f_internal(l,1,1))
                         end if
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         call comp_int(r_internal(i,1,1),
     &r_internal(l,1,2),f_internal(i,1,1),f_internal(l,1,2))
                       end do
                       !
                       ! overwrite any existing other values for this 
                       ! internal coordinate pair. Ensure that only the
                       ! ones above are used
                       !
                       r_internal(i,2,1) = 1000.0D0
                       r_internal(i,3,1) = 1000.0D0
                       f_internal(i,2,1) = 1000.0D0
                       f_internal(i,3,1) = 1000.0D0
                     !
                     ! If counter > 1: Multiple bond. Use force
                     ! perpendicular to bond axis, along (bond-center)-FOD
                     ! axis. 
                     ! Bond-center = center of all bond-FODs!!
                     ! Also, move along bond as well if heteronuclear
                     !
                     ! For homonuclear bonds -> keep it only perpendicular
                     ! For heteronuclear bonds -> include components along
                     ! the bond axis
                     !
                     else
                       !
                       ! Average forces for this structural element. 
                       ! More consistent
                       ! Analyze which FODs are in this bond
                       ! Use force component along the bond_center-FOD axis.
                       ! -> use: F_i \cdot (a_i-bond_center) / |a_i-bond_center|
                       !
                       len1 = 0.0D0
                       do k = 1, fod_UP
                         dist_tmp=
     &sqrt(sum((r_FOD(k,:)-bond_center(:))**2))
                         if ((abs(dist_tmp-dist2) < 0.05D0)) then
                           dist_tmp=
     &sqrt(sum((r_FOD(k,:)-r_atoms(ind1,:))**2))
                           if ((abs(dist_tmp-dist1) < 0.05D0)) then
                             len2 = dot_product(f_FOD(k,:),
     &r_FOD(k,:)-fod_center(:))/
     &sqrt(sum(((r_FOD(k,:)-fod_center(:))**2)))
                             !
                             ! add it up
                             !
                             len1 = len1 + len2
                           end if
                         end if
                       end do                     
                       len1 = len1/counter
                       !
                       ! Store internal force
                       !
                       f_internal(i,1,1) = len1
                       !
                       ! Store corresponding internal coordinate
                       ! (same for all points in this structural motif)
                       !
                       len1 = sqrt(sum(((r_FOD(i,:)-fod_center(:))**2)))
                       r_internal(i,1,1) = len1
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       ! 
                       ! UP Channel
                       do l = 1, fod_UP
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,1),
     &r_internal(l,1,1),f_internal(i,1,1),f_internal(l,1,1))
                         end if
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         call comp_int(r_internal(i,1,1),
     &r_internal(l,1,2),f_internal(i,1,1),f_internal(l,1,2))
                       end do
                       !
                       ! IF we have a homonuclear bonding situation
                       ! -> ONLY move FODs perpendicular to bond axis
                       !
                       if (species(ind1)==species(ind2)) then
                         r_internal(i,2,1) = 1000.0D0
                         f_internal(i,2,1) = 1000.0D0
                       ! 
                       ! IF heternonuclear bonding situation
                       ! -> take component along the bond axis
                       ! as well. Define one vector including
                       ! both components
                       ! 
                       else
                         !
                         ! Average forces for this structural element. 
                         ! More consistent
                         ! Analyze which FODs are in this bond
                         ! Use force component along the bond_center-atom
                         ! axis.
                         ! -> use: F_i \cdot (bond_center-atom)/|bond_center-atom|
                         !
                         len2 = 0.0D0
                         do k = 1, fod_UP
                           dist_tmp=
     &sqrt(sum((r_FOD(k,:)-bond_center(:))**2))
                           if ((abs(dist_tmp-dist2) < 0.05D0)) then
                             dist_tmp=
     &sqrt(sum((r_FOD(k,:)-r_atoms(ind1,:))**2))
                             if ((abs(dist_tmp-dist1) < 0.05D0)) then
                               len3 = dot_product(f_FOD(k,:),
     &fod_center(:)-r_atoms(ind1,:))/   ! fod_center -> atom
     &sqrt(sum(((fod_center(:)-r_atoms(ind1,:))**2)))
                               !
                               ! add it up
                               !
                               len2 = len2 + len3
                             end if
                           end if
                         end do
                         len2 = len2/counter
                         !
                         ! store internal coordinate
                         !
                         f_internal(i,2,1) = len2
                         !
                         ! Store corresponding internal coordinate
                         ! same for all points in this structural motif
                         !
                         len2=
     &sqrt(sum(((fod_center(:)-r_atoms(ind1,:))**2)))
                         r_internal(i,2,1) = len2
                         !
                         ! If a coordinate and a force similar to these ones
                         ! has been found before -> use that one instead !
                         ! Go through the other internal coordinates and
                         ! check. UP or DN channel!
                         !
                         ! UP Channel 
                         do l = 1, fod_UP
                           if (l.ne.i) then
                             call comp_int(r_internal(i,2,1),
     &r_internal(l,2,1),f_internal(i,2,1),f_internal(l,2,1))
                           end if
                         end do
                         ! DN Channel 
                         do l = 1, fod_DN
                           call comp_int(r_internal(i,2,1),
     &r_internal(l,2,2),f_internal(i,2,1),f_internal(l,2,2))
                         end do
                       end if
                       r_internal(i,3,1) = 1000.0D0
                       f_internal(i,3,1) = 1000.0D0
                     end if
                   end if
                 end if
               end do

               if (counter > 0) then
                 write(junk1,'(I4.4)') i
                 write(junk2,'(I4.4)') counter
                 if (counter == 1) then ! single FOD
                   write(6,777) 'There is   ',trim(junk2),' UP-bond 
     &FOD  in the structural motif of UP-bond FOD ',trim(junk1)
                 else                   ! more FODs
                   write(6,777) 'There are  ',trim(junk2),' UP-bond 
     &FODs in the structural motif of UP-bond FOD ',trim(junk1)
                 end if
               end if
             end if
           end do




           !!!!!!!!!!!!!
           ! Now, do all remaining FODs (lone/core)
           !!!!!!!!!!!!!!1
           n_motifs(:) = 0
           do i = 1, fod_UP
             dist1 = 1000.0                  ! initial value for the distance. Going to find minimum for all atoms
             do j = 1, n_atoms                ! Loop over atoms
               dist_tmp = sqrt(sum((r_FOD(i,:)-r_atoms(j,:))**2))
               if ((dist_tmp < dist1) .and. (dist_tmp < cutoff)) then  ! if distance to atom b is smaller than dist1 AND is is smaller than a cutoff radius
                 dist1 = dist_tmp                                      ! the cutoff radius ensures that we are only taking FODs for one atom 
                 ind1 = j
               end if
             end do
             !
             ! fix1s -> set 1s forces to 0
             !
             if (dist1 < 0.01) then
               if (constraint == 4) then
                 if (freeFOD(i).eqv..false.) then
                   f_FOD(i,:) = (/ 0.0D0, 0.0D0, 0.0D0 /)
                   write(junk1,'(I4.4)') i
                   write(6,778) 'UP-FOD no. ',trim(junk1),
     &' will be fixed (1s FOD)'
                 end if
               else
                 f_FOD(i,:) = (/ 0.0D0, 0.0D0, 0.0D0 /)
                 write(junk1,'(I4.4)') i
                 write(6,778) 'UP-FOD no. ',trim(junk1),
     &' will be fixed (1s FOD)'
               end if
             else
               !
               ! Check whether it is a bond or not
               ! If not -> keep going
               ! Make sure it is not a bond FOD!
               !
               if (is_bond(i,1).eqv..false.) then
                 !
                 ! Only if the fullForce contraint is suppossed to be used: Continue here
                 !
                 if ((constraint >= 2)) then 
                   !
                   ! If core/lone FOD -> symmetrize the core FODs accordingly
                   !
                   write(junk1,'(I4.4)') i
                   write(junk2,'(I4.4)') ind1
                   write(6,777) 'UP-FOD no. ',trim(junk1),
     &' is core/lone FOD in atom ',trim(junk2)
                   !
                   ! If the structural motif has already been found (distance similar to the one found before within 0.033 bohr) -> don't do everything again
                   !
                   exist = .false.
                   do l = 1, 100
                     if (abs(dist1 - r_int1(ind1,l,1)) < 0.033) then ! any distance FOR THIS ATOM. SPIN UP. 
                       exist = .true.
                       exit
                     end if
                   end do
                   !
                   ! If not -> new structural element
                   ! If structural object is new AND it has not been assigned before (i.e.
                   ! it is not in a bond) 
                   !
                   if (.not.exist) then
                     n_motifs(ind1) = n_motifs(ind1) + 1
                     r_int1(ind1,n_motifs(ind1),1) = dist1    ! Store the current identifier for the structural motif.

                     r_internal(i,1:3,1) = r_FOD(i,1:3)
                     f_internal(i,1:3,1) = f_FOD(i,1:3)
                     !
                     ! Use total force on FOD_i for all other FODs
                     ! Store it as internal coordinate
                     ! Resymmetrize the other FODs after optimization
                     !
                     ! HERE: just count how many FODs are in that
                     ! structural motif
                     !
                     ! Exclude already assigned FODs
                     !
                     fod_center = (/ 0.0D0, 0.0D0, 0.0D0 /)
                     counter = 0
                     do k = 1, fod_UP
                       if (is_bond(k,1).eqv..false.) then               ! if this FOD has not been assigned yet (no bond)
                         dist2 = 
     &sqrt(sum((r_FOD(k,:)-r_atoms(ind1,:))**2))                        ! check distance to the atom in question
                         if (abs(dist1 - dist2) < 0.033) then           ! If the new distance is the same as the old one. 
                           counter = counter + 1
                           fod_center(:) = fod_center(:)+r_FOD(k,:)
                         end if
                       end if
                     end do
                     fod_center(:) = fod_center(:)/counter
                     !
                     ! write to screen
                     ! 
                     write(junk1,'(I4.4)') i
                     write(junk2,'(I4.4)') counter
                     if (counter == 1) then ! single FOD
                       write(6,777) 'There is   ',trim(junk2),
     &' core/lone UP-FOD  in the motif of UP-FOD ',trim(junk1)
                     else                   ! more FODs
                       write(6,777) 'There are  ',trim(junk2),
     &' core/lone UP-FODs in the motif of UP-FOD ',trim(junk1)
                     end if
                     !
                     ! If counter = 1: Single FOD (like single lone FOD). adjust force to
                     ! reflect symmetry of the system. Take atom-FOD
                     ! axis only and project FOD force on this axis!
                     ! Effectively: Project FOD force on FOD coordinate!
                     !
                     ! Use only ONE coordinates as internal one,
                     ! not three
                     !
                     if (counter == 1) then
                       !
                       ! Use force component along the atom-FOD axis.
                       ! -> use: F_i \cdot (a_i-R_i) / |a_i-R_i|
                       ! -> Store as internal coordinates
                       !
                       len1 = dot_product(r_FOD(i,:)-r_atoms(ind1,:),
     &f_FOD(i,:))/sqrt(sum(((r_FOD(i,:)-r_atoms(ind1,:))**2)))
                       !
                       f_internal(i,1,1) = len1
                       !
                       ! Use distance of FOD to the atom as coordinate
                       ! -> Store as internal coordinates
                       !
                       len1=sqrt(sum(((r_FOD(i,:)-r_atoms(ind1,:))**2)))
                       r_internal(i,1,1) = len1
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       !
                       ! UP Channel
                       do l = 1, fod_UP
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,1),
     &r_internal(l,1,1),f_internal(i,1,1),f_internal(l,1,1))
                         end if
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         call comp_int(r_internal(i,1,1),
     &r_internal(l,1,2),f_internal(i,1,1),f_internal(l,1,2))
                       end do
                       !
                       ! overwrite any existing other values for this 
                       ! internal coordinate pair. Ensure that only the
                       ! ones above are used
                       !
                       r_internal(i,2,1) = 1000.0D0
                       r_internal(i,3,1) = 1000.0D0
                       f_internal(i,2,1) = 1000.0D0
                       f_internal(i,3,1) = 1000.0D0
                     end if
                     !
                     ! If there are more than 1, but not 4
                     !
                     if ((counter > 1).and.(counter.ne.4)) then
                       !
                       ! Average forces for this structural element. 
                       ! More consistent
                       ! Analyze which FODs are in this lone region
                       ! Use force component along the bond_center-FOD axis.
                       ! -> use: F_i \cdot (a_i-bond_center) / |a_i-bond_center|
                       !
                       len1 = 0.0D0
                       do k = 1, fod_UP
                         if (is_bond(k,1).eqv..false.) then             ! if this FOD has not been assigned yet (no bond)
                           dist2 = 
     &sqrt(sum((r_FOD(k,:)-r_atoms(ind1,:))**2))                        ! check distance to the atom in question
                           if (abs(dist1 - dist2) < 0.033) then         ! If the new distance is the same as the old one.
                             len2 = dot_product(f_FOD(k,:),
     &r_FOD(k,:)-fod_center(:))/
     &sqrt(sum(((r_FOD(k,:)-fod_center(:))**2)))
                             !
                             ! add it up
                             !
                             len1 = len1 + len2
                           end if
                         end if
                       end do
                       len1 = len1/counter
                       !
                       ! Store internal force
                       !
                       f_internal(i,1,1) = len1
                       !
                       ! Store corresponding internal coordinate
                       ! (same for all points in this structural motif)
                       !
                       len1 = sqrt(sum(((r_FOD(i,:)-fod_center(:))**2)))
                       r_internal(i,1,1) = len1
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       ! 
                       ! UP Channel
                       do l = 1, fod_UP
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,1),
     &r_internal(l,1,1),f_internal(i,1,1),f_internal(l,1,1))
                         end if
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         call comp_int(r_internal(i,1,1),
     &r_internal(l,1,2),f_internal(i,1,1),f_internal(l,1,2))
                       end do
                       !
                       ! Average forces for this structural element. 
                       ! More consistent
                       ! Analyze which FODs are in this lone region
                       ! Use force component along the fod_center-atom
                       ! axis.
                       ! -> use: F_i \cdot (bond_center-atom)/|bond_center-atom|
                       !
                       ! Do this only is there are several atoms
                       !
                       if (n_atoms > 1) then
                         len2 = 0.0D0
                         do k = 1, fod_UP
                           if (is_bond(k,1).eqv..false.) then          ! if this FOD has not been assigned yet (no bond)
                             dist2=
     &sqrt(sum((r_FOD(k,:)-r_atoms(ind1,:))**2))                       ! check distance to the atom in question
                             if (abs(dist1 - dist2) < 0.033) then      ! If the new distance is the same as the old one 
                               len3 = dot_product(f_FOD(k,:),
     &fod_center(:)-r_atoms(ind1,:))/
     &sqrt(sum(((fod_center(:)-r_atoms(ind1,:))**2)))
                               !
                               ! add it up
                               !
                               len2 = len2 + len3
                             end if
                           end if
                         end do
                         len2 = len2/counter
                         !
                         ! Store internal force
                         !
                         f_internal(i,2,1) = len2
                         !
                         ! Store corresponding internal coordinate
                         ! (same for all points in this structural motif)
                         !
                         len2=
     &sqrt(sum(((r_atoms(ind1,:)-fod_center(:))**2)))
                         r_internal(i,2,1) = len2
                         !
                         ! If a coordinate and a force similar to these ones
                         ! has been found before -> use that one instead !
                         ! Go through the other internal coordinates and
                         ! check. UP or DN channel!
                         !
                         ! UP Channel 
                         do l = 1, fod_UP
                           if (l.ne.i) then
                             call comp_int(r_internal(i,2,1),
     &r_internal(l,2,1),f_internal(i,2,1),f_internal(l,2,1))
                           end if
                         end do
                         ! DN Channel 
                         do l = 1, fod_DN
                           call comp_int(r_internal(i,2,1),
     &r_internal(l,2,2),f_internal(i,2,1),f_internal(l,2,2))
                         end do
                       !
                       ! If only one atom -> ignore this
                       !
                       else
                         r_internal(i,2,1) = 1000.0D0
                         f_internal(i,2,1) = 1000.0D0
                       end if

                       r_internal(i,3,1) = 1000.0D0
                       f_internal(i,3,1) = 1000.0D0
                     end if
                     !
                     ! If there are four -> take breathing force and 
                     ! only the distance to the center as coordinates
                     ! Average all forces along their FODs to get 
                     ! Consistent force for all FODs
                     !
                     ! assign the internal coordinates
                     !
                     if (counter == 4) then
                       r_internal(i,1,1) = dist1
                       r_internal(i,2,1) = 1000.0D0
                       r_internal(i,3,1) = 1000.0D0
                       f_internal(i,1,1) = 0.0D0
                       f_internal(i,2,1) = 1000.0D0
                       f_internal(i,3,1) = 1000.0D0
                       do k = 1, fod_UP
                         if ((is_bond(k,1).eqv..false.)) then
                           dist2 = 
     &sqrt(sum((r_FOD(k,:)-r_atoms(ind1,:))**2))! check distance to the atom in question
                           if (abs(dist1 - dist2) < 0.033) then         ! If the new distance is the same as the old one
                             !
                             ! Take the force along the atom-FOD axis
                             !
                             len1 = dot_product(f_FOD(k,:),
     &r_FOD(k,:)-r_atoms(ind1,:))/
     &sqrt(sum(((r_FOD(k,:)-r_atoms(ind1,:))**2)))
                             !
                             ! add it up
                             !
                             f_internal(i,1,1) = f_internal(i,1,1)+len1
                           end if
                         end if
                       end do
                       f_internal(i,1,1) = f_internal(i,1,1)/4.0D0
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       ! Need to be averaged: TBD
                       !
                       ! UP Channel
                       do l = 1, fod_UP
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,1),
     &r_internal(l,1,1),f_internal(i,1,1),f_internal(l,1,1))
                         end if
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         call comp_int(r_internal(i,1,1),
     &r_internal(l,1,2),f_internal(i,1,1),f_internal(l,1,2))
                       end do


                     end if
                   end if   ! end check whether structural element was already found
                 end if     ! end check what contraint to be used
               end if       ! end check for 1s FODs
             end if         ! ech check whether it is a bond FOD or not
           end do         ! end 1. loop

!!!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!           !
!!!           ! If only fix1s -> write cooridnates and forces right here
!!!           ! Forces on 1s are 0
!!!           ! 
!!!           if (constraint == 1) then
!!!             do k = 1, fod_UP
!!!               !
!!!               ! Do not take the 1s into the optimization
!!!               !
!!!               if ((f_FOD(k,1).ne.0.0D0)
!!!     &.and.(f_FOD(k,2).ne.0.0D0).and.(f_FOD(k,3).ne.0.0D0)) then
!!!                 r_internal(k,1:3,1) = r_FOD(k,1:3)
!!!                 f_internal(k,1:3,1) = f_FOD(k,1:3)
!!!               end if
!!!             end do
!!!           end if
!!!
!!!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!           ! If freeFOD constraint -> restore the actual forces for the
!!!           ! free FODs
!!!           !!!!!!!!!!!!!!!1
!!!           if (constraint == 4) then
!!!             do k = 1, fod_UP
!!!               !
!!!               ! Do not take the 1s into the optimization
!!!               !
!!!               if ((freeFOD(k).eqv..true.)) then
!!!                 r_internal(k,1:3,1) = r_FOD(k,1:3)
!!!                 f_internal(k,1:3,1) = f_FOD(k,1:3)
!!!               end if
!!!             end do
!!!           end if
             








           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! 1. loop over FODs - DN channel !
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !
           ! First, get all bond FODs
           !
           do i = 1, fod_DN
             dist1 = 1000.0                  ! initial value for the distance. Going to find minimum for all atoms
             do j = 1, n_atoms                ! Loop over atoms
               dist_tmp = sqrt(sum((r_FOD(i+fod_UP,:)-r_atoms(j,:))**2))
               if ((dist_tmp < dist1) .and. (dist_tmp < cutoff)) then  ! if distance to atom b is smaller than dist1 AND is is smaller than a cutoff radius
                 dist1 = dist_tmp                                      ! the cutoff radius ensures that 
                 ind1 = j
               end if
             end do
             !
             ! For bond FODs -> get FOD information
             ! Count how many bond FODs there are
             ind2 = 10000
             counter = 0
             if (n_atoms > 1 .and. r_internal(i,3,2).ne.100.0D0) then ! exclude FODs which have already been assigned as a bond FOD
               bond_center(:) =(/ 0.0D0, 0.0D0, 0.0D0 /)
               fod_center(:) =(/ 0.0D0, 0.0D0, 0.0D0 /)
               do j = 1, n_atoms                ! Loop over atoms
                 dist_tmp=sqrt(sum((r_FOD(i+fod_UP,:)-r_atoms(j,:))**2))
                 !
                 ! For X-H bonds: see whether there are any H in the
                 ! bond. If so -> single bond FODs!
                 !
                 ! If either is H
                 ! atom j is NOT atom ind1
                 ! FOD i it is not already a bond
                 ! FOD is not a 1s for either of the atoms
                 !
                 if (((trim(adjustl(species(ind1)))=='H').or.
     &(trim(adjustl(species(j)))=='H'))
     &.and.(j/=ind1).and.(is_bond(i,2).eqv..false.)
     &.and.(dist1>0.01).and.(dist_tmp>0.01)) then
                   !
                   ! Check whether atom j is the closest atom to atom
                   ! ind1. If so, continue here. If not -> do not do
                   ! anything
                   !
                   counter2 = 0
                   do k = 1, n_atoms
                     if ((k.ne.j).and.(k.ne.ind1)) then
                       !
                       ! If any atom is closer to atom ind1
                       ! -> abort
                       !
                       if (sqrt(sum((r_atoms(k,:)-r_atoms(ind1,:))**2))<
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))) then
                         counter2 = 1
                       end if
                     end if
                   end do
                   !
                   ! Only if atom j is closest to atom ind1->continue 
                   !
                   if (counter2 == 0) then
                     !
                     ! Exclude evaluation where the FOD has a larger
                     ! distance to any of the two atoms than the distance
                     ! between the two atoms
                     !
                     if ((dist1 > 
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))).or.((dist_tmp > 
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))))) then
                     else
                       ! 
                       ! Condition for X-H bonds:
                       !  1. Vectors X-FOD and H-FOD are parallel to the
                       !     bond axis
                       !  2. The direction of these vectors is exactly
                       !     parallel to the bond axis (by construction)
                       !
                       ! Vector between atoms. Pointing from initial to j
                       ! Distance between atoms.
                       ! Normalize the vector
                       !
                       tmp_vec = r_atoms(ind1,:)-r_atoms(j,:)
                       len3=sqrt(sum((r_atoms(ind1,:)-r_atoms(j,:))**2))
                       tmp_vec = tmp_vec/len3
                       !
                       ! Vector between FOD and first atom
                       ! And its length
                       ! Normalize it
                       !
                       tmp_vec2 = r_FOD(i+fod_UP,:)-r_atoms(ind1,:)
                       len2=
     &sqrt(sum((r_FOD(i+fod_UP,:)-r_atoms(ind1,:))**2))
                       tmp_vec2 = tmp_vec2/len2
                       !
                       ! Form dot porduct to bond vector between atoms
                       ! Gives the cosine of the angle between them
                       !
                       len2 = dot_product(tmp_vec,tmp_vec2)
                       !
                       ! The cosine needs to be 1, or very close to it 
                       ! Avoid numercial noise, by setting the length to 1
                       ! if is is numerically larger
                       !
                       if (abs(len2) > 1.00D0) len2 = 1.0D0
                       !
                       ! DO the same for the second atom
                       !
                       tmp_vec2 = r_FOD(i+fod_UP,:)-r_atoms(j,:)
                       len3 = 
     &sqrt(sum((r_FOD(i+fod_UP,:)-r_atoms(j,:))**2))
                       tmp_vec2 = tmp_vec2/len3
                       !
                       len3 = dot_product(tmp_vec,tmp_vec2)
                       !
                       if (abs(len3) > 1.00D0) len3 = 1.0D0
                       !
                       ! The cosines needs to be 1, or very close to it 
                       !
                       if ((abs(len2)>0.99).and.(abs(len3)>0.99)) then
                         !
                         counter = 1
                         r_internal(i,3,2) = 100.0D0
                         is_bond(i,2) = .true.
                       end if

                       if (counter == 1) then
                         !
                         ! Use force component along the atom-FOD axis.
                         ! -> use: F_i \cdot (a_i-R_i) / |a_i-R_i|
                         ! -> Store as internal coordinates
                         !
                         len1=
     &dot_product(r_FOD(i+fod_UP,:)-r_atoms(ind1,:),
     &f_FOD(i+fod_UP,:))/
     &sqrt(sum(((r_FOD(i+fod_UP,:)-r_atoms(ind1,:))**2)))
                         !
                         ! If the force along the bond is tiny -> numerical
                         ! 'noise' ->> ignore. Or we are in a homonuclear
                         ! bonding situation!
                         !
                         if (abs(len1) < 1.0D-4) then
                           f_internal(i,1,2) = 0.0D0
                         else
                           f_internal(i,1,2) = len1
                         end if
                         !
                         ! Use distance of FOD to the atom as coordinate
                         ! -> Store as internal coordinates
                         !
                         len1 = 
     &sqrt(sum(((r_FOD(i+fod_UP,:)-r_atoms(ind1,:))**2)))
                         r_internal(i,1,2) = len1
                         !
                         ! If a coordinate and a force similar to these ones
                         ! has been found before -> use that one instead !
                         ! Go through the other internal coordinates and
                         ! check. UP or DN channel!
                         !
                         ! UP Channel
                         ! If force is 0 -> don't do this.   TBD
                         !
                         do l = 1, fod_UP
                           call comp_int(r_internal(i,1,2),
     &r_internal(l,1,1),f_internal(i,1,2),f_internal(l,1,1))
                         end do                   
                         ! DN Channel
                         do l = 1, fod_DN
                           if (l.ne.i) then
                             call comp_int(r_internal(i,1,2),
     &r_internal(l,1,2),f_internal(i,1,2),f_internal(l,1,2))
                           end if
                         end do
                         !
                         ! overwrite any existing other values for this 
                         ! internal coordinate pair. Ensure that only the
                         ! ones above are used
                         !
                         r_internal(i,2,2) = 1000.0D0
                         r_internal(i,3,2) = 1000.0D0
                         f_internal(i,2,2) = 1000.0D0
                         f_internal(i,3,2) = 1000.0D0
                       end if
                     end if
                   end if

                 !
                 ! For any other bonding situation
                 !
                 else
                   !
                   ! if distance to atom j is equal to dist1 : same distance to two atoms -> bond FOD
                   !
                   if ((abs(dist_tmp-dist1)<cut_bond).and.(j/=ind1)
     &.and.(is_bond(i,2).eqv..false.))then 
                     ind2 = j
                     counter = 0
                     !
                     ! Get center between atoms
                     !
                     bond_center(:) = 
     &(r_atoms(ind1,:)+r_atoms(ind2,:))/2.0D0
                     !
                     ! Get distance of initial FOD to this bond center
                     !
                     dist2 = sqrt(sum((r_FOD(i+fod_UP,:)-
     &bond_center(:))**2))
                     !
                     ! Think about a better way to do this !
                     ! Use the center between the atoms. Analyze
                     ! distance of the FODs to that center. They should
                     ! be very similar!
                     !
                     do k = 1, fod_DN
                       !
                       ! Check distance to FOD center
                       !
                       dist_tmp=sqrt(sum((r_FOD(k+fod_UP,:)-
     &bond_center(:))**2))
                       if ((abs(dist_tmp-dist2) < 0.05D0)) then
                         !
                         ! Check distance to corresponding atom
                         !
                         dist_tmp=sqrt(sum((r_FOD(k+fod_UP,:)-
     &r_atoms(ind1,:))**2))
                         if ((abs(dist_tmp-dist1) < 0.05D0)) then
                           !
                           ! Get information about the center of all bond
                           ! FODs
                           !
                           fod_center(:)=fod_center(:)+r_FOD(k+fod_up,:)
                           !
                           ! Increase the counter for the bond FODs
                           !
                           counter = counter + 1
                           !
                           ! Set internal z-coordinates for these FODs to 100.0D0
                           ! -> differentiate them later on. Do not evaluate
                           ! these FODs later again
                           !
                           r_internal(k,3,2) = 100.0D0
                           is_bond(k,2) = .true.
                         end if
                       end if
                     end do
                     ! 
                     ! Calculate the actual bond center
                     !
                     fod_center(:) = fod_center(:)/counter
                     !
                     ! If counter = 1: Single bond. Project force onto the
                     ! bond axis. Thus, only move the FOD on the axis
                     ! Effectively: Project FOD force on FOD coordinate !
                     !
                     if (counter == 1) then
                       !
                       ! Use force component along the atom-FOD axis.
                       ! -> use: F_i \cdot (a_i-R_i) / |a_i-R_i|
                       ! -> Store as internal coordinates
                       !
                       len1=
     &dot_product(r_FOD(i+fod_UP,:)-r_atoms(ind1,:),
     &f_FOD(i+fod_UP,:))/
     &sqrt(sum(((r_FOD(i+fod_UP,:)-r_atoms(ind1,:))**2)))
                       !
                       ! If the force along the bond is tiny -> numerical
                       ! 'noise' ->> ignore. Or we are in a homonuclear
                       ! bonding situation!
                       !
                       if (abs(len1) < 2.5D-4.or.
     &(species(ind1)==species(ind2))) then
                         f_internal(i,1,2) = 0.0D0
                       else
                         f_internal(i,1,2) = len1
                       end if
                       !
                       ! Use distance of FOD to the atom as coordinate
                       ! -> Store as internal coordinates
                       !
                       len1 = 
     &sqrt(sum(((r_FOD(i+fod_UP,:)-r_atoms(ind1,:))**2)))
                       r_internal(i,1,2) = len1
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       !
                       ! UP Channel
                       do l = 1, fod_UP
                         call comp_int(r_internal(i,1,2),
     &r_internal(l,1,1),f_internal(i,1,2),f_internal(l,1,1))
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,2),
     &r_internal(l,1,2),f_internal(i,1,2),f_internal(l,1,2))
                         end if
                       end do
                       !
                       ! overwrite any existing other values for this 
                       ! internal coordinate pair. Ensure that only the
                       ! ones above are used
                       !
                       r_internal(i,2,2) = 1000.0D0
                       r_internal(i,3,2) = 1000.0D0
                       f_internal(i,2,2) = 1000.0D0
                       f_internal(i,3,2) = 1000.0D0
                     !
                     ! If counter > 1: Multiple bond. Use force
                     ! perpendicular to bond axis, along (bond-center)-FOD
                     ! axis. 
                     ! Bond-center = center of all bond-FODs!!
                     ! Maybe: Move along bond as well? Might be better
                     ! for a consistent treatment
                     !
                     else
                       !
                       ! Average forces for this structural element. 
                       ! More consistent
                       ! Analyze which FODs are in this bond
                       ! Use force component along the bond_center-FOD axis.
                       ! -> use: F_i \cdot (a_i-bond_center) / |a_i-bond_center|
                       !
                       len1 = 0.0D0
                       do k = 1, fod_DN
                         dist_tmp=
     &sqrt(sum((r_FOD(k+fod_UP,:)-bond_center(:))**2))
                         if ((abs(dist_tmp-dist2) < 0.05D0)) then
                           dist_tmp=
     &sqrt(sum((r_FOD(k+fod_UP,:)-r_atoms(ind1,:))**2))
                           if ((abs(dist_tmp-dist1) < 0.05D0)) then
                             len2 = dot_product(f_FOD(k+fod_UP,:),
     &r_FOD(k+fod_UP,:)-fod_center(:))/
     &sqrt(sum(((r_FOD(k+fod_UP,:)-fod_center(:))**2)))
                             !
                             ! add it up
                             !
                             len1 = len1 + len2
                           end if
                         end if
                       end do
                       len1 = len1/counter
                       !
                       ! Store internal force
                       !
                       f_internal(i,1,2) = len1
                       !
                       ! Store corresponding internal coordinate
                       ! (same for all points in this structural motif)
                       !
                       len1 = 
     &sqrt(sum(((r_FOD(i+fod_UP,:)-fod_center(:))**2)))
                       r_internal(i,1,2) = len1
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       ! 
                       ! UP Channel
                       do l = 1, fod_UP
                         call comp_int(r_internal(i,1,2),
     &r_internal(l,1,1),f_internal(i,1,2),f_internal(l,1,1))
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,2),
     &r_internal(l,1,2),f_internal(i,1,2),f_internal(l,1,2))
                         end if
                       end do
                       !
                       ! IF we have a homonuclear bonding situation
                       ! -> ONLY move FODs perpendicular to bond axis
                       !
                       if (species(ind1)==species(ind2)) then
                         r_internal(i,2,2) = 1000.0D0
                         f_internal(i,2,2) = 1000.0D0
                       ! 
                       !
                       ! IF heternonuclear bonding situation
                       ! -> take component along the bond axis
                       ! as well. Define one vector including
                       ! both components
                       !
                       ! Use force component along the bond_center-atom
                       ! axis.
                       ! -> use: F_i \cdot (bond_center-atom)/|bond_center-atom|
                       !
                       else
                         !
                         ! Average forces for this structural element. 
                         ! More consistent
                         ! Analyze which FODs are in this bond
                         ! Use force component along the bond_center-atom
                         ! axis.
                         ! -> use: F_i \cdot (bond_center-atom)/|bond_center-atom|
                         !
                         len2 = 0.0D0
                         do k = 1, fod_DN
                           dist_tmp=
     &sqrt(sum((r_FOD(k+fod_UP,:)-bond_center(:))**2))
                           if ((abs(dist_tmp-dist2) < 0.05D0)) then
                             dist_tmp=
     &sqrt(sum((r_FOD(k+fod_UP,:)-r_atoms(ind1,:))**2))
                             if ((abs(dist_tmp-dist1) < 0.05D0)) then
                               len3 = dot_product(f_FOD(k+fod_UP,:),
     &fod_center(:)-r_atoms(ind1,:))/
     &sqrt(sum(((fod_center(:)-r_atoms(ind1,:))**2)))
                               !
                               ! add it up
                               !
                               len2 = len2 + len3
                             end if
                           end if
                         end do
                         len2 = len2/counter
                         !
                         ! store internal coordinate
                         !
                         f_internal(i,2,2) = len2
                         !
                         ! Store corresponding internal coordinate
                         ! same for all points in this structural motif
                         !
                         len2=
     &sqrt(sum(((fod_center(:)-r_atoms(ind1,:))**2)))
                         r_internal(i,2,2) = len2
                         !
                         ! If a coordinate and a force similar to these ones
                         ! has been found before -> use that one instead !
                         ! Go through the other internal coordinates and
                         ! check. UP or DN channel!
                         !
                         ! UP Channel 
                         do l = 1, fod_UP
                           call comp_int(r_internal(i,2,2),
     &r_internal(l,2,1),f_internal(i,2,2),f_internal(l,2,1))
                         end do
                         ! DN Channel 
                         do l = 1, fod_DN
                           if (l.ne.i) then
                             call comp_int(r_internal(i,2,2),
     &r_internal(l,2,2),f_internal(i,2,2),f_internal(l,2,2))
                           end if
                         end do
                       end if
                       r_internal(i,3,2) = 1000.0D0
                       f_internal(i,3,2) = 1000.0D0
                     end if

                   end if
                 end if
               end do

               if (counter > 0) then
                 write(junk1,'(I4.4)') i
                 write(junk2,'(I4.4)') counter
                 if (counter == 1) then ! single FOD
                   write(6,777) 'There is   ',trim(junk2),' DN-bond 
     &FOD  in the structural motif of DN-bond FOD ',trim(junk1)
                 else                   ! more FODs
                   write(6,777) 'There are  ',trim(junk2),' DN-bond 
     &FODs in the structural motif of DN-bond FOD ',trim(junk1)
                 end if
               end if
             end if
           end do



           !!!!!!!!!!!!!
           ! Now, do all remaining FODs (lone/core)
           !!!!!!!!!!!!!!1
           n_motifs(:) = 0
           do i = 1, fod_DN
             dist1 = 1000.0                  ! initial value for the distance. Going to find minimum for all atoms
             do j = 1, n_atoms                ! Loop over atoms
               dist_tmp = sqrt(sum((r_FOD(i+fod_UP,:)-r_atoms(j,:))**2))
               if ((dist_tmp < dist1) .and. (dist_tmp < cutoff)) then  ! if distance to atom b is smaller than dist1 AND is is smaller than a cutoff radius
                 dist1 = dist_tmp                                      ! the cutoff radius ensures that 
                 ind1 = j
               end if
             end do
             !
             ! fix1s -> set 1s forces to 0
             !
             if (dist1 < 0.01) then
               if (constraint == 4) then
                 if (freeFOD(i+fod_UP).eqv..false.) then
                   f_FOD(i+fod_UP,:) = (/ 0.0D0, 0.0D0, 0.0D0 /)
                   write(junk1,'(I4.4)') i
                   write(6,778) 'DN-FOD no. ',trim(junk1),
     &' will be fixed (1s FOD)'
                 end if
               else
                 f_FOD(i+fod_UP,:) = (/ 0.0D0, 0.0D0, 0.0D0 /)
                 write(junk1,'(I4.4)') i
                 write(6,778) 'DN-FOD no. ',trim(junk1),
     &' will be fixed (1s FOD)'
               end if
             else
               !
               ! Check whether it is a bond or not
               ! If not -> keep going
               ! Make sure it is not a bond FOD!
               !
               if (is_bond(i,2).eqv..false.) then
                 !
                 ! Only if the fullForce contraint is suppossed to be used: Continue here
                 !
                 if ((constraint >= 2)) then
                   !
                   ! If core FOD -> symmetrize the core FODs accordingly
                   !
                   write(junk1,'(I4.4)') i
                   write(junk2,'(I4.4)') ind1
                   write(6,777) 'DN-FOD no. ',trim(junk1),
     &' is core/lone FOD in atom ',trim(junk2)
                   !
                   ! If the structural motif has already been found (distance similar to the one found before within 0.033 bohr) -> don't do everything again
                   !
                   exist = .false.
                   do l = 1, 100
                     if (abs(dist1 - r_int1(ind1,l,2)) < 0.033) then    ! any distance FOR THIS ATOM. SPIN DN
                       exist = .true.
                       exit
                     end if
                   end do
                   !
                   ! If not -> new structural element
                   !
                   if (.not.exist) then
                     n_motifs(ind1) = n_motifs(ind1) + 1
                     r_int1(ind1,n_motifs(ind1),2) = dist1    ! Store the current identifier for the structural motif.

                     r_internal(i,1:3,2) = r_FOD(i+fod_UP,1:3)
                     f_internal(i,1:3,2) = f_FOD(i+fod_UP,1:3)
                     !
                     ! Use total force on FOD_i for all other FODs
                     ! Store it as internal coordinate
                     ! Resymmetrize the other FODs after optimization
                     !
                     ! Exclude already assigned FODs
                     !
                     fod_center = (/ 0.0D0, 0.0D0, 0.0D0 /)
                     counter = 0
                     do k = 1, fod_DN
                       if (is_bond(k,2).eqv..false.) then               ! if this FOD has not been assigned yet (no bond)
                         dist2 = 
     &sqrt(sum((r_FOD(k+fod_UP,:)-r_atoms(ind1,:))**2))! check distance to the atom in question
                         if (abs(dist1 - dist2) < 0.033) then           ! If the new distance is the same as the old one 
                           counter = counter + 1
                           fod_center(:) = 
     &fod_center(:)+r_FOD(k+fod_UP,:)
                         end if
                       end if
                     end do
                     fod_center(:) = fod_center(:)/counter

                     !
                     ! write to screen
                     ! 
                     write(junk1,'(I4.4)') i
                     write(junk2,'(I4.4)') counter
                     if (counter == 1) then ! single FOD
                       write(6,777) 'There is   ',trim(junk2),
     &' core/lone DN-FOD  in the motif of DN-FOD ',trim(junk1)
                     else                   ! more FODs
                       write(6,777) 'There are  ',trim(junk2),
     &' core/lone DN-FODs in the motif of DN-FOD ',trim(junk1)
                     end if
                     !
                     ! If counter = 1: Single FOD (like single lone FOD). adjust force to
                     ! reflect symmetry of the system. Take atom-FOD
                     ! axis only and project FOD force on this axis!
                     ! Effectively: Project FOD force on FOD coordinate!
                     !
                     if (counter == 1) then
                       !
                       ! Use force component along the atom-FOD axis.
                       ! -> use: F_i \cdot (a_i-R_i) / |a_i-R_i|
                       ! -> Store as internal coordinates
                       !
                       len1=
     &dot_product(r_FOD(i+fod_UP,:)-r_atoms(ind1,:),
     &f_FOD(i+fod_UP,:))/
     &sqrt(sum(((r_FOD(i+fod_UP,:)-r_atoms(ind1,:))**2)))
                       !
                       f_internal(i,1,2) = len1
                       !
                       ! Use distance of FOD to the atom as coordinate
                       ! -> Store as internal coordinates
                       !
                       len1 =
     &sqrt(sum(((r_FOD(i+fod_UP,:)-r_atoms(ind1,:))**2)))
                       r_internal(i,1,2) = len1
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       !
                       ! UP Channel
                       do l = 1, fod_UP
                         call comp_int(r_internal(i,1,2),
     &r_internal(l,1,1),f_internal(i,1,2),f_internal(l,1,1))
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,2),
     &r_internal(l,1,2),f_internal(i,1,2),f_internal(l,1,2))
                         end if
                       end do
                       !
                       ! overwrite any existing other values for this 
                       ! internal coordinate pair. Ensure that only the
                       ! ones above are used
                       !
                       r_internal(i,2,2) = 1000.0D0
                       r_internal(i,3,2) = 1000.0D0
                       f_internal(i,2,2) = 1000.0D0
                       f_internal(i,3,2) = 1000.0D0
                     end if
                     !
                     ! If there are more than 1, but not 4
                     !
                     if ((counter > 1).and.(counter.ne.4)) then
                       !
                       ! Average forces for this structural element. 
                       ! More consistent
                       ! Analyze which FODs are in this lone region
                       ! Use force component along the bond_center-FOD axis.
                       ! -> use: F_i \cdot (a_i-bond_center) / |a_i-bond_center|
                       !
                       len1 = 0.0D0
                       do k = 1, fod_DN
                         if (is_bond(k,2).eqv..false.) then            ! if this FOD has not been assigned yet (no bond)
                           dist2=
     &sqrt(sum((r_FOD(k+fod_UP,:)-r_atoms(ind1,:))**2))                ! check distance to the atom in question
                           if (abs(dist1 - dist2) < 0.033) then        ! If the new distance is the same as the old one
                             len2 = dot_product(f_FOD(k+fod_UP,:),
     &r_FOD(k+fod_UP,:)-fod_center(:))/
     &sqrt(sum(((r_FOD(k+fod_UP,:)-fod_center(:))**2)))
                             !
                             ! add it up
                             !
                             len1 = len1 + len2
                           end if
                         end if
                       end do
                       len1 = len1/counter
                       !
                       ! Store internal force
                       !
                       f_internal(i,1,2) = len1
                       !
                       ! Store corresponding internal coordinate
                       ! (same for all points in this structural motif)
                       !
                       len1 = 
     &sqrt(sum(((r_FOD(i+fod_UP,:)-fod_center(:))**2)))
                       r_internal(i,1,2) = len1
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       ! 
                       ! UP Channel
                       do l = 1, fod_UP
                         call comp_int(r_internal(i,1,2),
     &r_internal(l,1,1),f_internal(i,1,2),f_internal(l,1,1))
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,2),
     &r_internal(l,1,2),f_internal(i,1,2),f_internal(l,1,2))
                         end if
                       end do
                       !
                       ! Average forces for this structural element. 
                       ! More consistent
                       ! Analyze which FODs are in this lone region
                       ! Use force component along the fod_center-atom
                       ! axis.
                       ! -> use: F_i \cdot (bond_center-atom)/|bond_center-atom|
                       !
                       ! If there is more than one atom
                       !
                       if (n_atoms > 1) then
                         len2 = 0.0D0
                         do k = 1, fod_DN
                           if (is_bond(k,2).eqv..false.) then          ! if this FOD has not been assigned yet (no bond)
                             dist2=
     &sqrt(sum((r_FOD(k+fod_UP,:)-r_atoms(ind1,:))**2))  ! check distance to the atom in question
                             if (abs(dist1 - dist2) < 0.033) then       ! If the new distance is the same as the old one
                               len3 = dot_product(f_FOD(k+fod_UP,:),
     &fod_center(:)-r_atoms(ind1,:))/
     &sqrt(sum(((fod_center(:)-r_atoms(ind1,:))**2)))
                               !
                               ! add it up
                               !
                               len2 = len2 + len3
                             end if
                           end if
                         end do
                         len2 = len2/counter
                         !
                         ! Store internal force
                         !
                         f_internal(i,2,2) = len2
                         !
                         ! Store corresponding internal coordinate
                         ! (same for all points in this structural motif)
                         !
                         len2=
     &sqrt(sum(((fod_center(:)-r_atoms(ind1,:))**2)))
                         r_internal(i,2,2) = len2
                         !
                         ! If a coordinate and a force similar to these ones
                         ! has been found before -> use that one instead !
                         ! Go through the other internal coordinates and
                         ! check. UP or DN channel!
                         !
                         ! UP Channel 
                         do l = 1, fod_UP
                           call comp_int(r_internal(i,2,2),
     &r_internal(l,2,1),f_internal(i,2,2),f_internal(l,2,1))
                         end do
                         ! DN Channel 
                         do l = 1, fod_DN
                           if (l.ne.i) then
                             call comp_int(r_internal(i,2,2),
     &r_internal(l,2,2),f_internal(i,2,2),f_internal(l,2,2))
                           end if
                         end do
                       !
                       ! If only one atom -> ignore this
                       !
                       else
                         r_internal(i,2,1) = 1000.0D0
                         f_internal(i,2,1) = 1000.0D0
                       end if
                       r_internal(i,3,2) = 1000.0D0
                       f_internal(i,3,2) = 1000.0D0
                     end if
                     !
                     ! If there are four -> breathing mode
                     !
                     if (counter == 4) then
                       r_internal(i,1,2) = dist1
                       r_internal(i,2,2) = 1000.0D0
                       r_internal(i,3,2) = 1000.0D0
                       f_internal(i,1,2) = 0.0D0
                       f_internal(i,2,2) = 1000.0D0
                       f_internal(i,3,2) = 1000.0D0
                       do k = 1, fod_DN
                         if ((is_bond(k,2).eqv..false.)) then
!!!     &.and.(k.ne.i))then                                                ! if this FOD has not been assigned yet (no bond)
                           dist2=
     &sqrt(sum((r_FOD(k+fod_UP,:)-r_atoms(ind1,:))**2))                ! check distance to the atom in question
                           if (abs(dist1 - dist2) < 0.033) then         ! If the new distance is the same as the old one
                             !
                             ! Take the force along the atom-FOD axis
                             !
                             len1 = dot_product(f_FOD(k+fod_UP,:),
     &r_FOD(k+fod_UP,:)-r_atoms(ind1,:))/
     &sqrt(sum(((r_FOD(k+fod_UP,:)-r_atoms(ind1,:))**2)))
                             !
                             f_internal(i,1,2) = f_internal(i,1,2)+len1
                           end if
                         end if
                       end do 
                       f_internal(i,1,2) = f_internal(i,1,2)/4.0D0
                       !
                       ! If a coordinate and a force similar to these ones
                       ! has been found before -> use that one instead !
                       ! Go through the other internal coordinates and
                       ! check. UP or DN channel!
                       !
                       ! UP Channel
                       do l = 1, fod_UP
                         call comp_int(r_internal(i,1,2),
     &r_internal(l,1,1),f_internal(i,1,2),f_internal(l,1,1))
                       end do
                       ! DN Channel
                       do l = 1, fod_DN
                         if (l.ne.i) then
                           call comp_int(r_internal(i,1,2),
     &r_internal(l,1,2),f_internal(i,1,2),f_internal(l,1,2))
                         end if
                       end do


                     end if
                   end if   ! end check whether structural element was already found
                 end if     ! end check what contraint to be used
               end if       ! end check for 1s FODs
             end if         ! end check whether it is a bond FOD
           end do         ! end 1. loop
         end if             ! constraint = 0 or not

         !
         ! Some post-processing for given constraints
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! UP CHANNEL
         ! If only fix1s -> write cooridnates and forces right here
         ! Forces on 1s are 0
         !!!!!!!!!!!!!!!!!!!!!! 
         if (constraint == 1) then
           do k = 1, fod_UP
             !
             ! Do not take the 1s into the optimization
             !
             if ((f_FOD(k,1).ne.0.0D0)
     &.and.(f_FOD(k,2).ne.0.0D0).and.(f_FOD(k,3).ne.0.0D0)) then
               r_internal(k,1:3,1) = r_FOD(k,1:3)
               f_internal(k,1:3,1) = f_FOD(k,1:3)
             end if
           end do
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! If freeFOD constraint -> restore the actual forces for the
         ! free FODs
         !!!!!!!!!!!!!!!1
         if (constraint == 4) then
           do k = 1, fod_UP
             !
             ! Do not take the 1s into the optimization
             !
             if ((freeFOD(k).eqv..true.)) then
               r_internal(k,1:3,1) = r_FOD(k,1:3)
               f_internal(k,1:3,1) = f_FOD(k,1:3)
             end if
           end do
         end if
         

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! DN CHANNEL
         ! If only fix1s -> write cooridnates and forces right here
         ! Forces on 1s are 0
         ! 
         if (constraint == 1) then
           do k = 1, fod_DN
             !
             ! Do not take the 1s into the optimization
             !
             if ((f_FOD(k+fod_UP,1).ne.0.0D0)
     &.and.(f_FOD(k+fod_UP,2).ne.0.0D0)
     &.and.(f_FOD(k+fod_UP,3).ne.0.0D0)) then
               r_internal(k,1:3,2) = r_FOD(k+fod_UP,1:3)
               f_internal(k,1:3,2) = f_FOD(k+fod_UP,1:3)
             end if
           end do
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! If freeFOD constraint -> restore the actual forces for the
         ! free FODs
         !!!!!!!!!!!!!!!1
         if (constraint == 4) then
           do k = 1, fod_DN
             !
             ! for the free FODs -> take full coordinates
             !
             if ((freeFOD(k+fod_UP).eqv..true.)) then
               r_internal(k,1:3,2) = r_FOD(k+fod_UP,1:3)
               f_internal(k,1:3,2) = f_FOD(k+fod_UP,1:3)
             end if
           end do
         end if


!
! DONE WITH FORCES
!

!
! Write fforce.dat file with updated forces (1s is zero or not)
!
         !
         ! Maybe determine f_max before this step -> better comparsion 
         ! to e.g. unconstrained optimization?
         ! Maybe not, because the constraint says it is zero. 
         !
         open(90,file='fforce.dat',status='old',action='write')
         do i = 1, fod_UP + fod_DN
           if ((constraint >= 2)) then
             write(90,fmt='(3F12.5)') f_FOD(i,1:3)
           else
             write(90,*) f_FOD(i,1:3)
           end if
         end do
         close(90)
!
! Write records file. If this file doesn't exist -> write it
!
         inquire(file='records',exist=exist)
         if (.not.exist) then
           !
           ! If it doesn't exist -> create it
           !
           open(90,file='records',status='new',action='write')
           write(90,*) ' '
           write(90,*) energy
           write(90,*) fod_UP, fod_DN
           do i = 1, fod_UP+fod_DN
             if ((constraint >= 2)) then
               write(90,fmt='(3F12.7)') r_FOD(i,1:3)
             else
               write(90,*) r_FOD(i,1:3)
             end if
           end do
           do i = 1, fod_UP+fod_DN
             if ((constraint >= 2)) then
               write(90,fmt='(3F12.5)') f_FOD(i,1:3)
             else
               write(90,*) f_FOD(i,1:3)
             end if
           end do
           close(90)
         else
           !
           ! If it exists -> write into the existing file
           !
           open(90,file='records',status='old',
     & position='append',action='write')
           write(90,*) ' '
           write(90,*) energy
           write(90,*) fod_UP, fod_DN
           do i = 1, fod_UP+fod_DN
             if ((constraint >= 2)) then
               write(90,fmt='(3F12.7)') r_FOD(i,1:3)
             else
               write(90,*) r_FOD(i,1:3)
             end if
           end do
           do i = 1, fod_UP+fod_DN
             if ((constraint >= 2)) then
               write(90,fmt='(3F12.5)') f_FOD(i,1:3)
             else
               write(90,*) f_FOD(i,1:3)
             end if
           end do
           close(90)
         end if
!
! Write fande.out. Define r_all and f_all. No need for fande.dat
!
         !
         ! Store internal coordinates in the same format as initially
         ! done, i.e. F12.7 for the coordinates and F12.5 for the forces
         ! -> Do that right here
         !
         if ((constraint >= 2)) then
           open(unit=91,file='internal',status='new',action='write')
           do i = 1, fod_UP
             if (r_internal(i,1,1).ne.1000.0D0) then
               write(91,fmt='(3F12.7)') r_internal(i,:,1)
               write(91,fmt='(3F12.5)') f_internal(i,:,1)
             end if
           end do
           do i = 1, fod_DN
             if (r_internal(i,1,2).ne.1000.0D0) then
               write(91,fmt='(3F12.7)') r_internal(i,:,2)
               write(91,fmt='(3F12.5)') f_internal(i,:,2)
             end if
           end do
           close(91)
           !
           ! Open file again and read coordinates/forces in the correct
           ! format
           !
           open(unit=91,file='internal',status='old',action='read')
           do i = 1, fod_UP
             if (r_internal(i,1,1).ne.1000.0D0) then
               read(91,fmt='(3F12.7)') r_internal(i,:,1)
               read(91,fmt='(3F12.5)') f_internal(i,:,1)
             end if
           end do
           do i = 1, fod_DN
             if (r_internal(i,1,2).ne.1000.0D0) then
               read(91,fmt='(3F12.7)') r_internal(i,:,2)
               read(91,fmt='(3F12.5)') f_internal(i,:,2)
             end if
           end do         
           close(91,status='delete')
         end if
         !
         ! Determine fmax - from internal coordinates!
         !
         fmax = 0.0D0
         do i = 1, fod_UP
           if (r_internal(i,1,1).ne.1000.0D0) then
             fmax = max(fmax,sqrt(f_internal(i,1,1)**2))
           end if
           if (r_internal(i,2,1).ne.1000.0D0) then
             fmax = max(fmax,sqrt(f_internal(i,2,1)**2))
           end if
           if ((r_internal(i,3,1).ne.1000.0D0).and.
     &(r_internal(i,3,1).ne.100.0D0)) then
             fmax = max(fmax,sqrt(f_internal(i,3,1)**2))
           end if
         end do
         do i = 1, fod_DN
           if (r_internal(i,1,2).ne.1000.0D0) then
             fmax = max(fmax,sqrt(f_internal(i,1,2)**2))
           end if
           if (r_internal(i,2,2).ne.1000.0D0) then
             fmax = max(fmax,sqrt(f_internal(i,2,2)**2))
           end if
           if ((r_internal(i,3,2).ne.1000.0D0).and.
     &(r_internal(i,3,2).ne.100.0D0)) then
             fmax = max(fmax,sqrt(f_internal(i,3,2)**2))
           end if
         end do


         ! 
         ! Idea: go through internal coordinates.
         ! If coordinates/forces are identical ->
         ! only use one for the optimizer.
         ! E.g., overwrite the other one with corresponding 
         ! Index. Then, after the optimizer is done we 
         ! can re-write the internal coordinates correctly
         ! GO through UP channel. if any other internal in the UP
         ! Channel is similar -> use only the first one.
         ! Then go through DN channel and have a look whether there
         ! is a corresponding internal coordinate in the UP channel.
         ! If so -> use the UP channel one.
         ! NOTE: One might want to go through the DN channel and check 
         ! against other DN channel FODs. might be more consistent
         !
         if ((constraint>=2)) then
           do i = 1, fod_UP-1     ! all internal coordinate besides the last one
             do j = i+1, fod_UP   ! all following internal coordinates
               do k = 1, 3        ! x, y, z
                 if ((r_internal(i,k,1).ne.1000.0D0).and.
     &(r_internal(j,k,1).ne.1000.0D0)) then
                   !
                   ! If coordinates and forces are equal -> use the first one
                   !
                   if(abs(r_internal(i,k,1)-r_internal(j,k,1))
     &<1.0D-4) then
                     if(abs(f_internal(i,k,1)-f_internal(j,k,1))
     &<1.0D-4) then
                       !
                       ! Overwrite internal forces of the second internal
                       ! coordinate for now. Write the index of the first one
                       counter2 = 0
                       ! Go through the indices of the FODs
                       do l = 1, fod_UP
                         ! if the force on the current FOD is simply the index
                         ! of another FOD -> Ignore
                         if (f_internal(j,k,1).eq.l) then
                           counter2 = 1
                         end if
                       end do
                       if (counter2 == 0) then
                         !
                         ! Do not do this is one of the FODs is free
                         !
                         if (constraint == 4) then
                           if ((freeFOD(i).eqv..true.)
     &.or.(freeFOD(j).eqv..true.)) then
                           else
                             f_internal(j,k,1) = i
                           end if
                         else
                           f_internal(j,k,1) = i
                         end if
                       end if
                     end if
                   end if
                 end if
               end do
             end do
           end do
           !
           ! Go through DN channel and do the same thing!
           !
           do i = 1, fod_UP        ! all internal coordinate for the UP channel
             do j = 1, fod_DN      ! all internal coordinate for the DN channel
               do k = 1, 3         ! x, y, z
                 if ((r_internal(i,k,1).ne.1000.0D0).and.
     &(r_internal(j,k,2).ne.1000.0D0)) then
                   !
                   ! If coordinates and forces are equal -> use the first one
                   !
                   if(abs(r_internal(i,k,1)-r_internal(j,k,2))
     &<1.0D-4) then
                     if(abs(f_internal(i,k,1)-f_internal(j,k,2))
     &<1.0D-4) then
                       !
                       ! Overwrite internal forces of the second internal
                       ! coordinate for now. Write the index of the first one
                       counter2 = 0
                       ! Go through the indices of the FODs
                       do l = 1, fod_UP
                         ! if the force on the current FOD is simply the index
                         ! of another FOD -> Ignore
                         if (f_internal(j,k,2).eq.l) then
                           counter2 = 1
                         end if
                       end do
                       if (counter2 == 0) then
                         !
                         ! Do not do this is one of the FODs is free
                         !
                         if (constraint == 4) then
                           if ((freeFOD(i).eqv..true.)
     &.or.(freeFOD(j).eqv..true.)) then
                           else
                             f_internal(j,k,2) = i
                           end if
                         else
                           f_internal(j,k,2) = i
                         end if
                       end if
                     end if
                   end if
                 end if
               end do
             end do
           end do
         end if
         !
         ! DONE
         !

         !
         ! Write r_all and f_all from the internal coordinates
         !
         ! First, allocate the arrays according to the number of
         ! internal coordinates
         !
         ! If a component is 1000.0D0 -> ignore. 
         ! Take all other components ! 
         ! With that, we can have single coordinates as input
         ! e.g. for breathing of tetrahedra
         !
         ! Take into account that some internal coordinates have been
         ! overwritten. Check for that (see last do loops)
         !
         ! To be simplified. Just take loop of 1,3

         counter = 0
         ! UP channel
         do i = 1, fod_UP
           do k = 1, 3
             if ((r_internal(i,k,1).ne.1000.0D0).and.
     &(r_internal(i,k,1).ne.100.0D0)) then
               ! If we use some constraint (other than fix1s)
               if (constraint >= 2) then
                 counter2 = 0
                 ! Go through the indices of the FODs
                 do j = 1, fod_UP
                   ! if the force on the current FOD is simply the index
                   ! of another FOD -> Ignore
                   if (f_internal(i,k,1).eq.j) then
                     counter2 = 1
                   end if
                 end do
                 if (counter2 == 0) then
                   counter = counter + 1
                 end if
               ! For fix1s or unconstrained
               ! -> always increase counter
               else
                 counter = counter + 1
               end if
             end if
           end do
         end do

         ! DN channel
         do i = 1, fod_DN
           do k = 1, 3
             if ((r_internal(i,k,2).ne.1000.0D0).and.
     &(r_internal(i,k,2).ne.100.0D0)) then
               ! If we use some constraint (other than fix1s)
               if (constraint >= 2) then
                 counter2 = 0
                 ! Go through the indices of the FODs
                 do j = 1, fod_UP
                   ! if the force on the current FOD is simply the index
                   ! of another FOD -> Ignore
                   if (f_internal(i,k,2).eq.j) then
                     counter2 = 1
                   end if
                 end do
                 if (counter2 == 0) then
                   counter = counter + 1
                 end if
               ! For fix1s or unconstrained
               ! -> always increase counter
               else
                 counter = counter + 1
               end if
             end if
           end do
         end do


!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,1,1)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1   ! x,y,z coordinates
!!!             end if
!!!           end if
!!!           if (r_internal(i,2,1).ne.1000.0D0) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,2,1)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1   ! x,y,z coordinates
!!!             end if
!!!           end if
!!!           if ((r_internal(i,3,1).ne.1000.0D0).and.
!!!     &(r_internal(i,3,1).ne.100.0D0)) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,3,1)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1   ! x,y,z coordinates
!!!             end if
!!!           end if
!!!         end do
!!!         do i = 1, fod_DN
!!!           if (r_internal(i,1,2).ne.1000.0D0) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,1,2)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1   ! x,y,z coordinates
!!!             end if
!!!           end if
!!!           if (r_internal(i,2,2).ne.1000.0D0) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,2,2)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1   ! x,y,z coordinates
!!!             end if
!!!           end if
!!!           if ((r_internal(i,3,2).ne.1000.0D0).and.
!!!     &(r_internal(i,3,2).ne.100.0D0)) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,3,2)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1   ! x,y,z coordinates
!!!             end if
!!!           end if
!!!         end do
         allocate(r_all(counter))
         allocate(f_all(counter)) 
         !
         ! Get internal coordinates for optimizer
         !
         ! Write internal coordinates (scaled or not) into the arrays
         ! which go into the optimizer
         !
         counter = 0         ! number of internal coordinates used
         count_hess = 0      ! entry of hess_D array (For scaling only)
         do i = 1, fod_UP
           do k = 1, 3
             if ((r_internal(i,k,1).ne.1000.0D0).and.
     &(r_internal(i,k,1).ne.100.0D0)) then
               ! If we use some constraint (other than fix1s)
               if (constraint >= 2) then
                 counter2 = 0
                 ! Go through the indices of the FODs
                 do j = 1, fod_UP
                   ! if the force on the current FOD is simply the index
                   ! of another FOD -> Ignore
                   if (f_internal(i,k,1).eq.j) then
                     counter2 = 1
                   end if
                 end do
                 if (counter2 == 0) then
                   counter = counter + 1
                   r_all(counter) = r_internal(i,k,1)
                   f_all(counter) = f_internal(i,k,1)
                   ! If scaled -> scale the coordinates and forces
                   if (scaled) then
                     count_hess = count_hess + 1
                     r_all(counter) = r_all(counter)/hess_D(count_hess)
                     f_all(counter) = f_all(counter)*hess_D(count_hess)
                   end if
                 ! If not taken into account -> increase counter_hess
                 ! anyway
                 else
                   count_hess = count_hess + 1
                 end if
               ! For fix1s or unconstrained
               ! -> always increase counter
               else
                 counter = counter + 1
                 r_all(counter) = r_internal(i,k,1)
                 f_all(counter) = f_internal(i,k,1)
                 ! If scaled -> scale the coordinates and forces
                 if (scaled) then
                   count_hess = count_hess + 1
                   r_all(counter) = r_all(counter)/hess_D(count_hess)
                   f_all(counter) = f_all(counter)*hess_D(count_hess)
                 end if                 
               end if
             ! If component is not used -> increase counter_hess
             else
               count_hess = count_hess + 1
             end if
           end do
         end do
                
         ! DN channel       
         do i = 1, fod_DN
           do k = 1, 3
             if ((r_internal(i,k,2).ne.1000.0D0).and.
     &(r_internal(i,k,2).ne.100.0D0)) then
               ! If we use some constraint (other than fix1s)
               if (constraint >= 2) then
                 counter2 = 0
                 ! Go through the indices of the FODs
                 do j = 1, fod_UP
                   ! if the force on the current FOD is simply the index
                   ! of another FOD -> Ignore
                   if (f_internal(i,k,2).eq.j) then
                     counter2 = 1
                   end if
                 end do
                 if (counter2 == 0) then
                   counter = counter + 1
                   r_all(counter) = r_internal(i,k,2)
                   f_all(counter) = f_internal(i,k,2)
                   ! If scaled -> scale the coordinates and forces
                   if (scaled) then
                     count_hess = count_hess + 1
                     r_all(counter) = r_all(counter)/hess_D(count_hess)
                     f_all(counter) = f_all(counter)*hess_D(count_hess)
                   end if
                 ! If not taken into account -> increase counter_hess
                 ! anyway
                 else
                   count_hess = count_hess + 1
                 end if
               ! For fix1s or unconstrained
               ! -> always increase counter
               else
                 counter = counter + 1
                 r_all(counter) = r_internal(i,k,2)
                 f_all(counter) = f_internal(i,k,2)
                 ! If scaled -> scale the coordinates and forces
                 if (scaled) then
                   count_hess = count_hess + 1
                   r_all(counter) = r_all(counter)/hess_D(count_hess)
                   f_all(counter) = f_all(counter)*hess_D(count_hess)
                 end if                 
               end if
             ! If component is not used -> increase counter_hess
             else
               count_hess = count_hess + 1
             end if
           end do
         end do



!!!           if (r_internal(i,1,1).ne.1000.0D0) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,1,1)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1
!!!               r_all(counter) = r_internal(i,1,1)  
!!!               f_all(counter) = f_internal(i,1,1)
!!!               !
!!!               ! If scaled -> scale the coordinates and forces
!!!               !
!!!               if (scaled) then
!!!                 count_hess = count_hess + 1
!!!                 r_all(counter) = r_all(counter)/hess_D(count_hess)
!!!                 f_all(counter) = f_all(counter)*hess_D(count_hess)
!!!               end if
!!!             else
!!!               count_hess = count_hess + 1
!!!             end if
!!!           !
!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!           ! to make sure we have the right entry at the right time
!!!           !
!!!           else
!!!             count_hess = count_hess + 1
!!!           end if
!!!
!!!           if (r_internal(i,2,1).ne.1000.0D0) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,2,1)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1
!!!               r_all(counter) = r_internal(i,2,1)  
!!!               f_all(counter) = f_internal(i,2,1)
!!!               !
!!!               ! If scaled -> scale the coordinates and forces
!!!               !
!!!               if (scaled) then
!!!                 count_hess = count_hess + 1
!!!                 r_all(counter) = r_all(counter)/hess_D(count_hess)
!!!                 f_all(counter) = f_all(counter)*hess_D(count_hess)
!!!               end if
!!!             else
!!!               count_hess = count_hess + 1
!!!             end if
!!!           !
!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!           ! to make sure we have the right entry at the right time
!!!           !
!!!           else
!!!             count_hess = count_hess + 1
!!!           end if
!!!
!!!           if ((r_internal(i,3,1).ne.1000.0D0).and.
!!!     &(r_internal(i,3,1).ne.100.0D0)) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,3,1)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1
!!!               r_all(counter) = r_internal(i,3,1)  
!!!               f_all(counter) = f_internal(i,3,1)
!!!               !
!!!               ! If scaled -> scale the coordinates and forces
!!!               !
!!!               if (scaled) then
!!!                 count_hess = count_hess + 1
!!!                 r_all(counter) = r_all(counter)/hess_D(count_hess)
!!!                 f_all(counter) = f_all(counter)*hess_D(count_hess)
!!!               end if
!!!             else
!!!               count_hess = count_hess + 1
!!!             end if
!!!           !
!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!           ! to make sure we have the right entry at the right time
!!!           !
!!!           else
!!!             count_hess = count_hess + 1
!!!           end if
!!!         end do
!!!         do i = 1, fod_DN
!!!           !
!!!           ! If the internal coordinate has been assigned
!!!           ! -> use it for optimization
!!!           !
!!!           if (r_internal(i,1,2).ne.1000.0D0) then           
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,1,2)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1
!!!               r_all(counter) = r_internal(i,1,2)  
!!!               f_all(counter) = f_internal(i,1,2)
!!!               !
!!!               ! If scaled -> scale
!!!               !
!!!               if (scaled) then
!!!                 count_hess = count_hess + 1
!!!                 r_all(counter) = r_all(counter)/hess_D(count_hess)
!!!                 f_all(counter) = f_all(counter)*hess_D(count_hess)
!!!               end if
!!!             else
!!!               count_hess = count_hess + 1
!!!             end if
!!!           !
!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!           ! to make sure we have the right entry at the right time
!!!           !
!!!           else
!!!             count_hess = count_hess + 1
!!!           end if
!!!
!!!           if (r_internal(i,2,2).ne.1000.0D0) then           
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,2,2)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1
!!!               r_all(counter) = r_internal(i,2,2)  
!!!               f_all(counter) = f_internal(i,2,2)
!!!               !
!!!               ! If scaled -> scale
!!!               !
!!!               if (scaled) then
!!!                 count_hess = count_hess + 1
!!!                 r_all(counter) = r_all(counter)/hess_D(count_hess)
!!!                 f_all(counter) = f_all(counter)*hess_D(count_hess)
!!!               end if
!!!             else
!!!               count_hess = count_hess + 1
!!!             end if
!!!           !
!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!           ! to make sure we have the right entry at the right time
!!!           !
!!!           else
!!!             count_hess = count_hess + 1
!!!           end if
!!!
!!!           if ((r_internal(i,3,2).ne.1000.0D0).and.
!!!     &(r_internal(i,3,2).ne.100.0D0)) then
!!!             if((((constraint >= 2))
!!!     &.and.(f_internal(i,3,2)<1.0D0)).or.
!!!     &((constraint < 2))) then
!!!               counter = counter + 1
!!!               r_all(counter) = r_internal(i,3,2)  
!!!               f_all(counter) = f_internal(i,3,2)
!!!               !
!!!               ! If scaled -> scale
!!!               !
!!!               if (scaled) then
!!!                 count_hess = count_hess + 1
!!!                 r_all(counter) = r_all(counter)/hess_D(count_hess)
!!!                 f_all(counter) = f_all(counter)*hess_D(count_hess)
!!!               end if
!!!             else
!!!               count_hess = count_hess + 1
!!!             end if
!!!           !
!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!           ! to make sure we have the right entry at the right time
!!!           !
!!!           else
!!!             count_hess = count_hess + 1
!!!           end if
!!!         end do

!
! Store internal coordinates temporarily
! Use them after optimization to place all 
! FODs per structural motifs symmetrically
! (steps after the optimization itself)
!
         allocate(r_org_int_UP(fod_UP,3))         
         allocate(r_org_int_DN(fod_DN,3))
         allocate(f_org_int_UP(fod_UP,3))         
         allocate(f_org_int_DN(fod_DN,3))
         do i = 1, fod_UP
           r_org_int_UP(i,1:3) = r_internal(i,1:3,1)
           f_org_int_UP(i,1:3) = f_internal(i,1:3,1)
           write(6,*) 'intR UP ',i,r_org_int_UP(i,1:3)
           write(6,*) 'intF UP ',i,f_org_int_UP(i,1:3)
         end do
         do i = 1, fod_DN
           r_org_int_DN(i,1:3) = r_internal(i,1:3,2)
           f_org_int_DN(i,1:3) = f_internal(i,1:3,2)
           write(6,*) 'intR DN ',i,r_org_int_DN(i,1:3)
           write(6,*) 'intF DN ',i,f_org_int_DN(i,1:3)
         end do


         !
         ! Check whether fande.out exists. Then write it accordingly
         ! 
         inquire(file='fande.out',exist=exist)
         if (.not.exist) then
         ! If it doesn't exist -> create it
           open(90,file='fande.out',status='new',action='write')
           counter = 1
           write(90,80) counter,energy,fmax !,lbfgs_step
           close(90)
         else
         ! If it does exist -> write into the existing one
         ! Open for reading
           open(90,file='fande.out',status='old',
     & position='append',action='read')
           !
           ! Find last line -> read last optimization step number 
           !
           backspace(90)       ! Go up one line -> then read
           read(90,*) counter
           close(90)
         ! Open for writing
           open(90,file='fande.out',status='old',
     & position='append',action='write')
           counter = counter + 1
           write(90,80) counter,energy,fmax !,lbfgs_step
           close(90)
         end if
 80      format(i5,f20.12,f20.12,f20.12)




         !!!! ShellOPT for structural objects. Individually !!!!!!!
         ! HERE: Choose which internal coordinate to optimize
         ! Go through internal coordinates -> Choose largest force OR 
         ! choose the structural object that has been optimized
         ! previously.  Optimize structural objects until they
         ! are converged. Then, go to the next object. 
         !
         !
         ! First, see whether there is an optimization going on. If so,
         ! Use this structural object. IF its force is larger than
         ! 10**--4
         !
         object = 0
         object_file = 0

         if (constraint == 3) then
           r_all_tmp(1) = 0.0D0
           f_all_tmp(1) = 0.0D0

           inquire(file='shellOPT',exist=exist)
           if (exist) then
             open(90,file='shellOPT',status='old',action='read')
             !
             ! Read the index of the structural motif, and the internal
             ! coordinate corresponding to it 
             !
             read(90,*) object_file, len1
             !
             ! if object_file > size(r_all) -> reset and move on
             !
             if (object_file > size(r_all)) then
               !
               ! Reset object identifier and remove shellOPT, FDIAG,
               ! FSEARCH and FSTEP
               !
               object = 0
               close(90,status='delete')
               if (object_file < 10) then
                 write(fdiag,'(A5,I1,A4)')  'FDIAG',object_file,'.LBF'
                 write(fsearch,'(A7,I1,A4)')'FSEARCH',object_file,'.LBF'
                 write(fstep,'(A5,I1,A4)')  'FSTEP',object_file,'.LBF'
               else if (object_file > 10) then
                 write(fdiag,'(A5,I2,A4)')  'FDIAG',object_file,'.LBF'
                 write(fsearch,'(A7,I2,A4)')'FSEARCH',object_file,'.LBF'
                 write(fstep,'(A5,I2,A4)')  'FSTEP',object_file,'.LBF'
               end if
               open(unit=90,file=fdiag,status='old')
               close(90,status='delete')
               open(unit=90,file=fsearch,status='old')
               close(90,status='delete')
               open(unit=90,file=fstep,status='old')
               close(90,status='delete')

             else
               !
               ! Check that the current internal coordinate at that index is
               ! the same as the one stored -> should be the same object
               !
               if (abs(r_all(object_file)-len1) < 1.0D-1) then
                 !
                 ! Check the the current force
                 !
                 if (abs(f_all(object_file)) > 1.0D-4) then
                   object = object_file
                   r_all_tmp(1) = r_all(object_file)
                   f_all_tmp(1) = f_all(object_file)
                 else
                   !
                   ! Reset object identifier and remove shellOPT, FDIAG,
                   ! FSEARCH and FSTEP
                   !
                   object = 0
                   close(90,status='delete')
                   if (object_file < 10) then
                     write(fdiag,'(A5,I1,A4)')  
     &'FDIAG',object_file,'.LBF'
                     write(fsearch,'(A7,I1,A4)')
     &'FSEARCH',object_file,'.LBF'
                     write(fstep,'(A5,I1,A4)')  
     &'FSTEP',object_file,'.LBF'
                   else if (object_file > 10) then
                     write(fdiag,'(A5,I2,A4)')  
     &'FDIAG',object_file,'.LBF'
                     write(fsearch,'(A7,I2,A4)')
     &'FSEARCH',object_file,'.LBF'
                     write(fstep,'(A5,I2,A4)')  
     &'FSTEP',object_file,'.LBF'
                   end if
                   open(unit=90,file=fdiag,status='old')
                   close(90,status='delete')
                   open(unit=90,file=fsearch,status='old')
                   close(90,status='delete')
                   open(unit=90,file=fstep,status='old')
                   close(90,status='delete')
                 end if
               !
               ! check index the coordinate has now
               ! if different from before
               ! THIS NEEDS TO BE ROBUST !!
               !
               else
                 do i = 1, size(r_all)
                   if (abs(r_all(i)-len1) < 1.0D-1) then
                     !
                     ! Check corresponding force
                     !
                     if (abs(f_all(i)) > 1.0D-4) then
                       object = i
                       r_all_tmp(1) = r_all(i)
                       f_all_tmp(1) = f_all(i)
                     else
                       !
                       ! Reset object identifier and remove shellOPT, FDIAG,
                       ! FSEARCH and FSTEP
                       !
                       object = 0
                       close(90,status='delete')
                       if (object_file < 10) then
                         write(fdiag,'(A5,I1,A4)')  
     &'FDIAG',object_file,'.LBF'
                         write(fsearch,'(A7,I1,A4)')
     &'FSEARCH',object_file,'.LBF'
                         write(fstep,'(A5,I1,A4)')  
     &'FSTEP',object_file,'.LBF'
                       else if (object_file > 10) then
                         write(fdiag,'(A5,I2,A4)')  
     &'FDIAG',object_file,'.LBF'
                         write(fsearch,'(A7,I2,A4)')
     &'FSEARCH',object_file,'.LBF'
                         write(fstep,'(A5,I2,A4)')  
     &'FSTEP',object_file,'.LBF'
                       end if
                       open(unit=90,file=fdiag,status='old')
                       close(90,status='delete')
                       open(unit=90,file=fsearch,status='old')
                       close(90,status='delete')
                       open(unit=90,file=fstep,status='old')
                       close(90,status='delete')
                     end if
                   end if
                 end do
               end if
             end if
           end if
           !
           ! If object == 0, i.e. need new optimization.
           !
           if (object == 0) then
             do i = 1, size(r_all)
               !
               ! If any force is sufficiently large
               !
               if (abs(f_all(i)) > 1.0D-4) then
                 !
                 ! Check whether it is larger than the other forces
                 !
                 if (abs(f_all(i)) > abs(f_all_tmp(1))) then
                   !
                   ! If so -> store
                   !
                   object = i
                   object_file = object
                   r_all_tmp(1) = r_all(i)
                   f_all_tmp(1) =  f_all(i)
                 end if
               end if
             end do
             ! If object = 0 -> converged. Take all internal coordinates for
             ! now and optimize them all together
             if (object == 0) then
               write(6,*) 'All forces below or equal to 1.0D-4'
               write(6,*) '-> We should stop right here (TBD)'
               write(6,*) '-> Right now: all are optimized'
               !
               ! Adjust LBFGS step size. Do NOT move that much anymore
               !
               lbfgs_step = 0.01D0
             end if
           end if
         end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now: Call the optimizer. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (do_lbfgs) then    ! if true -> do LBFGS
           ! 
           ! Initialize variables
           ! Number of internal coordinates!

           ! probably if ((object == 0).and.(constraint.ne.3)). For now,
           ! it is ok
           if (object == 0) then       ! If object == 0: either fullForce OR shellOPT with all coordinates
             n_variables = size(r_all) !  3*(intern_UP+intern_DN)   ! Number of variables
           else
             n_variables = size(r_all_tmp) !  ShellOPT
           end if
           write(6,*) 'n_variables ', n_variables
           allocate(work_space(n_variables*(2*n_corrections+1)+    ! Workspace array for the LBFGS
     &2*n_corrections))
           allocate(diag(n_variables))       ! DIAG array for the LBFGS
           work_space(:) = 0.D0
           diag(:) = 0.D0
           !
           ! Write stuff to screen
           !
           write(6,*) 'Internal coordinates and forces. BEFORE'
           do i = 1, size(r_all)
             write(6,*) r_all(i), f_all(i)
           end do
           ! For ShellOPT
           if (object.ne.0) then
             write(6,*) 'The one that will be optimized is:'
             write(6,*) r_all_tmp(1), f_all_tmp(1)
           end if


           if (object_file < 10) then
             write(fdiag,'(A5,I1,A4)')   'FDIAG',object_file,'.LBF'
             write(fsearch,'(A7,I1,A4)') 'FSEARCH',object_file,'.LBF'
             write(fstep,'(A5,I1,A4)')   'FSTEP',object_file,'.LBF'
           else if (object_file > 10) then
             write(fdiag,'(A5,I2,A4)') 'FDIAG',object_file,'.LBF'
             write(fsearch,'(A7,I2,A4)') 'FSEARCH',object_file,'.LBF'
             write(fstep,'(A5,I2,A4)')   'FSTEP',object_file,'.LBF'
           end if


           !!!!!!!!!!!!!!!!!!!!!
           ! Analyze whether LBFGS is stuck.
           ! If last three max_forces are the same ->
           ! Reset LBFGS
           ! current optimization step is stored in the variable
           ! 'counter'
           !!!!!!!!!!!!!!!!!!!!
           reset = .false.
           if (counter > 3) then
             open(90,file='fande.out',status='old',action='read')
             !
             ! skip all but the last three entries
             !
             do j = 1, counter-3
               read(90,*)
             end do
             !
             ! Read the last three entries 
             !
             read(90,*) i, len4, len1 
             read(90,*) i, len5, len2 
             read(90,*) i, len6, len3
             close(90)
             !
             ! Compare the energies and forces.
             ! If energies are within 1.0D-6, and
             ! forces are within 1.0D-4 -> LFBGS is 
             ! stuck (or it is simply converged).
             ! This way or another -> reset!
             !
             if ((abs(len6-len5)<1.0D-6)
     &.and.(abs(len5-len4)<1.0D-6).and.(abs(len6-len4)<1.0D-6)) then
               !
               ! If forces are greater than 5.0D-5
               ! - Rather abort if forces are 0?
               !
               if ((len1.ge.5.0D-5)
     &.and.(len2.ge.5.0D-5).and.(len3.ge.5.0D-5)) then
                 if (abs(len3-len2) < 1.0D-4) then
                   if (abs(len2-len1) < 1.0D-4) then
                     if (abs(len3-len1) < 1.0D-4) then
                       reset = .true.
                     end if
                   end if
                 end if
               end if
             end if
           end if      
           
           ! KAJ: check for convergence -- is the max force on an FOD less than
           ! the tolerance?
           !
           fmax = 0.0d0
           print *, 'fod_UP and fod_DN',fod_UP, fod_DN
           open(90,file='fforce.dat',form='formatted')
           do i=1,fod_UP + fod_DN
             read(90,*) (f_fod(i,j),j=1,3)
             print *, (f_fod(i,j),j=1,3)
             ftest = f_fod(i,1)**2 + f_fod(i,2)**2 + f_fod(i,3)**2
             ftest = sqrt(ftest)
             if(ftest.gt.fmax) fmax = ftest
           end do
           close(90)
           print *, 'fmax', fmax
           if(fmax.lt.0.0001) then       !set force tolerance here
             fod_converge=.true.
             print *, 'Converged in ELECTRONIC GEOMETRY deleting LBF
     &         files'
             open(90,file=fdiag,form='formatted',
     & status='unknown')
             close(90,status='delete')
             open(90,file=fsearch,form='formatted',
     & status='unknown')
             close(90,status='delete')
             open(90,file=fstep,form='formatted',
     & status='unknown')
             close(90,status='delete')
             return
           end if
           !
           ! If LBFGS is not stuck -> just proceed as usual
           !
           if (.not.reset) then
             write(6,*) 'LBFGS is not stuck'
             !
             ! Read information from previous LBF files if available
             !
             inquire(file=fdiag,exist=exist)  ! 'FDIAG.LBF'
             if(exist) then
               !
               ! Error analysis. If fewer/more internal
               ! coordinates than before -> RESET LBFGS, as
               ! number of coordinates might have increased OR
               ! decreased
               ! 
               ind1 = 0
               open(unit=90,file=fdiag,status='old',action='read')
               read(90,*) status_error
               read(90,*)
               loop321: do i = 1, n_variables
                 read(90,*,iostat=ind1) diag(i)
                 !
                 ! If there are not enough coordinates to read -> reset
                 ! LBFGS
                 !
                 if (ind1.ne.0) then
                   status_error = 0
                   write(6,*) 'Number of variables is larger than 
     &before -> Reset LBFGS'
                   exit loop321
                 end if
               end do loop321
               !
               ! If there was no error in reading the coordinates
               ! (diag), check that the next entry is NOT a number.
               ! Because if so, the current number of coordinates is
               ! smaller than what it was before. Thus, reset LBFGS
               !
               ind2 = 0
               if (ind1 == 0) then
                 !
                 ! Try to read a floating point number. If this works ->
                 ! RESET LBFGS, because we now have fewer coordinates
                 !
                 read(90,*,iostat=ind2) len1
                 if (ind2 == 0) then
                   status_error = 0
                   write(6,*) 'Number of variables is smaller than 
     &before -> Reset LBFGS'
                 !
                 ! In any other case, just read the file as usual
                 !
                 else
                   do i = 1, n_variables*(2*n_corrections+1)+
     &2*n_corrections
                     read(90,*) work_space(i)
                   end do
                 end if
               end if
               close(90)
             end if
           !
           ! If LBFGS is stuck -> reset it here
           !
           else
             write(6,*) 'LBFGS is stuck - resetting'
             status_error = 0                    ! 
           end if
           print *, 'status_error set here', status_error
           if(first) status_error = 0
           first = .false.

           !
           ! With all that information, call the LBFGS subroutine
           !
           ! Here, individual internal coordinate. ShellOPT
           !
           if (object.ne.0) then
                   print *, 'status_error', status_error
             call folbfgs(n_variables,n_corrections,
     &          r_all_tmp(1:n_variables), energy, 
     &          f_all_tmp(1:n_variables), .false., diag,
     &          additional_info, lbfgs_epsilon, lbfgs_precision,
     &          work_space, status_error, lbfgs_step, object_file)  ! last variable declares which internal coordinate to optimize. Use to write correct file names
           !
           ! Here, all internal coordinates
           !
           else
                   print *, 'what is in call to folbfgs'
c                  print *, 'n_variables n_corrections', n_variables,
c    &               n_corrections
c                  print *, 'r', (r_all(i),i=1,n_variables)
c                  print *, 'f', (f_all(i),i=1,n_variables)
c                  print *, 'diag',(diag(i),i=1,n_variables)
c                  print *, 'add info', (additional_info(i),i=1,2)
c                  print *, 'eps', lbfgs_epsilon
c                  print *, 'prec',lbfgs_precision 
                   print *, 'status_error', status_error
c                  print *, 'lbfgs_step', lbfgs_step

             call folbfgs(n_variables,n_corrections,
     &          r_all(1:n_variables), energy, 
     &          f_all(1:n_variables), .false., diag,
     &          additional_info, lbfgs_epsilon, lbfgs_precision,
     &          work_space, status_error, lbfgs_step, object_file)  ! last variable declares which internal coordinate to optimize. Here: 0 -> optimize all
           end if
          ! 
          ! After calling, r_all has the updated FOD positions
          !
           if ((status_error == 0).or.(reset)) then             ! If the LBFGS has terminated without errors or it is stuck
             open(90,file=fdiag,form='formatted',
     & status='unknown')
             close(90,status='delete')
             open(90,file=fsearch,form='formatted',
     & status='unknown')
             close(90,status='delete')
             open(90,file=fstep,form='formatted',
     & status='unknown')
             close(90,status='delete')
           end if
           if ((status_error > 0).and.(.not.reset)) then         ! If the LBFGS did not detect errors, but needs to continue
             open(unit=90,file=fdiag,status='unknown')
             write(90,*) status_error, "  IFLAG"
             write(90,*) "=====DIAG===="
             do i = 1, n_variables
              write(90,fmt='(F20.10)') diag(i)
             end do
             write(90,*) "====WORK===="
             do i = 1, n_variables*(2*n_corrections+1)+2*n_corrections
               write(90,*) work_space(i)
             end do
             close(90)
           end if 
           !
           ! Deallocate arrays
           !
           deallocate(work_space)
           deallocate(diag)

           ! KT test
           if (object == 0) then
             write(6,*) 'Internal coordinates. AFTER'
             do i = 1, size(r_all)
               write(6,*) r_all(i)
             end do
           ! KT test. ShellOPT
           else if (object.ne.0) then
             write(6,*) 'The one that was optimized is:'
             write(6,*) r_all_tmp(1)
           end if
         ! 
         ! In case you want to do CG
         ! NO SHELL-OPT YET !!!
         else
           n_variables = size(r_all)
           gtol = 0.000001d0
           ftol = 0.000001d0
           mopt = 1000 !1000 => 3000
           call fodcgrad(n_variables,mopt,energy,r_all(1:n_variables),
     &f_all(1:n_variables),gtol,ftol,scrv,istat) 
         end if
         ! 
         ! END OF OPTIMIZATION
         ! 


         !
         ! For shellOPT: use optimized internal coordinate to
         ! reconstruct everything else
         !
         ! Maybe write all optimization steps in shellOPT. Could be
         ! useful?
         !
         if (constraint == 3) then
           if (object.ne.0) then
             open(unit=90,file='shellOPT',status='unknown',
     &action='write')
             write(90,*) object_file, r_all_tmp(1)
             close(90)
             r_all(object) = r_all_tmp(1)
           end if
         end if
         !
         ! DONE 
         !




!
! RE-SYMMETRIZE FOD POSITIONS
!
! Overwrite original internal positions with optimized ones. 
! Afterwards, re-symmetrize FOD positions according to the structural
! motifs
!
         !
         ! Write optimized internal coordinates (scaled or not) 
         ! which come from the optimizer
         !
         counter = 0         ! number of internal coordinates used
         count_hess = 0      ! entry of hess_D array (For scaling only)
         ! UP Channel
         do i = 1, fod_UP
           do k = 1, 3
             !
             ! If the internal coordinate has been assigned
             ! -> use it for optimization
             !
             if ((r_org_int_UP(i,k).ne.1000.0D0).and.
     &(r_org_int_UP(i,k).ne.100.0D0)) then
               ! If we use some constraint (other than fix1s)
               if (constraint >= 2) then
                 counter2 = 0
                 ! Go through the indices of the FODs
                 do j = 1, fod_UP
                   ! if the force on the current FOD is simply the index
                   ! of another FOD -> Ignore
                   if (f_org_int_UP(i,k).eq.j) then
                     counter2 = 1
                   end if
                 end do
                 if (counter2 == 0) then
                   counter = counter + 1
                   if (scaled) then
                     count_hess = count_hess + 1
                     r_all(counter) = r_all(counter)*hess_D(count_hess)
                     f_all(counter) = f_all(counter)/hess_D(count_hess)
                   end if
                   r_internal(i,k,1) = r_all(counter)
                   f_internal(i,k,1) = f_all(counter)
                 else
                   count_hess = count_hess + 1
                 end if
               ! For fix1s or unconstrained
               ! -> always increase counter
               else
                 counter = counter + 1
                 if (scaled) then
                   count_hess = count_hess + 1
                   r_all(counter) = r_all(counter)*hess_D(count_hess)
                   f_all(counter) = f_all(counter)/hess_D(count_hess)
                 end if
                 r_internal(i,k,1) = r_all(counter)
                 f_internal(i,k,1) = f_all(counter)
               end if
             ! If component is not used -> increase counter_hess
             else
               count_hess = count_hess + 1
             end if
           end do
         end do

         ! DN Channel
         do i = 1, fod_DN
           do k = 1, 3
             !
             ! If the internal coordinate has been assigned
             ! -> use it for optimization
             !
             if ((r_org_int_DN(i,k).ne.1000.0D0).and.
     &(r_org_int_DN(i,k).ne.100.0D0)) then
               ! If we use some constraint (other than fix1s)
               if (constraint >= 2) then
                 counter2 = 0
                 ! Go through the indices of the FODs
                 do j = 1, fod_UP
                   ! if the force on the current FOD is simply the index
                   ! of another FOD -> Ignore
                   if (f_org_int_DN(i,k).eq.j) then
                     counter2 = 1
                   end if
                 end do
                 if (counter2 == 0) then
                   counter = counter + 1
                   if (scaled) then
                     count_hess = count_hess + 1
                     r_all(counter) = r_all(counter)*hess_D(count_hess)
                     f_all(counter) = f_all(counter)/hess_D(count_hess)
                   end if
                   r_internal(i,k,2) = r_all(counter)
                   f_internal(i,k,2) = f_all(counter)
                 else
                   count_hess = count_hess + 1
                 end if
               ! For fix1s or unconstrained
               ! -> always increase counter
               else
                 counter = counter + 1
                 if (scaled) then
                   count_hess = count_hess + 1
                   r_all(counter) = r_all(counter)*hess_D(count_hess)
                   f_all(counter) = f_all(counter)/hess_D(count_hess)
                 end if
                 r_internal(i,k,2) = r_all(counter)
                 f_internal(i,k,2) = f_all(counter)
               end if
             ! If component is not used -> increase counter_hess
             else
               count_hess = count_hess + 1
             end if
           end do
         end do


!!!!!           if (r_org_int_UP(i,1).ne.1000.0D0) then
!!!!!             if((((constraint >= 2))
!!!!!     &.and.(f_org_int_UP(i,1)<1.0D0)).or.
!!!!!     &((constraint < 2))) then
!!!!!               !
!!!!!               ! If scaled -> scale the coordinates and forces
!!!!!               !
!!!!!               counter = counter + 1
!!!!!               if (scaled) then
!!!!!                 count_hess = count_hess + 1
!!!!!                 r_all(counter) = r_all(counter)*hess_D(count_hess)
!!!!!                 f_all(counter) = f_all(counter)/hess_D(count_hess)
!!!!!               end if
!!!!!               r_internal(i,1,1) = r_all(counter)
!!!!!               f_internal(i,1,1) = f_all(counter)
!!!!!             else
!!!!!               count_hess = count_hess + 1
!!!!!             end if
!!!!!           !
!!!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!!!           ! to make sure we have the right entry at the right time
!!!!!           !
!!!!!           else
!!!!!             count_hess = count_hess + 1
!!!!!           end if
!!!!!
!!!!!           if (r_org_int_UP(i,2).ne.1000.0D0) then
!!!!!             if((((constraint >= 2))
!!!!!     &.and.(f_org_int_UP(i,2)<1.0D0)).or.
!!!!!     &((constraint < 2))) then
!!!!!               !
!!!!!               ! If scaled -> scale the coordinates and forces
!!!!!               !
!!!!!               counter = counter + 1
!!!!!               if (scaled) then
!!!!!                 count_hess = count_hess + 1
!!!!!                 r_all(counter) = r_all(counter)*hess_D(count_hess)
!!!!!                 f_all(counter) = f_all(counter)/hess_D(count_hess)
!!!!!               end if
!!!!!               r_internal(i,2,1) = r_all(counter)
!!!!!               f_internal(i,2,1) = f_all(counter)
!!!!!             else
!!!!!               count_hess = count_hess + 1
!!!!!             end if
!!!!!           !
!!!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!!!           ! to make sure we have the right entry at the right time
!!!!!           !
!!!!!           else
!!!!!             count_hess = count_hess + 1
!!!!!           end if
!!!!!
!!!!!           if ((r_org_int_UP(i,3).ne.1000.0D0).and.
!!!!!     &(r_org_int_UP(i,3).ne.100.0D0)) then
!!!!!             if((((constraint >= 2))
!!!!!     &.and.(f_org_int_UP(i,3)<1.0D0)).or.
!!!!!     &((constraint < 2))) then
!!!!!               !
!!!!!               ! If scaled -> scale the coordinates and forces
!!!!!               !
!!!!!               counter = counter + 1
!!!!!               if (scaled) then
!!!!!                 count_hess = count_hess + 1
!!!!!                 r_all(counter) = r_all(counter)*hess_D(count_hess)
!!!!!                 f_all(counter) = f_all(counter)/hess_D(count_hess)
!!!!!               end if
!!!!!               r_internal(i,3,1) = r_all(counter)
!!!!!               f_internal(i,3,1) = f_all(counter)
!!!!!             else
!!!!!               count_hess = count_hess + 1
!!!!!             end if
!!!!!           !
!!!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!!!           ! to make sure we have the right entry at the right time
!!!!!           !
!!!!!           else
!!!!!             count_hess = count_hess + 1
!!!!!           end if
!!!!!         end do
!!!!!
!!!!!         do i = 1, fod_DN
!!!!!           !
!!!!!           ! If the internal coordinate has been assigned
!!!!!           ! -> use it for optimization
!!!!!           !
!!!!!           if (r_org_int_DN(i,1).ne.1000.0D0) then
!!!!!             if((((constraint >= 2))
!!!!!     &.and.(f_org_int_DN(i,1)<1.0D0)).or.
!!!!!     &((constraint < 2))) then
!!!!!               !
!!!!!               ! If scaled -> scale the coordinates and forces
!!!!!               !
!!!!!               counter = counter + 1
!!!!!               if (scaled) then
!!!!!                 count_hess = count_hess + 1
!!!!!                 r_all(counter) = r_all(counter)*hess_D(count_hess)
!!!!!                 f_all(counter) = f_all(counter)/hess_D(count_hess)
!!!!!               end if
!!!!!               r_internal(i,1,2) = r_all(counter)
!!!!!               f_internal(i,1,2) = f_all(counter)
!!!!!             else
!!!!!               count_hess = count_hess + 1
!!!!!             end if
!!!!!           !
!!!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!!!           ! to make sure we have the right entry at the right time
!!!!!           !
!!!!!           else
!!!!!             count_hess = count_hess + 1
!!!!!           end if
!!!!!
!!!!!           if (r_org_int_DN(i,2).ne.1000.0D0) then
!!!!!             if((((constraint >= 2))
!!!!!     &.and.(f_org_int_DN(i,2)<1.0D0)).or.
!!!!!     &((constraint < 2))) then
!!!!!               !
!!!!!               ! If scaled -> scale the coordinates and forces
!!!!!               !
!!!!!               counter = counter + 1
!!!!!               if (scaled) then
!!!!!                 count_hess = count_hess + 1
!!!!!                 r_all(counter) = r_all(counter)*hess_D(count_hess)
!!!!!                 f_all(counter) = f_all(counter)/hess_D(count_hess)
!!!!!               end if
!!!!!               r_internal(i,2,2) = r_all(counter)
!!!!!               f_internal(i,2,2) = f_all(counter)
!!!!!             else
!!!!!               count_hess = count_hess + 1
!!!!!             end if
!!!!!           !
!!!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!!!           ! to make sure we have the right entry at the right time
!!!!!           !
!!!!!           else
!!!!!             count_hess = count_hess + 1
!!!!!           end if
!!!!!
!!!!!           if ((r_org_int_DN(i,3).ne.1000.0D0).and.
!!!!!     &(r_org_int_DN(i,3).ne.100.0D0)) then
!!!!!             if((((constraint >= 2))
!!!!!     &.and.(f_org_int_DN(i,3)<1.0D0)).or.
!!!!!     &((constraint < 2))) then
!!!!!               !
!!!!!               ! If scaled -> scale the coordinates and forces
!!!!!               !
!!!!!               counter = counter + 1
!!!!!               if (scaled) then
!!!!!                 count_hess = count_hess + 1
!!!!!                 r_all(counter) = r_all(counter)*hess_D(count_hess)
!!!!!                 f_all(counter) = f_all(counter)/hess_D(count_hess)
!!!!!               end if
!!!!!               r_internal(i,3,2) = r_all(counter)
!!!!!               f_internal(i,3,2) = f_all(counter)
!!!!!             else
!!!!!               count_hess = count_hess + 1
!!!!!             end if
!!!!!           !
!!!!!           ! If no internal coordinate is there -> add 3 to count_hess
!!!!!           ! to make sure we have the right entry at the right time
!!!!!           !
!!!!!           else
!!!!!             count_hess = count_hess + 1
!!!!!           end if
!!!!!         end do



         !!!!!!!!!!!!!! 
         ! Get back all internal coordinates from the optimized
         ! one. Use the same logic as before, just overwrite the new
         ! internal coordinates. Use the original internal coordinates
         !!!!!!!!!!!!!
         if ((constraint >= 2)) then
           do i = 1, fod_UP-1     ! all internal coordinate besides the last one
             do j = i+1, fod_UP   ! all following internal coordinates
               do k = 1, 3        ! x, y, z
                 !
                 ! If index is found -> write new internal coords
                 !
                 if (f_org_int_UP(j,k) == i) then
                   r_internal(j,k,1) = r_internal(i,k,1)
                   f_internal(j,k,1) = f_internal(i,k,1)
                 end if
               end do
             end do
           end do
           !
           ! Go through DN channel and do the same thing!
           !
           do i = 1, fod_UP        ! all internal coordinate for the UP channel
             do j = 1, fod_DN      ! all internal coordinate for the DN channel
               do k = 1, 3         ! x, y, z
                 !
                 ! If index is found -> write new internal coords
                 !
                 if (f_org_int_DN(j,k) == i) then
                   r_internal(j,k,2) = r_internal(i,k,1)
                   f_internal(j,k,2) = f_internal(i,k,1)
                 end if
               end do
             end do
           end do
         end if
         !
         ! DONE
         !







!
! separate things for unconstrained, fix1s and fullForce
!
         !
         ! If unconstrained -> just take all coordinates and forces
         !
         if (constraint == 0) then
           do i = 1, fod_UP
             r_FOD(i,1:3) = r_internal(i,1:3,1)
!             f_FOD(i,1:3) = f_internal(i,1:3,1)
           end do
           do i = 1, fod_DN
             r_FOD(i+fod_UP,1:3) = r_internal(i,1:3,2)
!             f_FOD(i+fod_UP,1:3) = f_internal(i,1:3,2)
           end do
         !
         ! In any other case (fix1s, fullForce)
         !
         else
!
! Symmetrically place new FOD position after optimization. Take first FOD position per structural element and rotate with respect
! to the initial positions (angles between them)
! Rotation matrix between old and new position of FOD i -> use this to place the other FODs
!
           deallocate(is_bond)
           !
           ! Use this array again
           !
           allocate(is_bond(fod_UP+fod_DN,2))             ! initialize size, for each spin (1,2)
           is_bond(:,:) = .false.
!
! r_tmp = original FOD positions
! r_FOD will be overwritten
!
           r_int1(:,:,:)   = 100.0D0  ! reset this array
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! 1. loop over FODs - UP channel !
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !
           ! determine which atom the FOD is closest to
           ! 
           ! Bond FODs first
           !
           do i = 1, fod_UP
             dist1 = 1000.0                                            ! initial value for the distance. Going to find minimum for all atoms
             do j = 1, n_atoms                                          ! Loop over atoms
               dist_tmp = sqrt(sum((r_tmp(i,:)-r_atoms(j,:))**2))
               if ((dist_tmp < dist1) .and. (dist_tmp < cutoff)) then  ! if distance to atom b is smaller than dist1 AND is is smaller than a cutoff radius
                 dist1 = dist_tmp                                      ! the cutoff radius ensures that we are only taking FODs for one atom 
                 ind1 = j
               end if
             end do
             !
             ! For bond FODs -> get FOD information
             ! Count how many bond FODs there are
             ind2 = 10000
             counter = 0                     ! count the number of bond FODs
             if (n_atoms > 1 .and. r_internal(i,3,1).ne.100.0D0) then
               bond_center(:) =(/ 0.0D0, 0.0D0, 0.0D0 /)
               fod_center(:) =(/ 0.0D0, 0.0D0, 0.0D0 /)
               do j = 1, n_atoms                ! Loop over atoms
                 dist_tmp = sqrt(sum((r_tmp(i,:)-r_atoms(j,:))**2))
                 !
                 ! For X-H bonds: see whether there are any H in the
                 ! bond. If so -> single bond FODs!
                 !
                 ! If either is H
                 ! atom j is NOT atom ind1
                 ! FOD i it is not already a bond
                 ! FOD is not a 1s for either of the atoms
                 !
                 if (((trim(adjustl(species(ind1)))=='H').or.
     &(trim(adjustl(species(j)))=='H'))
     &.and.(j/=ind1).and.(is_bond(i,1).eqv..false.)
     &.and.(dist1>0.01).and.(dist_tmp>0.01)) then
                   !
                   ! Check whether atom j is the closest atom to atom
                   ! ind1. If so, continue here. If not -> do not do
                   ! anything
                   !
                   counter2 = 0
                   do k = 1, n_atoms
                     if ((k.ne.j).and.(k.ne.ind1)) then
                       !
                       ! If any atom is closer to atom ind1
                       ! -> abort
                       !
                       if (sqrt(sum((r_atoms(k,:)-r_atoms(ind1,:))**2))<
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))) then
                         counter2 = 1
                       end if
                     end if
                   end do
                   !
                   ! Only if atom j is closest to atom ind1->continue 
                   !
                   if (counter2 == 0) then
                     !
                     ! Exclude evaluation where the FOD has a larger
                     ! distance to any of the two atoms than the distance
                     ! between the two atoms
                     !
                     if ((dist1 > 
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))).or.((dist_tmp > 
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))))) then
                     else
                       !
                       ! get correct new internal coordinate
                       !
                       len1 = r_internal(i,1,1)

                       ! 
                       ! Condition for X-H bonds:
                       !  1. Vectors X-FOD and H-FOD are parallel to the
                       !     bond axis
                       !  2. The direction of these vectors is exactly
                       !     parallel to the bond axis (by construction)
                       !
                       ! Vector between atoms. Pointing from initial to j
                       ! Distance between atoms.
                       ! Normalize the vector
                       !
                       tmp_vec = r_atoms(ind1,:)-r_atoms(j,:)
                       len3=sqrt(sum((r_atoms(ind1,:)-r_atoms(j,:))**2))
                       tmp_vec = tmp_vec/len3
                       !
                       ! Vector between FOD and first atom
                       ! And its length
                       ! Normalize it
                       !
                       tmp_vec2 = r_tmp(i,:)-r_atoms(ind1,:)
                       len2 = sqrt(sum((r_tmp(i,:)-r_atoms(ind1,:))**2))
                       tmp_vec2 = tmp_vec2/len2
                       !
                       ! Form dot porduct to bond vector between atoms
                       ! Gives the cosine of the angle between them
                       !
                       len2 = dot_product(tmp_vec,tmp_vec2)
                       !
                       ! The cosine needs to be 1, or very close to it 
                       ! Avoid numercial noise, by setting the length to 1
                       ! if is is numerically larger
                       !
                       if (abs(len2) > 1.00D0) len2 = 1.0D0
                       !
                       ! DO the same for the second atom
                       !
                       tmp_vec2 = r_tmp(i,:)-r_atoms(j,:)
                       len3 = sqrt(sum((r_tmp(i,:)-r_atoms(j,:))**2))
                       tmp_vec2 = tmp_vec2/len3
                       !
                       len3 = dot_product(tmp_vec,tmp_vec2)
                       !
                       if (abs(len3) > 1.00D0) len3 = 1.0D0
                       !
                       ! The cosines needs to be 1, or very close to it 
                       !
                       if ((abs(len2)>0.99).and.(abs(len3)>0.99)) then
                         is_bond(i,1) = .true.
                         r_FOD(i,:) = (r_tmp(i,:)-r_atoms(ind1,:))/
     &(sqrt(sum((r_atoms(ind1,:)-r_tmp(i,:))**2)))*len1+r_atoms(ind1,:)
                       end if
                     end if
                   end if
                 !
                 ! For any other bonding situation
                 !
                 else
                   if ((abs(dist_tmp-dist1)<cut_bond).and.(j/=ind1)
     &.and.(is_bond(i,1).eqv..false.)) then           ! if distance to atom j is equal to dist1 : same distance to two atoms -> bond FOD
                     ind2 = j
                     counter = 0
                     !
                     ! Get center between atoms
                     !
                     bond_center(:) = 
     &(r_atoms(ind1,:)+r_atoms(ind2,:))/2.0D0
                     !
                     ! Get distance of initial FOD to this bond center
                     !
                     dist2 = sqrt(sum((r_tmp(i,:)-bond_center(:))**2))
                     !
                     ! Think about a better way to do this !
                     ! Use the center between the atoms. Analyze
                     ! distance of the FODs to that center. They should
                     ! be very similar!
                     !
                     do k = 1, fod_UP
                       !
                       ! Check distance to FOD center
                       !
                       dist_tmp=
     &sqrt(sum((r_tmp(k,:)-bond_center(:))**2))
                       if ((abs(dist_tmp-dist2) < 0.05D0)) then
                         !
                         ! Check distance to corresponding atom
                         !
                         dist_tmp=
     &sqrt(sum((r_tmp(k,:)-r_atoms(ind1,:))**2))
                         if ((abs(dist_tmp-dist1) < 0.05D0)) then
                           !
                           ! Get information about the center of all bond
                           ! FODs
                           !
                           fod_center(:) = fod_center(:)+r_tmp(k,:)
                           !
                           ! Increase the counter for the bond FODs
                           !
                           counter = counter + 1
                           !
                           is_bond(k,1) = .true.
                         end if
                       end if
                     end do
                     ! 
                     ! Calculate the actual bond center
                     !
                     fod_center(:) = fod_center(:)/counter
                     !
                     ! get correct new internal coordinate
                     !
                     len1 = r_internal(i,1,1)
                     len2 = r_internal(i,2,1)
                     !
                     ! If counter = 1: Single bond. Use new internal
                     ! coordinate to wriute new FOD
                     !
                     if (counter == 1) then
                       r_FOD(i,:) = (r_tmp(i,:)-r_atoms(ind1,:))/
     &(sqrt(sum((r_tmp(i,:)-r_atoms(ind1,:))**2)))*len1+r_atoms(ind1,:)
      
                     !!!!!!!!!
                     ! For other bond situations
                     !!!!!!
                     else
                       !
                       ! Now, get the right coordinates for all FODs
                       !
                       do k = 1, fod_UP
                         !
                         ! Check distance to FOD center
                         !
                         dist_tmp=
     &sqrt(sum((r_tmp(k,:)-bond_center(:))**2))
                         if ((abs(dist_tmp-dist2) < 0.05D0)) then
                           !
                           ! Check distance to corresponding atom
                           !
                           dist_tmp=
     &sqrt(sum((r_tmp(k,:)-r_atoms(ind1,:))**2))
                           if ((abs(dist_tmp-dist1) < 0.05D0)) then
                             !
                             ! For now: Use ratio between new internal
                             ! coordinate and old coordinate .. Do some more
                             !
                             ! Define ratio
                             !
                             !
                             ! For homonuclear bonding situations!
                             !
                             if (species(ind1)==species(ind2)) then
                               r_FOD(k,:) = (r_tmp(k,:)-fod_center(:))/
     &(sqrt(sum((r_tmp(k,:)-fod_center(:))**2)))*len1 + fod_center(:)
                             !
                             ! For heteronuclear bonding situations
                             ! Need new_FOD-bond_center vector.
                             ! Rotate this vector into the other FODs.
                             ! Rotate !around! the bond axis
                             !
                             else
                               !
                               ! get new coordinate
                               ! 
                               r_FOD(k,:) = (r_tmp(k,:)-fod_center(:))/
     &(sqrt(sum((r_tmp(k,:)-fod_center(:))**2)))*len1 !! + fod_center(:)


                               r_FOD(k,:) = r_FOD(k,:) + 
     &(fod_center(:)-r_atoms(ind1,:))/
     &(sqrt(sum((fod_center(:)-r_atoms(ind1,:))**2)))*len2 +
     &r_atoms(ind1,:)


                             end if
                           end if
                         end if
                       end do
                     end if
                   end if
                 end if
               end do
             end if
           end do



           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Now, do all remaining FODs
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           n_motifs(:) = 0
           do i = 1, fod_UP
             dist1 = 1000.0                                            ! initial value for the distance. Going to find minimum for all atoms
             do j = 1, n_atoms                                          ! Loop over atoms
               dist_tmp = sqrt(sum((r_tmp(i,:)-r_atoms(j,:))**2))
               if ((dist_tmp < dist1) .and. (dist_tmp < cutoff)) then  ! if distance to atom b is smaller than dist1 AND is is smaller than a cutoff radius
                 dist1 = dist_tmp                                      ! the cutoff radius ensures that we are only taking FODs for one atom 
                 ind1 = j
               end if
             end do
             !
             ! If it is not a bond FOD
             !
             if (is_bond(i,1).eqv..false.) then
               !
               ! Only if the fullForce contraint is suppossed to be used: Continue here
               !
               if ((constraint >= 2)) then     
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                 ! Use center of lone FODs as bond_center -> keep lone
                 ! FODs on planes if they start on planes/lines (3 or 2
                 ! FODs). Tetrahedra will be conserved like this as well
                 !
                 fod_center = (/ 0.0D0, 0.0D0, 0.0D0 /)
                 !
                 ! If the structural motif has already been found (distance similar to the one found before within 0.033 bohr) -> don't do everything again
                 !
                 exist = .false.
                 do l = 1, 100
                   if (abs(dist1 - r_int1(ind1,l,1)) < 0.033) then    ! any distance FOR THIS ATOM. SPIN UP
                     exist = .true.
                     exit
                   end if
                 end do
                 !
                 ! If not -> new structural element
                 !
                 if ((.not.exist).and.
     &(r_internal(i,1,1) .ne. 1000.0D0)) then
                   counter = 0       ! count FODs in this structural element -> define 'bond_center'
                   do k = 1, fod_UP
                     if (is_bond(k,1).eqv..false.) then              ! do not do bond FODs
                       dist2 = 
     &sqrt(sum((r_tmp(k,:)-r_atoms(ind1,:))**2))                     ! check distance to the atom in question
                       if (abs(dist1 - dist2) < 0.033) then            ! If the new distance is the same as the old one               
                         counter = counter + 1
                         fod_center(:) = fod_center(:) + r_tmp(k,:)  ! bond center of the original positions!
                       end if
                     end if
                   end do
                   fod_center(:) = fod_center(:)/counter
                 end if
               end if


               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Resymmetrize lone/core FODs
               !!!!
               ! Only if the fullForce contraint is suppossed to be used: Continue here
               !
               if ((constraint >= 2)) then     
                 !
                 ! If the structural motif has already been found (distance similar to the one found before within 0.033 bohr) -> don't do everything again
                 !
                 exist = .false.
                 do l = 1, 100
                   if (abs(dist1 - r_int1(ind1,l,1)) < 0.033) then    ! any distance FOR THIS ATOM. SPIN UP
                     exist = .true.
                     exit
                   end if
                 end do
                 !
                 ! If not -> new structural element
                 !
                 if ((.not.exist).and. 
     &(r_internal(i,1,1) .ne. 1000.0D0)) then
                   n_motifs(ind1) = n_motifs(ind1) + 1
                   r_int1(ind1,n_motifs(ind1),1) = dist1    ! Store the current identifier for the structural motif.
                   !
                   ! Use internal FOD position for all other FODs
                   ! NEED: rotation matrix between old internal r
                   ! (r_org_UP/r_org_DN) and new internal r 
                   !       -> use this properly to rotate new r_i towards new r_j
                   !
                   ! constrain symmetry to structural motifs
                   !
                   !
                   ! Idea: use center of the FODs to re-symmetrize, such
                   ! like for bond FODs? might give a better picture
                   !
                   len1 = r_internal(i,1,1)
                   len2 = r_internal(i,2,1)
                   ! 
                   ! if counter == 1 -> just use internal
                   ! cooridnate. Rescale old position to new size
                   !
                   if (counter==1) then
                     r_FOD(i,:) = (r_tmp(i,:)-r_atoms(ind1,:))/
     &(sqrt(sum((r_tmp(i,:)-r_atoms(ind1,:))**2)))*len1+r_atoms(ind1,:)
                   !
                   ! If NOT a tetrahedron
                   !
                   else if (counter.ne.4) then
                     !
                     ! get new coords
                     !
                     do k = 1, fod_UP
                       if (is_bond(k,1).eqv..false.) then              ! do not do bond FODs
                         dist2 = 
     &sqrt(sum((r_tmp(k,:)-r_atoms(ind1,:))**2))                       ! check distance to the atom in question
                         if (abs(dist1 - dist2) < 0.033) then          ! If the new distance is the same as the old one
                           r_FOD(k,:) = (r_tmp(k,:)-fod_center(:))/
     &(sqrt(sum((r_tmp(k,:)-fod_center(:))**2)))*len1 !! + fod_center(:)
                           !
                           ! If more than one atom
                           !
                           if (n_atoms > 1) then
                             r_FOD(k,:) = r_FOD(k,:) +
     &(fod_center(:)-r_atoms(ind1,:))/
     &(sqrt(sum((fod_center(:)-r_atoms(ind1,:))**2)))*len2 +
     &r_atoms(ind1,:)
                           end if
                         end if
                       end if
                     end do

                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ! For four points (tetrahedron)
                   ! Do breathing mode
                   ! 
                   else
                     len1 = r_internal(i,1,1)        ! New internal FOD position of FOD i
                     !
                     ! Apply new distance to all FODs
                     !
                     do k = 1, fod_UP                    ! Find the other internal coordinate
                       if (is_bond(k,1).eqv..false.) then
!!!     &.and.(k.ne.i))then                                              ! if this FOD has not been assigned yet (no bond)
                         dist2=
     &sqrt(sum((r_atoms(ind1,:)-r_tmp(k,:))**2))                ! check distance to the atom in question
                         if (abs(dist1 - dist2) < 0.033) then   ! If the new distance is the same as the old one
                           r_FOD(k,:) =
     &(r_tmp(k,:)-r_atoms(ind1,:))/dist2*len1 + r_atoms(ind1,:)
                         end if
                       end if
                     end do

                   end if
                 end if   ! end check whether structural element was already found
               end if       ! check which contraint to use
             end if         ! check that it isn't a bond FOD
           end do         ! end 1. loop

!
! DO DN-CHANNEL!
!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! 1. loop over FODs - DN channel !
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !
           ! determine which atom the FOD is closest to
           ! 
           ! Bond FODs first
           !
           n_motifs(:) = 0
           do i = 1, fod_DN
             dist1 = 1000.0                                            ! initial value for the distance. Going to find minimum for all atoms
             do j = 1, n_atoms                                          ! Loop over atoms
               dist_tmp = sqrt(sum((r_tmp(i+fod_UP,:)-r_atoms(j,:))**2))
               if ((dist_tmp < dist1) .and. (dist_tmp < cutoff)) then  ! if distance to atom b is smaller than dist1 AND is is smaller than a cutoff radius
                 dist1 = dist_tmp                                      ! the cutoff radius ensures that we are only taking FODs for one atom 
                 ind1 = j
               end if
             end do
             !
             ! For bond FODs -> get FOD information
             ! Count how many bond FODs there are
             ind2 = 10000
             counter = 0                     ! count the number of bond FODs
             if (n_atoms > 1 .and. r_internal(i,3,2).ne.100.0D0) then
               bond_center(:) =(/ 0.0D0, 0.0D0, 0.0D0 /)
               fod_center(:) =(/ 0.0D0, 0.0D0, 0.0D0 /)
               do j = 1, n_atoms                ! Loop over atoms
                 dist_tmp=sqrt(sum((r_tmp(i+fod_UP,:)-r_atoms(j,:))**2))
                 !
                 ! For X-H bonds: see whether there are any H in the
                 ! bond. If so -> single bond FODs!
                 !
                 ! If either is H
                 ! atom j is NOT atom ind1
                 ! FOD i it is not already a bond
                 ! FOD is not a 1s for either of the atoms
                 !
                 if (((trim(adjustl(species(ind1)))=='H').or.
     &(trim(adjustl(species(j)))=='H'))
     &.and.(j/=ind1).and.(is_bond(i,2).eqv..false.)
     &.and.(dist1>0.01).and.(dist_tmp>0.01)) then
                   !
                   ! Check whether atom j is the closest atom to atom
                   ! ind1. If so, continue here. If not -> do not do
                   ! anything
                   !
                   counter2 = 0
                   do k = 1, n_atoms
                     if ((k.ne.j).and.(k.ne.ind1)) then
                       !
                       ! If any atom is closer to atom ind1
                       ! -> abort
                       !
                       if (sqrt(sum((r_atoms(k,:)-r_atoms(ind1,:))**2))<
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))) then
                         counter2 = 1
                       end if
                     end if
                   end do
                   !
                   ! Only if atom j is closest to atom ind1->continue 
                   !
                   if (counter2 == 0) then
                     !
                     ! Exclude evaluation where the FOD has a larger
                     ! distance to any of the two atoms than the distance
                     ! between the two atoms
                     !
                     if ((dist1 > 
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))).or.((dist_tmp > 
     &sqrt(sum((r_atoms(j,:)-r_atoms(ind1,:))**2))))) then
                     else
                       !
                       ! get correct new internal coordinate
                       !
                       len1 = r_internal(i,1,2)

                       ! 
                       ! Condition for X-H bonds:
                       !  1. Vectors X-FOD and H-FOD are parallel to the
                       !     bond axis
                       !  2. The direction of these vectors is exactly
                       !     parallel to the bond axis (by construction)
                       !
                       ! Vector between atoms. Pointing from initial to j
                       ! Distance between atoms.
                       ! Normalize the vector
                       !
                       tmp_vec = r_atoms(ind1,:)-r_atoms(j,:)
                       len3=sqrt(sum((r_atoms(ind1,:)-r_atoms(j,:))**2))
                       tmp_vec = tmp_vec/len3
                       !
                       ! Vector between FOD and first atom
                       ! And its length
                       ! Normalize it
                       !
                       tmp_vec2 = r_tmp(i+fod_UP,:)-r_atoms(ind1,:)
                       len2=
     &sqrt(sum((r_tmp(i+fod_UP,:)-r_atoms(ind1,:))**2))
                       tmp_vec2 = tmp_vec2/len2
                       !
                       ! Form dot porduct to bond vector between atoms
                       ! Gives the cosine of the angle between them
                       !
                       len2 = dot_product(tmp_vec,tmp_vec2)
                       !
                       ! The cosine needs to be 1, or very close to it 
                       ! Avoid numercial noise, by setting the length to 1
                       ! if is is numerically larger
                       !
                       if (abs(len2) > 1.00D0) len2 = 1.0D0
                       !
                       ! DO the same for the second atom
                       !
                       tmp_vec2 = r_tmp(i+fod_UP,:)-r_atoms(j,:)
                       len3 = 
     &sqrt(sum((r_tmp(i+fod_UP,:)-r_atoms(j,:))**2))
                       tmp_vec2 = tmp_vec2/len3
                       !
                       len3 = dot_product(tmp_vec,tmp_vec2)
                       !
                       if (abs(len3) > 1.00D0) len3 = 1.0D0
                       !
                       ! The cosines needs to be 1, or very close to it 
                       !
                       if ((abs(len2)>0.99).and.(abs(len3)>0.99)) then
                         is_bond(i,2) = .true. 
                         r_FOD(i+fod_UP,:)=
     &(r_tmp(i+fod_UP,:)-r_atoms(ind1,:))/
     &(sqrt(sum((r_tmp(i+fod_UP,:)-r_atoms(ind1,:))**2)))*len1+
     &r_atoms(ind1,:)
                       end if
                     end if
                   end if
                 !
                 ! For any other bonding situation
                 !
                 else 
                   if ((abs(dist_tmp-dist1)<cut_bond).and.(j/=ind1)
     &.and.(is_bond(i,2).eqv..false.)) then           ! if distance to atom j is equal to dist1 : same distance to two atoms -> bond FOD
                     ind2 = j
                     counter = 0
                     !
                     ! Get center between atoms
                     !
                     bond_center(:) = 
     &(r_atoms(ind1,:)+r_atoms(ind2,:))/2.0D0
                     !
                     ! Get distance of initial FOD to this bond center
                     !
                     dist2 = sqrt(sum((r_tmp(i+fod_UP,:)-
     &bond_center(:))**2))
                     !
                     ! Think about a better way to do this !
                     ! Use the center between the atoms. Analyze
                     ! distance of the FODs to that center. They should
                     ! be very similar!
                     !
                     do k = 1, fod_DN
                       !
                       ! Check distance to FOD center
                       !
                       dist_tmp=sqrt(sum((r_tmp(k+fod_UP,:)-
     &bond_center(:))**2))
                       if ((abs(dist_tmp-dist2) < 0.05D0)) then
                         !
                         ! Check distance to corresponding atom
                         !
                         dist_tmp=sqrt(sum((r_tmp(k+fod_UP,:)-
     &r_atoms(ind1,:))**2))
                         if ((abs(dist_tmp-dist1) < 0.05D0)) then
                           !
                           ! Get information about the center of all bond
                           ! FODs
                           !
                           fod_center(:)=fod_center(:)+r_tmp(k+fod_UP,:)
                           !
                           ! Increase the counter for the bond FODs
                           !
                           counter = counter + 1
                           !
                           is_bond(k,2) = .true. 
                         end if
                       end if
                     end do
                     ! 
                     ! Calculate the actual bond center
                     !
                     fod_center(:) = fod_center(:)/counter
                     !
                     ! get correct new internal coordinate
                     !
                     len1 = r_internal(i,1,2)
                     len2 = r_internal(i,2,2)
                     !
                     ! If counter = 1: Single bond. Use new internal
                     ! coordinate to wriute new FOD
                     !
                     if (counter == 1) then
                     !
                       r_FOD(i+fod_UP,:) = 
     &(r_tmp(i+fod_UP,:)-r_atoms(ind1,:))/
     &(sqrt(sum((r_tmp(i+fod_UP,:)-r_atoms(ind1,:))**2)))*len1 +
     &r_atoms(ind1,:)

                     !!!!!!!!!
                     ! For other bond situations
                     !!!!!!
                     else
                       ! 
                       ! Now, get the new FOD positions
                       !
                       do k = 1, fod_DN
                         !
                         ! Check distance to FOD center
                         !
                         dist_tmp=sqrt(sum((r_tmp(k+fod_UP,:)-
     &bond_center(:))**2))
                         if ((abs(dist_tmp-dist2) < 0.05D0)) then
                           !
                           ! Check distance to corresponding atom
                           !
                           dist_tmp=sqrt(sum((r_tmp(k+fod_UP,:)-
     &r_atoms(ind1,:))**2))
                           if ((abs(dist_tmp-dist1) < 0.05D0)) then
                             !
                             ! For now: Use ratio between new internal
                             ! coordinate and old coordinate .. Do some more
                             !
                             ! Define ratio
                             !
                             !
                             ! For homonuclear bonding situations!
                             !
                             if (species(ind1)==species(ind2)) then
                               r_FOD(k+fod_UP,:) = 
     &(r_tmp(k+fod_UP,:)-fod_center(:))/
     &(sqrt(sum((r_tmp(k+fod_UP,:)-fod_center(:))**2)))*len1 +
     &fod_center(:)

                             !
                             ! For heteronuclear bonding situations
                             ! Need new_FOD-bond_center vector.
                             ! Rotate this vector into the other FODs.
                             ! Rotate !around! the bond axis
                             !
                             else
                               !
                               ! get new coordinate
                               ! 
                               r_FOD(k+fod_UP,:) = 
     &(r_tmp(k+fod_UP,:)-fod_center(:))/
     &(sqrt(sum((r_tmp(k+fod_UP,:)-fod_center(:))**2)))*len1 !!+ 
!!     &fod_center(:)


                               r_FOD(k+fod_UP,:) = r_FOD(k+fod_UP,:) +
     &(fod_center(:)-r_atoms(ind1,:))/
     &(sqrt(sum((fod_center(:)-r_atoms(ind1,:))**2)))*len2 +
     &r_atoms(ind1,:)

                             end if
                           end if
                         end if
                       end do
                     end if
                   end if
                 end if
               end do
             end if
           end do





           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Now, do all remaining FODs
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
           n_motifs(:) = 0
           do i = 1, fod_DN
             dist1 = 1000.0                                            ! initial value for the distance. Going to find minimum for all atoms
             do j = 1, n_atoms                                          ! Loop over atoms
               dist_tmp = sqrt(sum((r_tmp(i+fod_UP,:)-r_atoms(j,:))**2))
               if ((dist_tmp < dist1) .and. (dist_tmp < cutoff)) then  ! if distance to atom b is smaller than dist1 AND is is smaller than a cutoff radius
                 dist1 = dist_tmp                                      ! the cutoff radius ensures that we are only taking FODs for one atom 
                 ind1 = j
               end if
             end do
             !
             ! If it is not a bond FOD
             !
             if (is_bond(i,2).eqv..false.) then
               !
               ! Only if the fullForce contraint is suppossed to be used: Continue here
               !
               if ((constraint >= 2)) then
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                 ! Use center of lone FODs as bond_center -> keep lone
                 ! FODs on planes if they start on planes/lines (3 or 2
                 ! FODs). Tetrahedra will be conserved like this as well
                 !
                 fod_center = (/ 0.0D0, 0.0D0, 0.0D0 /)
                 !
                 ! If the structural motif has already been found (distance similar to the one found before within 0.033 bohr) -> don't do everything again
                 !
                 exist = .false.
                 do l = 1, 100
                   if (abs(dist1 - r_int1(ind1,l,2)) < 0.033) then    ! any distance FOR THIS ATOM. SPIN UP
                     exist = .true.
                     exit
                   end if
                 end do
                 !
                 ! If not -> new structural element
                 !
                 if ((.not.exist).and.
     &(r_internal(i,1,2) .ne. 1000.0D0)) then
                   counter = 0       ! count FODs in this structural element -> define 'bond_center'
                   do k = 1, fod_DN
                     if (is_bond(k,2).eqv..false.) then              ! do not do bond FODs
                       dist2 = sqrt(sum((r_tmp(k+fod_UP,:)-
     &r_atoms(ind1,:))**2))! check distance to the atom in question
                       if (abs(dist1 - dist2) < 0.033) then            ! If the new distance is the same as the old one               
                         counter = counter + 1
                         fod_center(:)=fod_center(:)+r_tmp(k+fod_UP,:)! bond center of the original positions!
                       end if
                     end if
                   end do
                   fod_center(:) = fod_center(:)/counter
                 end if
               end if

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Resymmetrize lone/core FODs
               !!!!
               ! Only if the fullForce contraint is suppossed to be used: Continue here
               !
               if ((constraint >= 2)) then
                 !
                 ! If the structural motif has already been found (distance similar to the one found before within 0.033 bohr) -> don't do everything again
                 !
                 exist = .false.
                 do l = 1, 100
                   if (abs(dist1 - r_int1(ind1,l,2)) < 0.033) then    ! any distance FOR THIS ATOM. SPIN UP
                     exist = .true.
                     exit
                   end if
                 end do
                 !
                 ! If not -> new structural element
                 !
                 if ((.not.exist).and.
     &(r_internal(i,1,2) .ne. 1000.0D0)) then
                   n_motifs(ind1) = n_motifs(ind1) + 1
                   r_int1(ind1,n_motifs(ind1),2) = dist1    ! Store the current identifier for the structural motif.
                   !
                   ! Use internal FOD position for all other FODs
                   ! NEED: rotation matrix between old internal r
                   ! (r_org_UP/r_org_DN) and new internal r 
                   !       -> use this properly to rotate new r_i towards new r_j
                   !
                   ! constrain symmetry to structural motifs
                   !

                   !
                   ! Idea: use center of the FODs to re-symmetrize, such
                   ! like for bond FODs? might give a better picture
                   !
                   !
                   ! Use this internal coordinate
                   !
                   len1 = r_internal(i,1,2)
                   len2 = r_internal(i,2,2)
                   ! 
                   ! if counter == 1 -> just use internal
                   ! cooridnate. Rescale old position to new size
                   !
                   if (counter==1) then
                     r_FOD(i+fod_UP,:) = 
     &(r_tmp(i+fod_UP,:)-r_atoms(ind1,:))/
     &(sqrt(sum((r_tmp(i+fod_UP,:)-r_atoms(ind1,:))**2)))*len1 +
     &r_atoms(ind1,:)
                     
                   !
                   ! If NOT a tetrahedron
                   !
                   else if (counter.ne.4) then
                     !
                     ! get new coords
                     !
                     do k = 1, fod_DN
                       if (is_bond(k,2).eqv..false.) then              ! do not do bond FODs
                         dist2 = 
     &sqrt(sum((r_tmp(k+fod_UP,:)-r_atoms(ind1,:))**2))! check distance to the atom in question
                         if (abs(dist1 - dist2) < 0.033) then          ! If the new distance is the same as the old one
                           r_FOD(k+fod_UP,:) = 
     &(r_tmp(k+fod_UP,:)-fod_center(:))/
     &(sqrt(sum((r_tmp(k+fod_UP,:)-fod_center(:))**2)))*len1 !! + 
!!     &fod_center(:)
                           !
                           ! If more than one atom
                           !
                           if (n_atoms > 1) then
                             r_FOD(k+fod_UP,:) = r_FOD(k+fod_UP,:) +
     &(fod_center(:)-r_atoms(ind1,:))/
     &(sqrt(sum((fod_center(:)-r_atoms(ind1,:))**2)))*len2 +
     &r_atoms(ind1,:)
                           end if
                         end if
                       end if
                     end do
      


                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ! For four points (tetrahedron)
                   ! Do breathing
                   !
                   else
                     len1 = r_internal(i,1,2)        ! New internal FOD position of FOD i
                     !
                     ! Apply new distance to all FODs
                     !
                     do k = 1, fod_DN                    ! Find the other internal coordinate
                       if (is_bond(k,2).eqv..false.) then
!!!     &.and.(k.ne.i))then                                              ! if this FOD has not been assigned yet (no bond)
                         dist2=
     &sqrt(sum((r_tmp(k+fod_UP,:)-r_atoms(ind1,:))**2))                 ! check distance to the atom in question
                         if (abs(dist1 - dist2) < 0.033) then   ! If the new distance is the same as the old one 
                           r_FOD(k+fod_UP,:) =
     &(r_tmp(k+fod_UP,:)-r_atoms(ind1,:))/dist2*len1 + r_atoms(ind1,:)
                         end if
                       end if
                     end do

                   end if
                 end if   ! end check whether structural element was already found
               end if     ! check constraint
             end if          ! check whether it is a bond or not
           end do         ! end 1. loop
         end if           ! constraint == 0 or not


         !
         ! Some post-porcessing for some constraint
         !
         !!!!! UP CHANNEL
         ! If only fix1s -> get coordinates (and forces) right here
         !
         if (constraint == 1) then
           do k = 1, fod_UP
             !
             ! No 1s !
             !
             if ((f_FOD(k,1).ne.0.0D0)
     &.and.(f_FOD(k,2).ne.0.0D0).and.(f_FOD(k,3).ne.0.0D0)) then
               r_FOD(k,1:3) = r_internal(k,1:3,1)
!!             else
!!               r_FOD(k,1:3) = r_tmp(k,1:3)
             end if
           end do
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! If freeFOD constraint -> restore the actual forces for the
         ! free FODs
         !!!!!!!!!!!!!!!1
         if (constraint == 4) then
           do k = 1, fod_UP
             !
             ! Do not take the 1s into the optimization
             !
             if ((freeFOD(k).eqv..true.)) then
               r_FOD(k,1:3) = r_internal(k,1:3,1)
             end if
           end do
         end if
         
         ! !!!! DN CNHANNEL
         ! If only fix1s -> get coordinates (and forces) right here
         !
         if (constraint == 1) then
           do k = 1, fod_DN
             !
             ! No 1s !
             !
             if ((f_FOD(k+fod_UP,1).ne.0.0D0)
     &.and.(f_FOD(k+fod_UP,2).ne.0.0D0)
     &.and.(f_FOD(k+fod_UP,3).ne.0.0D0)) then
               r_FOD(k+fod_UP,1:3) = r_internal(k,1:3,2)
!!!             else
!!!               r_FOD(k+fod_UP,1:3) = r_tmp(k+fod_UP,1:3)
             end if
           end do
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! If freeFOD constraint -> restore the actual forces for the
         ! free FODs
         !!!!!!!!!!!!!!!1
         if (constraint == 4) then
           do k = 1, fod_DN
             !
             ! Do not take the 1s into the optimization
             !
             if ((freeFOD(k+fod_UP).eqv..true.)) then
               r_FOD(k+fod_UP,1:3) = r_internal(k,1:3,2)
             end if
           end do
         end if




         !!!!!!!!!!!!!!!!!!!!!!!!!
         ! Final symmetrizations !
         !!!!!!!!!!!!!!!!!!!!!!!!!
         if ((constraint >= 2).and.(n_atoms > 1)) then
           !
           ! Symmetrize symmetry-equivalent FODs at
           ! different atoms/bonds. Same spin
           ! Use center of the molecule -> see whether positions are
           ! 'equivalent' with respect to that (i.e. same distance, same
           ! components...)
           !
           bond_center = (/ 0.0D0, 0.0D0, 0.0D0 /)    ! center of all atoms
           do i = 1, n_atoms
             bond_center(:) = bond_center(:) + r_atoms(i,:)
           end do
           bond_center(:) = bond_center(:)/n_atoms
 
!!!           !
!!!           ! GET ALL FODs with similar components -> average
!!!           !
!!!           !
!!!           ! UP Channel
!!!           !
!!!           do i = 1, fod_UP
!!!             vec1(:) = r_FOD(i,:) - bond_center(:)
!!!             !
!!!             ! See how many similar FODs we have
!!!             ! Store all values in tmp_vec
!!!             !
!!!             counter = 1
!!!             tmp_vec(:) = abs(vec1(:))
!!!             do j = 1, fod_UP
!!!               !
!!!               ! No need to evaluate an FOD with itself
!!!               !
!!!               if (i.ne.j) then
!!!                 vec2(:) = r_FOD(j,:) - bond_center(:)
!!!                 !
!!!                 ! If the absolute components are very similar
!!!                 !
!!!                 if (abs(sum(abs(vec1(:))-abs(vec2(:)))) < 1.0D-2) then
!!!                   counter = counter + 1
!!!                   tmp_vec(:) = tmp_vec(:) + abs(vec2(:))
!!!                 end if
!!!               end if
!!!             end do
!!!             !
!!!             ! If  counter > 1 -> do something
!!!             !
!!!             if (counter > 1) then
!!!               !
!!!               ! After all is done -> average coordinates and apply them
!!!               !
!!!               tmp_vec(:) = tmp_vec(:)/counter
!!!               !
!!!               do j = 1, fod_UP
!!!                 !
!!!                 ! Include first vector as well!
!!!                 !
!!!                 vec2(:) = r_FOD(j,:) - bond_center(:)
!!!                 !
!!!                 ! If the absolute components are very similar
!!!                 !
!!!                 if (abs(sum(abs(vec1(:))-abs(vec2(:)))) < 1.0D-2) then
!!!                   !
!!!                   ! Use the average positions. Conserve the sign !!
!!!                   ! 
!!!                   r_FOD(j,:) = tmp_vec(:)*sign(1.0D0,r_FOD(j,:))
!!!                   r_FOD(j,:) = r_FOD(j,:) + bond_center(:)
!!!                 end if
!!!               end do
!!!             end if
!!!           end do
!!!             
!!!
!!!           !
!!!           ! DN Channel
!!!           !
!!!           do i = 1, fod_DN
!!!             vec1(:) = r_FOD(fod_UP+i,:) - bond_center(:)
!!!             !
!!!             ! See how many similar FODs we have
!!!             ! Store all values in tmp_vec
!!!             !
!!!             counter = 1
!!!             tmp_vec(:) = abs(vec1(:))
!!!             do j = 1, fod_DN
!!!               !
!!!               ! No need to evaluate an FOD with itself
!!!               !
!!!               if (i.ne.j) then
!!!                 vec2(:) = r_FOD(fod_UP+j,:) - bond_center(:)
!!!                 !
!!!                 ! If the absolute components are very similar
!!!                 !
!!!                 if (abs(sum(abs(vec1(:))-abs(vec2(:)))) < 1.0D-2) then
!!!                   counter = counter + 1
!!!                   tmp_vec(:) = tmp_vec(:) + abs(vec2(:))
!!!                 end if
!!!               end if
!!!             end do
!!!             !
!!!             ! If counter > 1
!!!             !
!!!             if (counter > 1) then
!!!               !
!!!               ! After all is done -> average coordinates and apply them
!!!               !
!!!               tmp_vec(:) = tmp_vec(:)/counter
!!!               !
!!!               do j = 1, fod_DN
!!!                 !
!!!                 ! Include first vector as well!
!!!                 !
!!!                 vec2(:) = r_FOD(fod_UP+j,:) - bond_center(:)
!!!                 !
!!!                 ! If the absolute components are very similar
!!!                 !
!!!                 if (abs(sum(abs(vec1(:))-abs(vec2(:)))) < 1.0D-2) then
!!!                   !
!!!                   ! Use the average positions. Conserve the sign !!
!!!                   ! 
!!!                   r_FOD(fod_UP+j,:) = 
!!!     &tmp_vec(:)*sign(1.0D0,r_FOD(fod_UP+j,:))
!!!                   r_FOD(fod_UP+j,:)=r_FOD(fod_UP+j,:) + bond_center(:)
!!!                 end if
!!!               end do
!!!             end if
!!!           end do

           !!!!!!!!!!!!!!!!!
           ! In case of spin-paired FODs -> put UP and DN at exactly the same
           ! position after optimization
           !!!!!!!!!!!!!!!!!
           do i = 1, fod_UP
             do j = 1, fod_DN
               !
               ! Do not modify free FODs that were being optimized
               !
               if ((constraint==4).and.
     &((freeFOD(i).eqv..true.).or.(freeFOD(j+fod_UP).eqv..true.))) then
               else
                 len1 = sqrt(sum((r_FOD(i,:)-r_FOD(fod_UP+j,:))**2))
                 if (len1 < 2.0D-1) then
                   !
                   ! Average positions
                   !
                   r_FOD(i,:) = (r_FOD(i,:)+r_FOD(fod_UP+j,:))/2.0
                   r_FOD(fod_UP+j,:) = r_FOD(i,:)
                 end if
               end if
             end do
           end do

         end if



         !
         ! Write results. Again, assume that only FRMORB exists
         ! These are already re-scaled if scaled-optimization is used
         !
         open(90,file=FODFILESTR,status='old',action='write')
         write(90,*) fod_UP, fod_DN
         do i = 1, fod_UP+fod_DN
           ! 
           ! Write new FOD positions to FRMORB
           !
           if (constraint >= 2) then
             write(90,99) r_FOD(i,1:3)
           else
             write(90,98) r_FOD(i,1:3)
           end if
         end do
         close(90)

         deallocate(n_motifs)
         deallocate(r_int1)
         deallocate(r_bond)
         deallocate(r_FOD)
         deallocate(f_FOD)
         deallocate(r_all)
         deallocate(f_all)
         deallocate(r_tmp)
         deallocate(r_org_int_UP)
         deallocate(r_org_int_DN)
         deallocate(r_internal)
         deallocate(f_internal)
         deallocate(r_atoms)
         deallocate(species)
         if (scaled) deallocate(hess_D)
         if (constraint==4) deallocate(freeFOD)

         !call stopit
 777     format(A,A,A,A)
 778     format(A,A,A)
 98      format(3F20.12)
 99      format(3F12.7)
         !
         ! Last call for timing
         !
         call cpu_time(finish)
         !
         ! Write timing
         !
         write(6,*) 'Total CPU time for optimizer ',finish-start,' s'
         call flush(6)
         end






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to generate rotation matrix !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine get_rotMat(rotMat,alpha,axis)
          real(8), intent(inout) :: rotMat(3,3)
          real(8), intent(in)    :: alpha
          real(8), intent(in)    :: axis(3)
          
          rotMat(1,1)=cos(alpha) + axis(1)**2*(1 - cos(alpha))
          rotMat(1,2)=axis(1)*axis(2)*(1-cos(alpha))-axis(3)*sin(alpha)
          rotMat(1,3)=axis(1)*axis(3)*(1-cos(alpha))+axis(2)*sin(alpha)
          rotMat(2,1)=axis(2)*axis(1)*(1-cos(alpha))+axis(3)*sin(alpha)
          rotMat(2,2)=cos(alpha) + axis(2)**2*(1 - cos(alpha))
          rotMat(2,3)=axis(2)*axis(3)*(1-cos(alpha))-axis(1)*sin(alpha)
          rotMat(3,1)=axis(3)*axis(1)*(1-cos(alpha))-axis(2)*sin(alpha)
          rotMat(3,2)=axis(3)*axis(2)*(1-cos(alpha))+axis(1)*sin(alpha)
          rotMat(3,3)=cos(alpha) + axis(3)**2*(1 - cos(alpha))

          return
        end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to compare internal coordinates. !
! And average them if necessary               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine comp_int(r1,r2,f1,f2)
          real(8), intent(inout) :: r1, r2, f1,f2
          real(8)                :: lenn

          !
          ! If force component is set to 0 -> leave it
          !
          if ((f1 == 0.0D0).or.(f2 == 0.0D0)) then
            if ((abs(r1-r2) < 5.0D-4)) then
              !
              ! Don't take internal coordinates that are not assigned
              ! (= 1000.0D0)
              !
              if ((r1<10.0D0).and.(r2<10.0D0)) then
                lenn = (r1+r2)/2.0D0
                r1 = lenn
                r2 = lenn
                f1 = 0.0D0
                f2 = 0.0D0
              end if
            end if
          
          else
            if ((abs(r1-r2) < 5.0D-4).and.(abs(f1-f2) < 5.0D-4)) then
              ! 
              ! If coordinate and force are already
              ! known as internal coordinates 
              ! -> reset the current one.
              ! Average the forces and coordinates
              !
              lenn = (r1+r2)/2.0D0
              r1 = lenn
              r2 = lenn
              lenn = (f1+f2)/2.0D0
              f1 = lenn
              f2 = lenn
            end if
          end if

          return
        end subroutine

C
C Originally by Kushantha Withanage (CMU)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine will write the inverse of the estimated 
c diagonal elements of the Hessian
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine write_HD_constr()
        real*8  rf,ra,HD,
     &  r,D,tab,DD,HS
        dimension rf(3,1000),ra(3,1000)
        dimension HD(1000),D(1000)
        dimension tab(36,5)
        integer nup,ndn,na,nz(1000),id
        character(6) :: FODFILESTR
        logical :: exist

        inquire(file='FRMIDT',exist=exist)
        if(.not.exist)then
          write(FODFILESTR,'(A)')'FRMORB'
        else
          write(FODFILESTR,'(A)')'FRMIDT'
        end if



        open(4,file=FODFILESTR)
        read(4,*) nup,ndn
        do i=1,nup+ndn
         read(4,*)(rf(j,i),j=1,3)
        end do
        close(4)
C tab contains the distance of the shells of FODs
C tab  is currntly read in as an input file Table based on previous
C FLOSIC calculations.
        open(unit=3,file='Table')
        do i=1,36
         read(3,*)(tab(i,j),j=1,5)
        end do
        close(3)
c reading atomic data
        open(unit=2,file='XMOL.DAT')
        read(2,*) na
        read(2,*)
        do i=1,na
         read(2,*)nz(i),(ra(j,i),j=1,3)
        end do
        close(2)

c converting atomic coordinates to bohr

        do i=1,na
         do j=1,3
          ra(j,i)=ra(j,i)*1.889725989d0
         end do
        end do
c start of the process of finding right inversed Hessian element
c finding the closest atom to each FOD
        do i=1,nup+ndn
         DD=1000.D0
         id=1
         do j=1,na
          D(j)=dsqrt((ra(1,j)-rf(1,i))**2 + (ra(2,j)-rf(2,i))**2
     &  +(ra(3,j)-rf(3,i))**2)
          if (D(j).lt.DD) then
           id=j
           DD=D(j)
          end if
         end do
         print *,DD
         call chooseit_constr(i,nz(id),HS,DD,tab)
         HD(i)=HS

        end do
c writing the  inversed Hessian diagonals
        open(unit=2,file='HESS_DIAG')
        do i=1,nup+ndn
         write(2,*) dsqrt(1.D0/HD(i))
         write(2,*) dsqrt(1.D0/HD(i))
         write(2,*) dsqrt(1.D0/HD(i))
        end do
        close(2)
        end  subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine chooseit_constr(i,nz,HS,DD,tab)
        real*8 HS,tab,DD,HESS_constr
        dimension tab(36,5)
        integer  nz,is,i

c lets find the shell
        do is=1,5
         if(DD.LT.tab(nz,is)) then
          HS= HESS_constr(nz,is)
          go to  100
         end if
        end do
 100    end subroutine

        function HESS_constr(nz,is)
        real*8 HESS_constr,z
        integer nz,is
        z=dble(nz)
        if(is.eq.1)then
         if(nz.lt.10) then
          HESS_constr=1.D0*z/10.D0
         else
          HESS_constr= 0.0029d0*z**3-0.0526d0*z**2 + 0.4*z-0.8762D0
         end if
        end if
        if(is.eq.2.and.nz.lt.10) HESS_constr=0.021D0*z/10.D0
        if(is.eq.2.and.nz.ge.10) then
         HESS_constr=0.0012d0*z**3 -0.0382d0*z**2+0.4146D0*z-1.442d0
        end if
        if(is.eq.3.and.nz.ge.18) then
         HESS_constr=0.00009d0*z**3-0.004d0*z**2+0.0623D0*z - 0.3225d0
        end if
        if(is.eq.3.and.nz.lt.18) HESS_constr=0.0124D0*z/18.D0
        if(is.eq.4)HESS_constr=0.0022D0*z/36.0D0
        end function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
C Subroutine to write Table file required for scaled LBFGS. 
C copied from write_table.f90
C
        subroutine write_table_constr
        logical :: EXIST

        inquire(file='Table',EXIST=EXIST)
        if(EXIST) then
          return
        else
          open(86,file='Table')
          write(86,'(A)') "0.50 100.00 101.00 102.00 103.00"
          write(86,'(A)') "0.50 100.00 101.00 102.00 103.00"
          write(86,'(A)') "0.50 4.00   101.00 102.00 103.00"
          write(86,'(A)') "0.50 4.00   101.00 102.00 103.00"
          write(86,'(A)') "0.40 4.00   101.00 102.00 103.00"
          write(86,'(A)') "0.30 4.00   101.00 102.00 103.00"
          write(86,'(A)') "0.30 3.50   101.00 102.00 103.00"
          write(86,'(A)') "0.30 2.50   101.00 102.00 103.00"
          write(86,'(A)') "0.30 2.00   101.00 102.00 103.00"
          write(86,'(A)') "0.30 2.50   101.00 102.00 103.00"
          write(86,'(A)') "0.20 2.09   8.00   102.00 103.00"
          write(86,'(A)') "0.20 1.83   8.00   102.00 103.00"
          write(86,'(A)') "0.10 1.33   5.00   102.00 103.00"
          write(86,'(A)') "0.10 1.19   5.00   102.00 103.00"
          write(86,'(A)') "0.05 1.05   4.00   102.00 103.00"
          write(86,'(A)') "0.05 0.88   4.00   102.00 103.00"
          write(86,'(A)') "0.05 0.81   3.00   102.00 103.00"
          write(86,'(A)') "0.04 0.72   3.00   102.00 103.00"
          write(86,'(A)') "0.03 0.66   2.43   7.00   103.00"
          write(86,'(A)') "0.03 0.63   2.23   7.00   103.00"
          write(86,'(A)') "0.03 0.60   2.28   7.00   103.00"
          write(86,'(A)') "0.02 0.54   2.24   7.00   103.00"
          write(86,'(A)') "0.02 0.53   2.00   7.00   103.00"
          write(86,'(A)') "0.02 0.48   2.21   7.00   103.00"
          write(86,'(A)') "0.02 0.46   2.46   7.00   103.00"
          write(86,'(A)') "0.02 0.78   2.14   7.00   103.00"
          write(86,'(A)') "0.02 0.61   2.11   7.00   103.00"
          write(86,'(A)') "0.02 0.39   1.86   7.00   103.00"
          write(86,'(A)') "0.02 0.55   2.34   8.00   103.00"
          write(86,'(A)') "0.01 0.37   1.94   7.00   103.00"
          write(86,'(A)') "0.01 0.36   1.75   7.00   103.00"
          write(86,'(A)') "0.01 0.34   1.49   7.00   103.00"
          write(86,'(A)') "0.01 0.33   1.41   7.00   103.00"
          write(86,'(A)') "0.01 0.31   1.38   6.00   103.00"
          write(86,'(A)') "0.01 0.30   1.22   6.00   103.00"
          write(86,'(A)') "0.01 0.28   1.14   6.00   103.00"
          close(86)
          return
        end if
        end subroutine write_table_constr


C##############################################################
C 9/15/2017
C  Slightly modified lbfgs routine for FOD position optimization.
C   -YY
C     ----------------------------------------------------------------------
C     This file contains the LBFGS algorithm and supporting routines
C
C     ****************
C     FOLBFGS SUBROUTINE
C     ****************
C
      SUBROUTINE FOLBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPSL,XTOL,W,
     & IFLAG,DGUESS,OBJECT)
C
      INTEGER N,M,IPRINT(2),IFLAG
      DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
      DOUBLE PRECISION F,EPSL,XTOL,DGUESS
      LOGICAL DIAGCO
      ! KT
      INTEGER OBJECT  ! determine which internal coordinate to optimize
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C 
C     This subroutine solves the unconstrained minimization problem
C 
C                      min F(x),    x= (x1,x2,...,xN),
C
C      using the limited memory BFGS method. The routine is especially
C      effective on problems involving a large number of variables. In
C      a typical iteration of this method an approximation Hk to the
C      inverse of the Hessian is obtained by applying M BFGS updates to
C      a diagonal matrix Hk0, using information from the previous M steps.
C      The user specifies the number M, which determines the amount of
C      storage required by the routine. The user may also provide the
C      diagonal matrices Hk0 if not satisfied with the default choice.
C      The algorithm is described in "On the limited memory BFGS method
C      for large scale optimization", by D. Liu and J. Nocedal,
C      Mathematical Programming B 45 (1989) 503-528.
C 
C      The user is required to calculate the function value F and its
C      gradient G. In order to allow the user complete control over
C      these computations, reverse  communication is used. The routine
C      must be called repeatedly under the control of the parameter
C      IFLAG. 
C
C      The steplength is determined at each iteration by means of the
C      line search routine MCVSRCH, which is a slight modification of
C      the routine CSRCH written by More' and Thuente.
C 
C      The calling statement is 
C 
C          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPSL,XTOL,W,IFLAG)
C 
C      where
C 
C     N       is an INTEGER variable that must be set by the user to the
C             number of variables. It is not altered by the routine.
C             Restriction: N>0.
C 
C     M       is an INTEGER variable that must be set by the user to
C             the number of corrections used in the BFGS update. It
C             is not altered by the routine. Values of M less than 3 are
C             not recommended; large values of M will result in excessive
C             computing time. 3<= M <=7 is recommended. Restriction: M>0.
C 
C     X       is a DOUBLE PRECISION array of length N. On initial entry
C             it must be set by the user to the values of the initial
C             estimate of the solution vector. On exit with IFLAG=0, it
C             contains the values of the variables at the best point
C             found (usually a solution).
C 
C     F       is a DOUBLE PRECISION variable. Before initial entry and on
C             a re-entry with IFLAG=1, it must be set by the user to
C             contain the value of the function F at the point X.
C 
C     G       is a DOUBLE PRECISION array of length N. Before initial
C             entry and on a re-entry with IFLAG=1, it must be set by
C             the user to contain the components of the gradient G at
C             the point X.
C 
C     DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
C             user  wishes to provide the diagonal matrix Hk0 at each
C             iteration. Otherwise it should be set to .FALSE., in which
C             case  LBFGS will use a default value described below. If
C             DIAGCO is set to .TRUE. the routine will return at each
C             iteration of the algorithm with IFLAG=2, and the diagonal
C              matrix Hk0  must be provided in the array DIAG.
C 
C 
C     DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE.,
C             then on initial entry or on re-entry with IFLAG=2, DIAG
C             it must be set by the user to contain the values of the 
C             diagonal matrix Hk0.  Restriction: all elements of DIAG
C             must be positive.
C 
C     IPRINT  is an INTEGER array of length two which must be set by the
C             user.
C 
C             IPRINT(1) specifies the frequency of the output:
C                IPRINT(1) < 0 : no output is generated,
C                IPRINT(1) = 0 : output only at first and last iteration,
C                IPRINT(1) > 0 : output every IPRINT(1) iterations.
C 
C             IPRINT(2) specifies the type of output generated:
C                IPRINT(2) = 0 : iteration count, number of function 
C                                evaluations, function value, norm of the
C                                gradient, and steplength,
C                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
C                                variables and  gradient vector at the
C                                initial point,
C                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
C                                variables,
C                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
C 
C 
C     EPSL     is a positive DOUBLE PRECISION variable that must be set by
C             the user, and determines the accuracy with which the solution
C             is to be found. The subroutine terminates when
C
C                         ||G|| < EPSL max(1,||X||),
C             Changed to RMS < EPSL by DJW.
C
C             where ||.|| denotes the Euclidean norm.
C 
C     XTOL    is a  positive DOUBLE PRECISION variable that must be set by
C             the user to an estimate of the machine precision (e.g.
C             10**(-16) on a SUN station 3/60). The line search routine will
C             terminate if the relative width of the interval of uncertainty
C             is less than XTOL.
C 
C     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
C             workspace for LBFGS. This array must not be altered by the
C             user.
C 
C     IFLAG   is an INTEGER variable that must be set to 0 on initial entry
C             to the subroutine. A return with IFLAG<0 indicates an error,
C             and IFLAG=0 indicates that the routine has terminated without
C             detecting errors. On a return with IFLAG=1, the user must
C             evaluate the function F and gradient G. On a return with
C             IFLAG=2, the user must provide the diagonal matrix Hk0.
C 
C             The following negative values of IFLAG, detecting an error,
C             are possible:
C 
C              IFLAG=-1  The line search routine MCSRCH failed. The
C                        parameter INFO provides more detailed information
C                        (see also the documentation of MCSRCH):
C
C                       INFO = 0  IMPROPER INPUT PARAMETERS.
C
C                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
C                                 UNCERTAINTY IS AT MOST XTOL.
C
C                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
C                                 REQUIRED AT THE PRESENT ITERATION.
C
C                       INFO = 4  THE STEP IS TOO SMALL.
C
C                       INFO = 5  THE STEP IS TOO LARGE.
C
C                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
C                                 THERE MAY NOT BE A STEP WHICH SATISFIES
C                                 THE SUFFICIENT DECREASE AND CURVATURE
C                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
C
C 
C              IFLAG=-2  The i-th diagonal element of the diagonal inverse
C                        Hessian approximation, given in DIAG, is not
C                        positive.
C           
C              IFLAG=-3  Improper input parameters for LBFGS (N or M are
C                        not positive).
C 
C
C
C    ON THE DRIVER:
C
C    The program that calls LBFGS must contain the declaration:
C
C                       EXTERNAL LB2
C
C    LB2 is a BLOCK DATA that defines the default values of several
C    parameters described in the COMMON section. 
C
C 
C 
C    COMMON:
C 
C     The subroutine contains one common area, which the user may wish to
C    reference:
C 
         COMMON /LB3/MP,LP,GLTOL,STPMIN,STPMAX
C 
C    MP  is an INTEGER variable with default value 6. It is used as the
C        unit number for the printing of the monitoring information
C        controlled by IPRINT.
C 
C    LP  is an INTEGER variable with default value 6. It is used as the
C        unit number for the printing of error messages. This printing
C        may be suppressed by setting LP to a non-positive value.
C 
C    GLTOL is a DOUBLE PRECISION variable with default value 0.9, which
C        controls the accuracy of the line search routine MCSRCH. If the
C        function and gradient evaluations are inexpensive with respect
C        to the cost of the iteration (which is sometimes the case when
C        solving very large problems) it may be advantageous to set GLTOL
C        to a small value. A typical small value is 0.1.  Restriction:
C        GLTOL should be greater than 1.D-04.
C 
C    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which
C        specify lower and uper bounds for the step in the line search.
C        Their default values are 1.D-20 and 1.D+20, respectively. These
C        values need not be modified unless the exponents are too large
C        for the machine being used, or unless the problem is extremely
C        badly scaled (in which case the exponents should be increased).
C 
C
C  MACHINE DEPENDENCIES
C
C        The only variables that are machine-dependent are XTOL,
C        STPMIN and STPMAX.
C 
C
C  GENERAL INFORMATION
C 
C    Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH
C 
C    Input/Output  :  No input; diagnostic messages on unit MP and
C                     error messages on unit LP.
C 
C 
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      DOUBLE PRECISION GLTOL,ONE,ZERO,GNORM,DDOT,STP1,FLTOL,STPMIN,
     &                 STPMAX,STP,YS,YY,SQ,YR,BETA,XNORM
      INTEGER MP,LP,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,INFO,
     &        BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN,L
      LOGICAL FINISH, EXIST
C
      SAVE
      DATA ONE,ZERO/1.0D+0,0.0D+0/
C
C     INITIALIZE
C     ----------
C
! KT

      character(len=15) :: fdiag, fsearch, fstep ! names of files to distinguish between objects
      !
      ! For shellwise OPT -> define file names
      !
      if (OBJECT < 10) then
        write(fdiag,'(A5,I1,A4)')   'FDIAG',OBJECT,'.LBF'
        write(fsearch,'(A7,I1,A4)') 'FSEARCH',OBJECT,'.LBF'
        write(fstep,'(A5,I1,A4)')   'FSTEP',OBJECT,'.LBF'
      else if (object > 10) then
        write(fdiag,'(A5,I2,A4)')   'FDIAG',OBJECT,'.LBF'
        write(fsearch,'(A7,I2,A4)') 'FSEARCH',OBJECT,'.LBF'
        write(fstep,'(A5,I2,A4)')   'FSTEP',OBJECT,'.LBF'
      end if

**********************CHANGES MADE HERE*********************
C IF GRADIENTS HAVE CONVERGED THEN EXIT
        IF (N .LT. 1) THEN
         GOTO 300
        ENDIF

        INQUIRE(FILE=fsearch,EXIST=EXIST)
        IF (EXIST) THEN
         OPEN(UNIT=7,FILE=fsearch,STATUS='OLD')
         READ(7,*)ITER
         READ(7,*)INFO
         READ(7,*)NFEV
         READ(7,*)ISPT
         READ(7,*)POINT
         READ(7,*)STP
         READ(7,*)IYPT
         READ(7,*)NPT
         READ(7,*)CP
         READ(7,*)NFUN
         FLTOL= 1.0D-4
         MAXFEV= 20
         CLOSE(UNIT=7)
*	 write(6,*)"FINISHED READING SRCH.DAT..."
        ENDIF
***************************TILL HERE**************************************** 

      IF(IFLAG.EQ.0) GO TO 10
      GO TO (172,100) IFLAG
  10  ITER= 0
      IF(N.LE.0.OR.M.LE.0) GO TO 196
      IF(GLTOL.LE.1.D-04) THEN
        IF(LP.GT.0) WRITE(LP,245)
        GLTOL=9.D-01
      ENDIF
      NFUN= 1
      POINT= 0
      FINISH= .FALSE.
      IF(DIAGCO) THEN
         DO 30 I=1,N
 30      IF (DIAG(I).LE.ZERO) GO TO 195
      ELSE
         DO 40 I=1,N
 40      DIAG(I)=DGUESS
      ENDIF
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
      ISPT= N+2*M
      IYPT= ISPT+N*M     
      DO 50 I=1,N
 50   W(ISPT+I)= -G(I)*DIAG(I)
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      STP1= ONE/GNORM
C
C     PARAMETERS FOR LINE SEARCH ROUTINE
C     
      FLTOL= 1.0D-4
      MAXFEV= 20
C
      IF(IPRINT(1).GE.0) CALL LB1(IPRINT,ITER,NFUN,
     &                     GNORM,N,M,X,F,G,STP,FINISH)
C
C    --------------------
C     MAIN ITERATION LOOP
C    --------------------
C
 80   ITER= ITER+1
      INFO=0
      BOUND=ITER-1
      IF(ITER.EQ.1) GO TO 165
      IF (ITER .GT. M)BOUND=M
C
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
      IF(.NOT.DIAGCO) THEN
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         DO 90 I=1,N
   90    DIAG(I)= YS/YY
      ELSE
         IFLAG=2
         RETURN
      ENDIF
 100  CONTINUE
      IF(DIAGCO) THEN
        DO 110 I=1,N
 110    IF (DIAG(I).LE.ZERO) GO TO 195
      ENDIF
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
      CP= POINT
      IF (POINT.EQ.0) CP=M
      W(N+CP)= ONE/YS
      DO 112 I=1,N
 112  W(I)= -G(I)
      CP= POINT
      DO 125 I= 1,BOUND
         CP=CP-1
         IF (CP.EQ. -1)CP=M-1
         SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
         INMC=N+M+CP+1
         IYCN=IYPT+CP*N
         W(INMC)= W(N+CP+1)*SQ
         CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
 125  CONTINUE
C
      DO 130 I=1,N
 130  W(I)=DIAG(I)*W(I)
C
      DO 145 I=1,BOUND
         YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
         BETA= W(N+CP+1)*YR
         INMC=N+M+CP+1
         BETA= W(INMC)-BETA
         ISCN=ISPT+CP*N
         CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
         CP=CP+1
         IF (CP.EQ.M)CP=0
 145  CONTINUE
C
C     STORE THE NEW SEARCH DIRECTION
C     ------------------------------
C
       DO 160 I=1,N
 160   W(ISPT+POINT*N+I)= W(I)
C
C     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION 
C     BY USING THE LINE SEARCH ROUTINE MCSRCH
C     ----------------------------------------------------
 165  NFEV=0
      STP=ONE
      IF (ITER.EQ.1) STP=STP1
      DO 170 I=1,N
 170  W(I)=G(I)
 172  CONTINUE
      CALL FOMCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FLTOL,
     &            XTOL,MAXFEV,INFO,NFEV,DIAG,OBJECT)
*****************CHANGES MADE HERE*****************
        OPEN(UNIT=7,FILE=fsearch,STATUS='UNKNOWN')
        WRITE(7,*)ITER, " ITER"
        WRITE(7,*)INFO, " INFO"
        WRITE(7,*)NFEV, " NFEV"  
        WRITE(7,*)ISPT, " ISPT"
        WRITE(7,*)POINT," POINT"
        WRITE(7,*)STP, " STP"
        WRITE(7,*)IYPT, " IYPT"
        WRITE(7,*)NPT, " NPT"
        WRITE(7,*)CP, " CP"
        WRITE(7,*)NFUN, " NFUN"
        CLOSE(UNIT=7)
*****************TILL HERE****************************
      IF (INFO .EQ. -1) THEN
        IFLAG=1
        RETURN
      ENDIF
      IF (INFO .NE. 1) GO TO 190
      NFUN= NFUN + NFEV
C
C     COMPUTE THE NEW STEP AND GRADIENT CHANGE 
C     -----------------------------------------
C
      NPT=POINT*N
      DO 175 I=1,N
      W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
 175  W(IYPT+NPT+I)= G(I)-W(I)
      POINT=POINT+1
      IF (POINT.EQ.M)POINT=0
C
C     TERMINATION TEST
C     ----------------
C
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      XNORM= DSQRT(DDOT(N,X,1,X,1))
      XNORM= DMAX1(1.0D0,XNORM)
C     IF (GNORM/XNORM .LE. EPSL) FINISH=.TRUE.
C  Changed the convergence criterion to RMS. DJW
      IF (GNORM/SQRT(1.0D0*N) .LE. EPSL) FINISH=.TRUE.
C
      IF(IPRINT(1).GE.0) CALL LB1(IPRINT,ITER,NFUN,
     &               GNORM,N,M,X,F,G,STP,FINISH)
      IF (FINISH) THEN
         IFLAG=0
***************CHANGE MADE HERE***********************
*         CLOSE(UNIT=7,STATUS='DELETE')
****************TILL HERE***************************
         RETURN
      ENDIF
      GO TO 80
C
C     ------------------------------------------------------------
C     END OF MAIN ITERATION LOOP. ERROR EXITS.
C     ------------------------------------------------------------
C
 190  IFLAG=-1
C     IF(LP.GT.0) WRITE(LP,200) INFO
      IF(LP.GT.0) WRITE(LP,200) INFO
      RETURN
 195  IFLAG=-2
      IF(LP.GT.0) WRITE(LP,235) I
      RETURN
 196  IFLAG= -3
      IF(LP.GT.0) WRITE(LP,240)
C
C     FORMATS
C     -------
C
 200  FORMAT(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE'
     @          ' DOCUMENTATION OF ROUTINE MCSRCH',/' ERROR RETURN'
     @          ' OF LINE SEARCH: INFO= ',I2,/
     @          ' POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT',/,
     @          ' OR INCORRECT TOLERANCES')
 235  FORMAT(/' IFLAG= -2',/' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     @       ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
 240  FORMAT(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N OR M',
     @       ' ARE NOT POSITIVE)')
 245  FORMAT(/'  GLTOL IS LESS THAN OR EQUAL TO 1.D-04',
     @       / ' IT HAS BEEN RESET TO 9.D-01')
      RETURN

 300   IF ((IFLAG .EQ. 0) .OR. (N .EQ. 0)) THEN
        OPEN(UNIT=4,FILE=fdiag,FORM='formatted',STATUS='unknown')
        CLOSE(4,STATUS='delete')
        OPEN(UNIT=7,FILE=fsearch,FORM='formatted',
     &                                 STATUS='unknown')
        CLOSE(7,STATUS='DELETE')
        OPEN(UNIT=8,FILE=fstep, FORM='formatted',STATUS='unknown')
        CLOSE(8,STATUS='DELETE')
       END IF
       RETURN
      END
C#################################################
C    ------------------------------------------------------------------
C
C     **************************
C     LINE SEARCH ROUTINE MCSRCH
C     **************************
C
      SUBROUTINE FOMCSRCH(N,X,F,G,S,STP,FLTOL,XTOL,MAXFEV,INFO,NFEV,WA,
     &OBJECT)
      INTEGER N,MAXFEV,INFO,NFEV
      DOUBLE PRECISION F,STP,FLTOL,GLTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
      COMMON /LB3/MP,LP,GLTOL,STPMIN,STPMAX
      ! KT
      INTEGER OBJECT
      SAVE
C
C                     SUBROUTINE MCSRCH
C                
C     A slight modification of the subroutine CSRCH of More' and Thuente.
C     The changes are to allow reverse communication, and do not affect
C     the performance of the routine. 
C
C     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
C     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
C
C     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
C     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
C     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
C     MINIMIZER OF THE MODIFIED FUNCTION
C
C          F(X+STP*S) - F(X) - FLTOL*STP*(GRADF(X)'S).
C
C     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
C     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
C     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
C     CONTAINS A MINIMIZER OF F(X+STP*S).
C
C     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
C     THE SUFFICIENT DECREASE CONDITION
C
C           F(X+STP*S) .LE. F(X) + FLTOL*STP*(GRADF(X)'S),
C
C     AND THE CURVATURE CONDITION
C
C           ABS(GRADF(X+STP*S)'S)) .LE. GLTOL*ABS(GRADF(X)'S).
C
C     IF FLTOL IS LESS THAN GLTOL AND IF, FOR EXAMPLE, THE FUNCTION
C     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
C     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
C     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
C     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
C     SATISFIES THE SUFFICIENT DECREASE CONDITION.
C
C     THE SUBROUTINE STATEMENT IS
C
C        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FLTOL,XTOL, MAXFEV,INFO,NFEV,WA)
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
C         X + STP*S.
C
C       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
C         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
C
C       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
C         OF F AT X + STP*S.
C
C       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
C         SEARCH DIRECTION.
C
C       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
C         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
C         STP CONTAINS THE FINAL ESTIMATE.
C
C       FLTOL AND GLTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
C         communication implementation GLTOL is defined in a COMMON
C         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
C         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
C         SATISFIED.
C
C       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
C         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C         IS AT MOST XTOL.
C
C       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
C         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
C         communication implementatin they are defined in a COMMON
C         statement).
C
C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
C         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
C         MAXFEV BY THE END OF AN ITERATION.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C
C         INFO = 0  IMPROPER INPUT PARAMETERS.
C
C         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
C
C         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
C                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
C
C         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C                   IS AT MOST XTOL.
C
C         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
C
C         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
C
C         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
C
C         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
C                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
C                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
C                   TOLERANCES MAY BE TOO SMALL.
C
C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
C         CALLS TO FCN.
C
C       WA IS A WORK ARRAY OF LENGTH N.
C
C     SUBPROGRAMS CALLED
C
C       MCSTEP
C
C       FORTRAN-SUPPLIED...ABS,MAX,MIN
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
C     JORGE J. MORE', DAVID J. THUENTE
C
C     **********
      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1,EXIST
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,
     &       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,
     &       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/

! KT
      character(len=15) :: fstep ! names of files to distinguish between objects
        !
        ! For shellwise OPT -> define file names
        !
        if (OBJECT < 10) then
          write(fstep,'(A5,I1,A4)')   'FSTEP',OBJECT,'.LBF'
        else if (object > 10) then
          write(fstep,'(A5,I2,A4)')   'FSTEP',OBJECT,'.LBF'
        end if

*********************CHANGES MADE HERE*****************************
        INQUIRE(FILE=fstep,EXIST=EXIST)
        IF (EXIST) THEN
        OPEN(UNIT=8,FILE=fstep,STATUS="UNKNOWN")
        READ(8,*)BRACKT
        READ(8,*)STAGE1
        READ(8,*)FINIT
        READ(8,*)DGTEST
        READ(8,*)WIDTH
        READ(8,*)WIDTH1
        READ(8,*)DGINIT
        READ(8,*)INFOC
        READ(8,*)STX
        READ(8,*)FX
        READ(8,*)DGX
        READ(8,*)STY
        READ(8,*)FY
        READ(8,*)DGY
        READ(8,*)STMIN
        READ(8,*)STMAX
        CLOSE(UNIT=8)
        ENDIF
************************TILL HERE***********************************

      IF(INFO.EQ.-1) GO TO 45
      INFOC = 1
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. FLTOL .LT. ZERO .OR.
     &    GLTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. STPMIN .LT. ZERO
     &    .OR. STPMAX .LT. STPMIN .OR. MAXFEV .LE. 0) RETURN
C
C     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
C     AND CHECK THAT S IS A DESCENT DIRECTION.
C
      DGINIT = ZERO
      DO 10 J = 1, N
         DGINIT = DGINIT + G(J)*S(J)
   10    CONTINUE
      IF (DGINIT .GE. ZERO) then
         write(LP,15)
   15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
         RETURN
         ENDIF
C
C     INITIALIZE LOCAL VARIABLES.
C
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = FLTOL*DGINIT
      WIDTH = STPMAX - STPMIN
      WIDTH1 = WIDTH/P5

	DO 20 J = 1, N
         WA(J) = X(J)
   20    CONTINUE
C
C     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
C     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
C     THE INTERVAL OF UNCERTAINTY.
C     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
C
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
C
C     START OF ITERATION.
C
   30 CONTINUE
C
C        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
C        TO THE PRESENT INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
         ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
            END IF
C
C        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
C
         STP = MAX(STP,STPMIN)
         STP = MIN(STP,STPMAX)
C
C        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
C        STP BE THE LOWEST POINT OBTAINED SO FAR.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     &      .OR. NFEV .GE. MAXFEV-1 .OR. INFOC .EQ. 0
     &      .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) STP = STX

*********************CHANGES MADE HERE*****************************
        OPEN(UNIT=8,FILE=fstep,STATUS="UNKNOWN")
        WRITE(8,*)BRACKT," BRACKT"
        WRITE(8,*)STAGE1," STAGE1"
        WRITE(8,*)FINIT," FINIT"
        WRITE(8,*)DGTEST," DGTEST"
        WRITE(8,*)WIDTH," WIDTH"
        WRITE(8,*)WIDTH1," WIDTH1"
        WRITE(8,*)DGINIT," DGINIT"
        WRITE(8,*)INFOC, " INFOC"
	WRITE(8,*)STX, " STX"
	WRITE(8,*)FX, " FX"
	WRITE(8,*)DGX," DGX"
	WRITE(8,*)STY," STY"
	WRITE(8,*)FY, " FY"
	WRITE(8,*)DGY," DGY"
	WRITE(8,*)STMIN, " STMIN"
	WRITE(8,*)STMAX, " STMAX"
        CLOSE(UNIT=8)
************************TILL HERE***********************************

C
C        EVALUATE THE FUNCTION AND GRADIENT AT STP
C        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
C        We return to main program to obtain F and G.
C
         DO 40 J = 1, N
            X(J) = WA(J) + STP*S(J)
   40       CONTINUE
         INFO=-1
         RETURN
C
   45    INFO=0
         NFEV = NFEV + 1
         DG = ZERO
         DO 50 J = 1, N
            DG = DG + G(J)*S(J)
   50       CONTINUE
         FTEST1 = FINIT + STP*DGTEST
C
C        TEST FOR CONVERGENCE.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     &      .OR. INFOC .EQ. 0) INFO = 6
         IF (STP .EQ. STPMAX .AND.
     &       F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
         IF (STP .EQ. STPMIN .AND.
     &       (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
         IF (NFEV .GE. MAXFEV) INFO = 3
         IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GLTOL*(-DGINIT)) INFO = 1
C
C        CHECK FOR TERMINATION.
C
         IF (INFO .NE. 0) RETURN
C
C        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
C
         IF (STAGE1 .AND. F .LE. FTEST1 .AND.
     &       DG .GE. MIN(FLTOL,GLTOL)*DGINIT) STAGE1 = .FALSE.
C
C        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
C        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
C        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
C        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
C
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
C
C           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
C
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST

           DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
C
C           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,
     &                 BRACKT,STMIN,STMAX,INFOC)
C
C           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
C
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
         ELSE
C
C           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,
     &                 BRACKT,STMIN,STMAX,INFOC)
            END IF
C
C        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
C        INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1)
     &         STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
            END IF
C
C        END OF ITERATION.
C
         GO TO 30
C
C     LAST LINE OF SUBROUTINE MCSRCH.
C
      END
