C UTEP Electronic Structure Lab (2020)
C> @file electronic_geometry.f
C> @note YY. This subroutine is a draft. Not meant for production.
C> FOD force calculation and optimization 
C> SIC
C> create an empty file, FODLBF, to use LBFGS
C> @param[out] energy
C> @param xvec : FOD positions
C> @param gvec : FOD forces
C> @param msite(2) : This array holds numbers of FODs for spin up/down. 
         subroutine fod_opt(energy,fod_converge)
!        subroutine electronic_geometry(energy)
         use global_inputs,only : fod_opt1,fod_opt2,fod_opt3,
     &                            symmetrymodule1
         use xmol,only : AU2ANG,NUM_ATMS,XMOL_LIST,GET_LETTER
         implicit real*8 (a-h,o-z)
         !INCLUDE  'PARAMA2' not yet.
         logical exist,reset,exist1
         logical fod_converge
         parameter(mxidt=500)
         dimension msite(2),nidt(2)
         dimension r(3,1000),f(3,1000),d2inv(1000),suggest(1000)
         dimension xvec(3000),gvec(3000)
         dimension scrv(18000)
         dimension nwht(mxidt*2),nequiv(120,mxidt*2)

C LBFGS
         double precision :: DGUESS,ACCSOLN,XTOL
         dimension :: IPRINT(2)
         integer :: MUPDATE
         dimension :: DIAG(3000) !DIAG(3*MXATMS) 
         !DIMENSION WORK(3*MXATMS*(2*NMUPMAX+1)+2*NMUPMAX)
         !dimension :: work(27008)   ! 3*1000* (2*4+1)+2*4 = 3000*9+8
         real*8, allocatable :: work(:)
         real*8 :: fmax,dist
         logical :: LFODLBF,scaled
         real(8) :: len1, len2, len3                ! lengths of forces etc.
         real(8) :: len4, len5, len6                ! evaluation of energy convergence -> reset LBFGS or not
         integer :: constraint,i,ig
         integer :: frozen(1000) ! 0: frozen 1: not frozen
         integer             :: n_atoms             ! number of atoms
         integer :: counter, count_hess          
         real(8),allocatable :: r_atoms(:,:)        ! atomic coodinates
         real(8), allocatable :: hess_D(:)          ! HESS_D used in the scaledLBFGS                     (former HESS_D)
         character(len=2),allocatable :: species(:) ! atomic species
         character(2)  :: LETTER
         allocate(r_atoms(NUM_ATMS,3))
         allocate(species(NUM_ATMS))


         do i=1,1000
           d2inv  (i)=1.0d0
           suggest(i)=1.0d0
         end do
         inquire(file='FRMIDT',exist=exist)
         if(exist)call symfrm(2) 
         open(90,file='temp_energy')
         write(90,*)energy
         close(90)
         call system('echo " "       >> records')
         call system('cat temp_energy >>records')
         call system('rm temp_energy')
         call system('cat FRMORB     >> records')
         call system('cat fforce.dat >> records')

         if(.not.exist)then
          open(90,file='FRMORB')
          open(91,file='fforce.dat')
          if(symmetrymodule1) then
            read(90,*)(nidt(i), msite(i),  i=1,2)
          else !legacy format
            read(90,*)msite(1),msite(2)
          endif
          nopt=0
c
c  KAJ 6-22-2023 -- Do we need FRMORB to have info about FOD symmetry?
c
          do if=1,msite(1)+msite(2)
!           read(90,*)(r(j,if),j=1,3)
            read(90,*)(r(j,if),j=1,3),nwht(if),
     &          (nequiv(i,if),i=1,   nwht(if))
           read(91,*)(f(j,if),j=1,3)
          end do
          close(90)
          close(91)
         else
          open(90,file='gforce.dat')
          read(90,*)msite(1),msite(2)
          nopt=0
          do if=1,msite(1)+msite(2)
           read(90,*)(r(j,if),j=1,3)
           read(90,*)(f(j,if),j=1,3)
          end do
          close(90)
         end if

         scaled     = .false.
         call check_inputs 
         LFODLBF=fod_opt1
          scaled = fod_opt2
          constraint = fod_opt3 
          do ig=1,msite(1)+msite(2)
           do i=1,NUM_ATMS
             CALL GET_LETTER(XMOL_LIST(i)%ANUM,LETTER)
             species(i)=LETTER
             r_atoms(i,1)=XMOL_LIST(i)%RX
             r_atoms(i,2)=XMOL_LIST(i)%Ry
             r_atoms(i,3)=XMOL_LIST(i)%Rz
             dist=(r(1,ig)-r_atoms(i,1))**2.0d0+(r(2,ig)-r_atoms(i,2))
     &       **2.0d0 +(r(3,ig)-r_atoms(i,3))**2.0d0
             dist=sqrt(dist)
             frozen(ig)=1
             if (constraint .eq. 1) then
             if (trim(adjustl(species(i))).ne.'H') then
             if  (dist .lt. 0.1d0 ) then
             frozen(ig)=0
             GO TO 227
             else 
             frozen(ig)=1
             end if
             end if
             end if
           end do
 227           CONTINUE 
          end do
!          open(3000,file='frozen.tmp',status='old')
!          rewind(3000)
!          do if=1,msite(1)+msite(2)
!           read(3000,*) frozen(if)
!          end do
!          close(3000,status='delete')

         gabs=0.0d0
         fmax=0.0d0
!!PB changed
           nopt=0
         do ig=1,msite(1)+msite(2)
          if (frozen(ig).ne. 0) then
           do j=1,3
             gabs=gabs+f(j,ig)**2
           end do
           fmax=max(fmax,f(1,ig)**2+f(2,ig)**2+f(3,ig)**2)
           do j=1,3
             nopt=nopt+1
             xvec(nopt)= r(j,ig)
             gvec(nopt)= f(j,ig) *dble(frozen(ig))
           end do
          end if 
         end do
         gabs1=gabs/nopt
         gabs1=sqrt(gabs1)
         gabs=sqrt(gabs)
         fmax=sqrt(fmax)
         open(80,file='fande.dat')
         nnn=0
         read(80,*,end=10)nnn
 10      continue
         rewind(80)
         nnn=nnn+1
         write(80,80)nnn,energy,gabs,fmax!,qmin,qmax,detr
         close(80)
         call system('cat fande.dat >>fande.out')
 80      format(i5,f20.12,4g20.12,f20.12)
         gtol=0.000001d0
         ftol=0.000001d0
         mopt=1000 !1000 => 3000

           ! KAJ: check for convergence -- is the max force on an FOD less than
           ! the tolerance?
           !
           if(fmax.lt.0.0001) then       !set force tolerance here
             fod_converge=.true.
            if (fod_opt1) then
           open(99,file='FDIAG.LBF',form='formatted',status='unknown')
           close(99,status='delete')
           open(99,file='FSEARCH.LBF',form='formatted',status='unknown')
           close(99,status='delete')
           open(99,file='FSTEP.LBF',form='formatted',status='unknown') 
           close(99,status='delete')
             print *, 'Converged in ELECTRONIC GEOMETRY deleting LBF
     &         files'
              else
           open(99,file='FGRAD',form='formatted',status='unknown')
           close(99,status='delete')
           open(99,file='FGRLOG',form='formatted',status='unknown')
           close(99,status='delete')
            end if
             return
           end if
!END OF COPIED LINES
!YY. In case you don't want to optimized FODs, change this to .true.
         if(.false.) return

!         inquire(file='FODLBF',exist=exist)
!         if (.not. exist) THEN
!          open(unit=1199,file='FODLBF',status='new')
!          write(1199,*) "F    0.1"
!         end if
         DGUESS = 0.1d0
!         open(unit=1199,file='FODLBF')
!         rewind(1199)
!         read(1199,*) LFODLBF,DGUESS
!         close(1199)
!PB changed for scaled LBFGS

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
           allocate(hess_D(3*(msite(1)+msite(2))))
           !
           ! open HESS_DIAG file
           !
           open(95,file='HESS_DIAG',status='old',action='read')
           !
           ! Write the corresponding array
           !
           count_hess = 0
           do ig=1,msite(1)+msite(2)
             do j = 1, 3
               count_hess = count_hess + 1
               read(95,*) hess_D(count_hess)
             end do
           end do
           close(95)
         end if
        
        if (scaled) then
         nopt=0
         count_hess = 0
         do ig=1,msite(1)+msite(2)
          do j = 1, 3
           count_hess = count_hess + 1
           if (frozen(ig).ne. 0) then
            nopt=nopt+1
            xvec(nopt)= r(j,ig)
            gvec(nopt)= f(j,ig)
            xvec(nopt)=xvec(nopt)/hess_D(count_hess)
            gvec(nopt)=gvec(nopt)*hess_D(count_hess)
            r(j,ig)=xvec(nopt)
            f(j,ig)=gvec(nopt)
            end if
           end do
          end do  
        end if
!END of scaled LBFGS part
         if(LFODLBF) then
          print * ,'LFODLBF',LFODLBF


c !YY LBFGS parameters
          !DGUESS=0.03d0
          !Default value used in DFT update is 0.1d0
          !Large -> More jumps in updating. 
          !Small -> slow convergence

          IPRINT(1)=-1
          IPRINT(2)=0
          MUPDATE=4
          ACCSOLN=1.0D-6
          XTOL=10D-16
          IFLAG=0
c !Preparation
          GMAX= 0.0d0
          GSUM= 0.0d0
!          NPAR= 3*(msite(1)+msite(2))
          ! 3(xyz)*(no. FOD UP + no. FOD DN)
          allocate(work(27008))
c Read
          inquire(file='FDIAG.LBF',exist=exist)
          if(exist) then
           open(unit=4,file='FDIAG.LBF',status='old')
           read(4,*)iflag
           read(4,*)
           do i=1,nopt
            read(4,*)diag(i)
           end do
           read(4,*)
           do i=1,nopt*(2*MUPDATE+1)+2*MUPDATE
            read(4,*)work(i)
           end do
           close(4)
          end if
!PB modification
         counter = 0
         OPEN (90, file = 'fande.out')
         DO
             READ (90,*, END=20)
             counter = counter + 1
         END DO
 20      CLOSE (90)

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
         if (reset) then
             write(6,*) 'LBFGS is stuck'
             iflag=0
         end if

          call folbfgs(nopt,mupdate,xvec(1:nopt),energy,gvec(1:nopt),
     &         .false.,diag,iprint,accsoln,xtol,work,iflag,dguess)
! xvec should be the updated FOD position
          if(iflag .eq. 0) then
           open(99,file='FDIAG.LBF',form='formatted',status='unknown')
           close(99,status='delete')
           open(99,file='FSEARCH.LBF',form='formatted',status='unknown')
           close(99,status='delete')
           open(99,file='FSTEP.LBF',form='formatted',status='unknown') 
           close(99,status='delete')
          else if (iflag .GT. 0) then
           open(unit=4,file='FDIAG.LBF',status='unknown') 
           write(4,*)iflag, "  IFLAG"
           write(4,*)"=====DIAG===="
           do i=1,nopt
             write(4,*)DIAG(I)
           END DO
           write(4,*)"====WORK===="
           do I=1,NOPT*(2*MUPDATE+1)+2*MUPDATE
            write(4,*)WORK(I)
           end do
           close(4)
          end if     

          deallocate(work)

         else
          call fodcgrad(nopt,mopt,energy,xvec,gvec,gtol,ftol,scrv,istat) 
         print *,'CONJUGATE GRADIENT'
         end if

!PB changed to rescale r and f 
         nopt=0
         count_hess = 0
         if (scaled) then
          do ig=1,msite(1)+msite(2)
           do j = 1, 3
            count_hess = count_hess + 1
            if (frozen(ig).ne. 0) then
             nopt=nopt+1
             xvec(nopt)=xvec(nopt)*hess_D(count_hess)
             gvec(nopt)=gvec(nopt)/hess_D(count_hess)
             r(j,ig)=xvec(nopt)
             f(j,ig)=gvec(nopt)
            end if
           end do
          end do
         end if  

c        print*,'istat:',istat
         inquire(file='FRMIDT',exist=exist) !YY recheck if FRMIDT is there 
         if(.not.exist)then
           open(90,file='FRMORB')
           rewind(90)
           if(symmetrymodule1) then
             write(90,*) (nidt(i),msite(i),i=1,2)
           else !legacy format
             write(90,*)msite(1),msite(2)
           endif
         else
           open(90,file='FRMIDT')
           rewind(90)
           write(90,*)msite(1),msite(2)
         end if

         nopt=0
         do ig=1,msite(1)+msite(2)
          if (frozen(ig).ne.0) then
           do j=1,3
             nopt=nopt+1
             r(j,ig)=xvec(nopt)
             f(j,ig)=gvec(nopt)
           end do
c           print*,'frozenif',frozen(if)
          end if
          if(.not.exist)then !FRMORB
c           write(90,*) (r(j,ig),j=1,3)
           write(90,101)(r(j,ig),j=1,3),nwht(ig),
     &             (nequiv(i,ig),i=1,   nwht(ig))
          else !FRMIDT
           write(90,90)(r(j,ig),j=1,3),int(frozen(ig)),
     &                d2inv(ig),suggest(ig)
           !print*,"New FOD position",(r(j,if),j=1,3)
          end if
         end do
         close(90)
c        call system('echo "0 1" >  RUNS')
c        call system('echo "4 4" >> RUNS')
c        call system('echo "0  " >> RUNS')
         !call stopit
 90      format(3d20.12,I3,'  D2INV=',F12.4,' SUGGEST= ',F12.4)
101      format(' ',3F20.10,' ',I4,120I4)


         end

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
     & IFLAG,DGUESS)
C
      INTEGER N,M,IPRINT(2),IFLAG
      DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
      DOUBLE PRECISION F,EPSL,XTOL,DGUESS
      LOGICAL DIAGCO
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
**********************CHANGES MADE HERE*********************
C IF GRADIENTS HAVE CONVERGED THEN EXIT
        IF (N .LT. 1) THEN
         GOTO 300
        ENDIF

        INQUIRE(FILE='FSEARCH.LBF',EXIST=EXIST)
        IF (EXIST) THEN
         OPEN(UNIT=7,FILE='FSEARCH.LBF',STATUS='OLD')
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
     &            XTOL,MAXFEV,INFO,NFEV,DIAG)
*****************CHANGES MADE HERE*****************
        OPEN(UNIT=7,FILE='FSEARCH.LBF',STATUS='UNKNOWN')
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
        OPEN(UNIT=4,FILE='FDIAG.LBF',FORM='formatted',STATUS='unknown')
        CLOSE(4,STATUS='delete')
        OPEN(UNIT=7,FILE='FSEARCH.LBF',FORM='formatted',
     &                                 STATUS='unknown')
        CLOSE(7,STATUS='DELETE')
        OPEN(UNIT=8,FILE='FSTEP.LBF', FORM='formatted',STATUS='unknown')
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
      SUBROUTINE FOMCSRCH(N,X,F,G,S,STP,FLTOL,XTOL,MAXFEV,INFO,NFEV,WA)
      INTEGER N,MAXFEV,INFO,NFEV
      DOUBLE PRECISION F,STP,FLTOL,GLTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
      COMMON /LB3/MP,LP,GLTOL,STPMIN,STPMAX
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

*********************CHANGES MADE HERE*****************************
        INQUIRE(FILE="FSTEP.LBF",EXIST=EXIST)
        IF (EXIST) THEN
        OPEN(UNIT=8,FILE="FSTEP.LBF",STATUS="UNKNOWN")
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
        OPEN(UNIT=8,FILE="FSTEP.LBF",STATUS="UNKNOWN")
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
