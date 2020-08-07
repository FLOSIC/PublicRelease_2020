! UTEP Electronic Structure Lab (2020)
 	implicit none	
	include 'mpif.h'
        integer :: Natom, Nmsh, Npt
        integer :: Ndipole, Nx, Ny, Nz, Mdipole
        integer :: nnsolute, iz, Mdip, niter
	integer:: i,j,k,iatom,jatom,M,p,iter, T
        real*8 :: sum, disti,distj,add,func,qtot
        real*8 :: funcold,factorlo, factorhi
        real*8 :: lamda,dx,dy,dz,dipole2
        real*8 :: epsl,xtol,dguess,gtol,mfrc,ftol
        real*8 :: dipolex,dipoley,dipolez
        real*8 :: comx,comy,comz
        real*8 :: xmin,xmax,ymin,ymax,zmin,zmax
        real*8 :: xbmin,xbmax,ybmin,ybmax,zbmin,zbmax
        real*8 :: xx,yy,zz
        real*8 :: spacing, distnn, dist, K_B, kT
        real*8 :: bufferx,buffery,bufferz,dsolvent
        real*8 :: ranx,error,anum
        real*8 :: pi, twopi, theta,phi
        real*8 :: old_mux,old_muy,old_muz
        real*8 :: Ene, Enew, Eold, Ebest, Eave
        real*8 :: deltaE, Boltz_prob
        real*8 :: rijx,rijy,rijz,pidotrij, pjdotrij, pidotpj
	integer, dimension (:), allocatable :: Mpt
	integer, dimension (2):: Iprint
	integer, dimension (:), allocatable  :: Iwork
	real*8, dimension (:), allocatable :: Solpot
	real*8, dimension (:), allocatable :: work, diag
	real*8, dimension (30) :: vdw
	real*8, dimension (3) :: dip
        INTEGER              :: isize,idate(8)
        INTEGER,ALLOCATABLE  :: iseed(:)
	character(len=30) :: filename
	integer :: rank,size,ierr,best_rank,step,istep
	integer, dimension (:), allocatable  :: node_list
	real*8  :: Global_Ebest
	character :: rank_txt*10
	character :: temp_txt*20
	character :: iter_txt*20
	character :: step_txt*20
	integer :: status

        Type grid
        real*8  x,y,z,wmsh
        end type grid


        Type ion
           real*8 Rx,Ry,Rz,charge
           integer Anum
        End type ion

        Type dgrid
        real*8  x,y,z
        real*8  mux,muy,muz
        end type dgrid


        Type (ion),dimension(:),allocatable :: atom
        Type (dgrid),dimension(:),allocatable :: dipole, best_config,ave_config
        Type (grid),dimension(:),allocatable :: Rmsh
        Logical Exist
        data vdw(1:30)/1.09,4*0.0d0,1.70,1.55,1.52,0.0d0,0.0d0,2.27,4*0.0d0,1.80,1.75,12*0.0d0,1.39/
        vdw= vdw/0.5292
        pi=4.0d0*datan(1.0d0)
        Twopi=2.0d0*pi
	K_B=8.617343E-5
	call MPI_Init(ierr)
	call MPI_Comm_size(MY_NEW_COM,size,ierr)
	call MPI_Comm_rank(MY_NEW_COM,rank,ierr)
               CALL DATE_AND_TIME(VALUES=idate)
               CALL RANDOM_SEED(SIZE=isize)
               ALLOCATE( iseed(isize) )
               CALL RANDOM_SEED(GET=iseed)
               iseed = iseed * (idate(8)-500)      ! idate(8) contains millisecond
               CALL RANDOM_SEED(PUT=iseed)
	allocate(node_list(size))
!   Read the  input parameters
!
!        Spacing : water grid spacing in Angstrom - cubic grid 
!        bufferx, buffery,bufferz in Angstrom: water layer thickness along x,y,z
!        dsolvent in Debye: permanent dipole moment of water
!        T : Temperature in Kelvin: 
!        Niter : number of iterations
!	step : number of steps that minimum energy is updated to every node (use -1 to desactivate)
!
	if(rank.eq.0) then
         open(5,file = 'input')

         read(5,*)spacing
         read(5,*)bufferx,buffery, bufferz
         read(5,*)dsolvent
         read(5,*)T
         read(5,*)niter
	 read(5,*)step
	 close(5)
	endif
	call MPI_Bcast(spacing,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
	call MPI_Bcast(bufferx,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
	call MPI_Bcast(buffery,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
	call MPI_Bcast(bufferz,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
	call MPI_Bcast(dsolvent,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
	call MPI_Bcast(kT,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
	call MPI_Bcast(T,1,MPI_INTEGER,0,MY_NEW_COM,ierr)
	call MPI_Bcast(niter,1,MPI_INTEGER,0,MY_NEW_COM,ierr)
	call MPI_Bcast(step,1,MPI_INTEGER,0,MY_NEW_COM,ierr)

	kT=K_B*T

!       convert to atomic units

        spacing=spacing/0.5292
        bufferx=bufferx/0.5292
        buffery=bufferx/0.5292
        bufferz=bufferx/0.5292
        dsolvent=dsolvent/2.5416d0


!       Create the water grid
        Nx=5
        Ny=1
        Nz=1
        Mdipole=Nx*Ny*Nz
        allocate(dipole(Mdipole))
        allocate(best_config(Mdipole))
        allocate(ave_config(Mdipole))
        write(6,*) 'Nx, Ny, Nz = ', Nx, Ny, Nz
        Ndipole=0

        do 120 i=1, Nx
             xx=(i-1)*spacing
             yy=0d0
	     zz=0d0
             Ndipole=Ndipole+1
             dipole(Ndipole)%x=xx
             dipole(Ndipole)%y=yy
             dipole(Ndipole)%z=zz

 120     continue

 50          format('Na',3F14.6)
 60          format('O',3F14.6)
!    	Get the initial structure.
       
         add=0.0d0
         do i = 1, Ndipole
           call random_number(ranx)
           theta  =ranx *pi
           call random_number(ranx)
           phi    =ranx*Twopi 
           dipole(i)%mux= dsolvent*sin(theta)*cos(phi)
           dipole(i)%muy= dsolvent*sin(theta)*sin(phi)
           dipole(i)%muz= dsolvent*cos(theta)
           dist=dipole(i)%mux**2+dipole(i)%muy**2+dipole(i)%muz**2
           dist=dsqrt(dist)
           add=add+dist
         end do
         Error=add-Ndipole*dsolvent
	write(6,*)'Error = ',Error,'Process = ',rank 

	write(rank_txt,70) rank
70	format(I0)
	rank_txt=adjustl(rank_txt)
	filename='watergrid-'//trim(rank_txt)//'.xyz'

        open(7, file = filename)
             write(7,*) Ndipole+Natom
             write(7,*) 
             do iatom=1, Ndipole
               write(7,400)  dipole(iatom)%x, dipole(iatom)%y, dipole(iatom)%z, &
	dipole(iatom)%mux,dipole(iatom)%muy,dipole(iatom)%muz
             end do
          
	close(7)

 
         call Energy(Natom, Ndipole,atom, dipole, Ene)
	 write(6,*) 'Energy = ',Ene,'process = ',rank
         Eold=Ene
         Ebest=Eold
         Eave=0.0d0
	 istep=1
         Do iter=1, Niter
           call random_number(ranx)
           Mdip = Ndipole*ranx+1
              call random_number(ranx)
              theta  =ranx *pi
              call random_number(ranx)
              phi    =ranx*Twopi 
              old_mux=dipole(Mdip)%mux
              old_muy=dipole(Mdip)%muy
              old_muz=dipole(Mdip)%muz
              dipole(Mdip)%mux= dsolvent*sin(theta)*cos(phi)
              dipole(Mdip)%muy= dsolvent*sin(theta)*sin(phi)
              dipole(Mdip)%muz= dsolvent*cos(theta)
              call Energy(Natom, Ndipole, atom, dipole, Ene)


         Enew=Ene
         deltaE= Enew-Eold
         if(Enew.lt. Eold) then
           Eold=Enew
         else
           Boltz_prob=exp(-deltaE/kT)
           call random_number(ranx)
           if ((Boltz_prob).gt.ranx) then
            Eold = Enew
           else 
           dipole(Mdip)%mux= old_mux
           dipole(Mdip)%muy= old_muy
           dipole(Mdip)%muz= old_muz
           end if 
         end if
	 write(6,300) iter, Mdip, Eold, Enew, Eave/iter,rank
         If (iter.gt.5000)then
            Eave=Eave+Eold
            do i=1,Ndipole
             ave_config(i)%mux=ave_config(i)%mux+dipole(i)%mux
             ave_config(i)%muy=ave_config(i)%muy+dipole(i)%muy
             ave_config(i)%muz=ave_config(i)%muz+dipole(i)%muz
            end do
         end if
         if (Ebest.gt.Eold)  then
            Ebest=Eold
            best_config=dipole
         end if
! Update minimum energy at given step
	 if(step.ne.-1)then
	  if(iter.eq.istep*step)then
	   do i=1,size
	    node_list(i)=0
	   end do
	   call MPI_Allreduce(Enew,Global_Ebest,1,MPI_DOUBLE_PRECISION,MPI_MIN,MY_NEW_COM,ierr)
	   if(Enew.eq.Global_Ebest)then
            node_list(1)=1
	   endif
	   call MPI_Allgather(node_list,1,MPI_INTEGER,node_list,1,MPI_INTEGER,MY_NEW_COM,ierr)
!	   write(6,*) 'node',rank,node_list(1),node_list(2),node_list(3),node_list(4),node_list(5),node_list(6)
	   do i=1,size
	    if(node_list(i).eq.1)exit
	   end do
	   write(step_txt,70)step*istep
	   step_txt=adjustl(step_txt)
	   write(temp_txt,70) T
	   temp_txt=adjustl(temp_txt)
	   filename='step-'//trim(step_txt)//'-'//trim(temp_txt)
	   if(rank+1.eq.i)then
            open(20, file=filename)
	    write(20,*)Enew
            do i=1, Ndipole
              write(20,400)dipole(i)%mux,dipole(i)%muy,dipole(i)%muz
            end do
            close(20)
	   endif
	   call MPI_Barrier(MY_NEW_COM,ierr)
	   if(rank.eq.0)then
            open(7,file=filename)
	      read(7,*)Ene
            do i=1,Ndipole
              read(7,*)dipole(i)%mux,dipole(i)%muy,dipole(i)%muz
            end do
            close(7)
	   endif
	   call MPI_Bcast(Ene,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
           do i=1,Ndipole
	    call MPI_Bcast(dipole(i)%mux,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
	    call MPI_Bcast(dipole(i)%muy,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
	    call MPI_Bcast(dipole(i)%muz,1,MPI_DOUBLE_PRECISION,0,MY_NEW_COM,ierr)
           end do
	   istep=istep+1
	  endif
	 endif

! end of update step
        end do
! end of iteration cycle
        call MPI_Barrier(MY_NEW_COM,ierr)
        write(6,*) 'Best energy =', Ebest,'process = ',rank
	write(rank_txt,70) rank
	rank_txt=adjustl(rank_txt)
	filename='best_config-'//trim(rank_txt)
        open(20, file=filename)
        write(20,*) Ndipole, Ebest
        do i=1, Ndipole
          write(20,400)best_config(i)%x,best_config(i)%y,best_config(i)%z,best_config(i)%mux,best_config(i)%muy,best_config(i)%muz
        end do
        close(20)
	call MPI_Allreduce(Ebest,Global_Ebest,1,MPI_DOUBLE_PRECISION,MPI_MIN,MY_NEW_COM,ierr)
	if(rank.eq.0)then
          write(6,*) 'Best Global Energy =', Global_Ebest
	endif
	if(Global_Ebest.eq.Ebest)then
	   write(rank_txt,70) rank
	   rank_txt=adjustl(rank_txt)
	   write(temp_txt,70) T
	   temp_txt=adjustl(temp_txt)
	   write(iter_txt,70) niter
	   iter_txt=adjustl(iter_txt)
	   filename='global-'//trim(iter_txt)//'-'//trim(temp_txt)//'-'//trim(rank_txt)
           open(20, file=filename)
           write(20,*) Ndipole, Ebest
           do i=1, Ndipole
             write(20,400)best_config(i)%x,best_config(i)%y,best_config(i)%z,best_config(i)%mux,&
	best_config(i)%muy,best_config(i)%muz
           end do
           close(20)
	endif
	write(rank_txt,70) rank
	rank_txt=adjustl(rank_txt)
	filename='current_config-'//trim(rank_txt)
         open(20, file=filename)
         write(20,*) Ndipole, Eold, Eave/Niter, Niter
         do i=1, Ndipole
           write(20,400)dipole(i)%x,dipole(i)%y,dipole(i)%z,dipole(i)%mux,dipole(i)%muy,dipole(i)%muz
         end do
         close(20)

	write(rank_txt,70) rank
	rank_txt=adjustl(rank_txt)
	filename='average_config-'//trim(rank_txt)
         open(20, file=filename)
         write(20,*) Ndipole, Eold, Eave/Niter, Niter
         do i=1, Ndipole
           write(20,400) dipole(i)%x, dipole(i)%y, dipole(i)%z, ave_config(i)%mux/Niter, &
         ave_config(i)%muy/Niter, ave_config(i)%muz/Niter
         end do
         close(20)

!        Write out the potential in the VMOLD grid
          
!         Solpot=0
!         Do i=1,Nmsh
!            Do j=1, Ndipole
!             rijx=dipole(j)%x - Rmsh(j)%x
!             rijy=dipole(j)%y - Rmsh(j)%y
!             rijz=dipole(j)%z - Rmsh(j)%z
!              dist=rijx**2+rijy**2+rijz**2
!              dist=dsqrt(dist)
!              anum=ave_config(j)%mux*rijx+ave_config(j)%muy*rijy+ave_config(j)%muz*rijz
!              solpot(i)=solpot(i)+anum/dist**3
!            end do
!          end do
              
         

 300     format(2I8, 3F18.6)
 400     format(6F14.6)
	call MPI_Finalize(ierr)
        stop
	end


          subroutine Energy(natm, Ndip, atm, dip, En)

          integer, intent(in) :: natm, ndip
          integer :: i,j
          real*8 :: rijx,rijy,rijz,dist
          real*8 :: pidotrij,pjdotrij,pidotpj
          real*8, intent(out) :: En
          type ion
           real*8  :: x,y,z,charge
           integer :: At 
          end type ion
          type grid
           real*8  :: x,y,z
           real*8  :: px,py,pz
          end type grid

          type(ion), dimension(natm),intent(in) ::  atm
          type(grid), dimension(ndip),intent(in) ::  dip


          En=0.0d0
          do i=1, Ndip
           do j=i+1,Ndip
           rijx=dip(i)%x - dip(j)%x
           rijy=dip(i)%y - dip(j)%y
           rijz=dip(i)%z - dip(j)%z
           dist=rijx**2+rijy**2+rijz**2
           dist=dsqrt(dist)
           pidotrij=dip(i)%px*rijx+dip(i)%py*rijy+dip(i)%pz*rijz
           pjdotrij=dip(j)%px*rijx+dip(j)%py*rijy+dip(j)%pz*rijz
           pidotpj =dip(i)%px*dip(j)%px+dip(i)%py*dip(j)%py+dip(i)%pz*dip(j)%pz
           if(dist.gt.1.0d-1) then
           En=En+pidotpj/dist**3 - 3.0d0*pidotrij*pjdotrij/dist**5
           end if
          end do 
         end do 

          return
          end

