! UTEP Electronic Structure Lab (2020)
!  This is a driver routine to call empirical disperion routines  
!  from Grimme's group. (DFT-D3)
!           Raja- El Paso, June, 2010

!        implicit none
        subroutine d_dftd3(max_ident,ident,ridt,edisp,ug)
        integer     :: i,j,inp,natoms,max_ident,ident,k
        integer, allocatable :: izat(:)
        real*8    :: aux,beriz,evdw,edisp,ridt(3,max_ident)
        real*8    :: fdisp(3), ug(3,max_ident)
        real*8    :: x,y,z,xp,yp,zp,zz,fdiff
        character*80    :: func
        real*8, allocatable :: xyz(:,:),g(:,:)
        logical   ::  exist
        zz = 0.0001d0
        INP = 99

          write(6,*) 'Welcome to D_DFTD3'
          call flush(6)
        INQUIRE(FILE='XMOL.DAT',EXIST=EXIST)
        if (EXIST) then
          write(6,*) 'Calculating vander Wall forces using DFT3D method'
          open(INP, file = 'XMOL.DAT')
        else 
          stop 'XMOL.DAT not found for DFT3D calculation'
        endif
          call flush(6)

         read(INP,*) natoms
            allocate(xyz(3,natoms))
            allocate(izat(natoms))
            allocate(g(3,natoms))
         read(INP,*) 
          do i=1,natoms
            read(INP,*) izat(i),xyz(1,i),xyz(2,i),xyz(3,i)
            xyz(:,i) = xyz(:,i)/0.529177249d0   ! Convert to atomic units.
          enddo

           func = 'pbe'  ! DFT+D3 parameters of PBE functionals will be used.
           edisp = 0.d0

          call   dftd3_grimme(func,natoms,izat,xyz,edisp,g)

          write(6,*) 'Dispersion energy = ',edisp     
          write(6,*) 'Writing forces from main'
          do i=1,natoms
!           write(6,10) izat(i),xyz(1,i),xyz(2,i),xyz(3,i)
!           write(6,11) g(1,i),g(2,i),g(3,i)
          enddo
!     ify symmetry unique atoms and return the forces on them
!   to add to Hellman-Feynman and Pulay forces.

               write(6,*) 'WRITING POSITIONS FROM DISP DRIVER'
          do i =1, ident
                 write(6,11) ridt(1,i), ridt(2,i), ridt(3,i)
          enddo
               call flush(6)

               k = 0;     
          do i =1, ident
                xp=ridt(1,i); yp=ridt(2,i); zp=ridt(3,i);
             do j=1,natoms
                 x=xyz(1,j); y=xyz(2,j); z=xyz(3,j);
                if ( abs(x-xp) < zz  .and.  abs(y-yp) < zz .and.   abs(z-zp) < zz) then
                  k = k+1;
                 !ug(1,i) = x; ug(2,i) = y; ug(3,i) = z;
                 ug(1,i) = g(1,j); ug(2,i) = g(2,j); ug(3,i) = g(3,j);
                endif
             enddo
           enddo
            write(6,*) 'No. of inequivalent atoms found',k

            do i=1,k
                 write(6,11) ug(1,i), ug(2,i), ug(3,i)
                fdisp(:) = ug(:,i)
            CALL FRCSYM(RIDT(1,i),Fdisp,FDIFF)
            IF (FDIFF .GT. 0.001) THEN
             PRINT 12,IID,FDIFF
 12   FORMAT('WARNING: Dispersion foRCE OF ATOM',I4,'VIOLATES SYMMETRY BY ',D12.4)
            ENDIF
            enddo
            call flush(6)
          
                      


 10     format(1x,i3,1x,3(f13.6,1x))
 11     format(3x,3(f13.6,1x))
         return
         end

          

