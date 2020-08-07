C UTEP Electronic Structure Lab (2020)
C> @file electronic_geometry.f

C> FOD force calculation and optimization 
C> SIC
C> @param[out] energy
         subroutine electronic_geometry(energy)
         implicit real*8 (a-h,o-z)
         logical exist
         dimension msite(2)
         dimension r(3,1000),f(3,1000),d2inv(1000),suggest(1000)
         dimension xvec(3000),gvec(3000)
         dimension scrv(18000)

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
          read(90,*)msite(1),msite(2)
          nopt=0
          do if=1,msite(1)+msite(2)
           read(90,*)(r(j,if),j=1,3)
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
         gabs=0.0d0
         fmax=0.0d0
         do if=1,msite(1)+msite(2)
           do j=1,3
             gabs=gabs+f(j,if)**2
           end do
           fmax=max(fmax,f(1,if)**2+f(2,if)**2+f(3,if)**2)
c          print*,(f(j,if),j=1,3)
           do j=1,3
             nopt=nopt+1
             xvec(nopt)= r(j,if)
             gvec(nopt)= f(j,if)
           end do
         end do
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
c        print*,'nopt:',nopt
         mopt=1000

         call fodcgrad(nopt,mopt,energy,xvec,gvec,gtol,ftol,scrv,istat) 

c        print*,'istat:',istat
         inquire(file='FRMIDT',exist=exist) !YY recheck if FRMIDT is there 
         if(.not.exist)then
           open(90,file='FRMORB')
         else
           open(90,file='FRMIDT')
         end if
         rewind(90)
         write(90,*)msite(1),msite(2)
         nopt=0
         do if=1,msite(1)+msite(2)
           do j=1,3
             nopt=nopt+1
             r(j,if)=xvec(nopt)
           end do
           write(90,90)(r(j,if),j=1,3),d2inv(if),suggest(if)
         end do
         close(90)
c        call system('echo "0 1" >  RUNS')
c        call system('echo "4 4" >> RUNS')
c        call system('echo "0  " >> RUNS')
         !call stopit
 90      format(3d20.12,' D2INV=',F12.4,' SUGGEST= ',F12.4)
         end
