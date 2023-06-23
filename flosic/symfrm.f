C UTEP Electronic Structure Lab (2020)
C> @file symfrm.f
C> A subroutine that read FRMIDT and FRMGRP
C> Also used in electronic_geometry
        subroutine symfrm(mode)
        use global_inputs,only : symmetrymodule1
        implicit real*8 (a-h,o-z)
        character*100 line
        logical exist
        integer froze,ioerr
        parameter(mxidt=500)
        dimension nidt(2),nfrm(2)
        dimension g(3,3,121),fodidt(3,mxidt,2),fordrv(3,mxidt,2)
        dimension r(3),f(3),h(3)
        dimension nwht(mxidt,2),nequiv(120,mxidt,2)
        inquire(file='FRMIDT',exist=exist)
        if(.not.exist)return
        inquire(file='FRMGRP',exist=exist)
        if(symmetrymodule1) then
          call system('cp GRPMAT FRMGRP')
        elseif(.not.exist) then
          !make FRMGRP identity matrix
          open(40,file='FRMGRP',status='new')
            write(40,*) '1'
            write(40,*) '1  0  0'
            write(40,*) '0  1  0'
            write(40,*) '0  0  1'
          close(40)
        else
          continue ! use user provided FRMGRP
        end if
        inquire(file='frozen.tmp',exist=exist)
        if(exist)call system('rm  frozen.tmp') 
!       open(3000,file='frozen.tmp',status='new')
        open(300,file='FRMIDT')
        read(300,*)nidt(1),nidt(2)
        do ispn=1,2   
          do idnt=1,nidt(ispn)
            froze=1
            read(300,50)line
 50         format(a100)
! Freeze/not freeze FOD position
!           read(line,*,IOSTAT=ioerr)(fodidt(j,idnt,ispn),j=1,3),froze
! Previous statement:
            read(line,*,ERR=51)(fodidt(j,idnt,ispn),j=1,3)
 51         continue
            print  103,(fodidt(j,idnt,ispn),j=1,3),froze
 103        format(' ',3F20.10,I3)
!           write(3000,*) froze
          end do
        end do
        close(300)
!       close(3000)
        open(300,file='FRMGRP')
        read(300,*)ngp
        do igp=1,ngp
          read(300,*)((g(j,i,igp),j=1,3),i=1,3)
        end do
c close the group:
        do icls=1,5
          mgp=ngp
          do igp=1,ngp
            do jgp=1,ngp
              mgp=mgp+1
              if(mgp.gt.121)stop 'symfrm mgp'
              do i=1,3
                do j=1,3
                  g(j,i,mgp)=0.0d0
                  do k=1,3
                    g(j,i,mgp)=g(j,i,mgp)+g(j,k,igp)*g(k,i,jgp)
                  end do
                end do
              end do
              errmin=1.0d30
              do kgp=1,mgp-1
                err=0.0d0
                do i=1,3
                  do j=1,3
                    err=err+abs(g(j,i,kgp)-g(j,i,mgp))
                  end do
                end do
                errmin=min(err,errmin)
              end do
              if(errmin.lt.1.0d-5)mgp=mgp-1
            end do
          end do
          print*,icls,ngp,mgp
          ngp=mgp
        end do
        rewind(300)
        write(300,*)ngp
        do igp=1,ngp
          write(300,100)((g(j,i,igp),j=1,3),i=1,3)
          write(300,*)' '
        end do
 100    format(' ',3F20.10)
 101    format(' ',3F20.10,' ',I4,120I4)
 102    format( 'H ',4F20.10)
        close(300)
c end of group closure
c generate all inequivalent Fermi Orbital Discriptors:
        open(300,file='FRMORB')
        if(mode.eq.1)then
          do ispn=1,2    
            nfrm(ispn)=nidt(ispn)
            do idnt=1,nidt(ispn)  
             nwht(idnt,ispn)=1
             nequiv(nwht(idnt,ispn),idnt,ispn)=idnt       
             do igp=1,ngp
              do j=1,3
               fodidt(j,nfrm(ispn)+1,ispn)=0.0d0
               do k=1,3
                fodidt(j,nfrm(ispn)+1,ispn)=
     &                   fodidt(j,nfrm(ispn)+1,ispn)+
     &                   g(j,k,igp)*fodidt(k,idnt,ispn)
               end do
              end do
              errmin=1.0d30
              do ifrm=1,nfrm(ispn)
               err=0.0d0
               do j=1,3
                err=err+abs(fodidt(j,nfrm(ispn)+1,ispn)
     &                     -fodidt(j,ifrm,ispn))
               end do
               errmin=min(err,errmin)
              end do
              if(errmin.gt.0.001)then
               nfrm(ispn)=nfrm(ispn)+1
               nwht(idnt,ispn)=nwht(idnt,ispn)+1
               print*,'nfrm(ispn):',nfrm(ispn)
               nequiv(nwht(idnt,ispn),idnt,ispn)=nfrm(ispn)
              end if
             end do
            end do
          end do
          if(symmetrymodule1) then
            write(300,20)(nidt(ispn),nfrm(ispn),ispn=1,2)
          else
            write(300,20)nfrm,nidt
          endif
          do ispn=1,2          
            if(ispn.eq.1)open(99,file='XMOL_FRM.UP')
            if(ispn.eq.2)open(99,file='XMOL_FRM.DN')
            write(99,*)nfrm(ispn)
            write(99,*)nfrm,nidt
            do ifrm=1,nfrm(ispn)
             if(ifrm.le.nidt(ispn))then
              do jfrm=1,nwht(ifrm,ispn)
               kfrm=nequiv(jfrm,ifrm,ispn)
               rn=fodidt(1,ifrm,ispn)**2+
     &            fodidt(2,ifrm,ispn)**2+
     &            fodidt(3,ifrm,ispn)**2
               rn=sqrt(rn)
               write(99,102)(fodidt(j,kfrm,ispn)*0.529,j=1,3),rn
              end do
              write(300,101)(fodidt(j,ifrm,ispn),j=1,3),nwht(ifrm,ispn),
     &                      (nequiv(i,ifrm,ispn),i=1,nwht(ifrm,ispn))
             else
              if(ifrm.eq.nidt(ispn)+1)write(300,*)' '
              write(300,100)(fodidt(j,ifrm,ispn),j=1,3)
             end if
            end do
            write(300,*)' '
            write(300,*)' '
            close(99)
            if(ispn.eq.1)call system('cat XMOL_FRM.UP >> movie.up')
            if(ispn.eq.2)call system('cat XMOL_FRM.DN >> movie.dn')
          end do
          close(300)
!if mode is .ne. 1
        else 
          open(300,file='FRMORB')
          open(301,file='fforce.dat')
          open(302,file='gforce.dat')
          if(symmetrymodule1) then
            read(300,*)(nidt(ispn),nfrm(ispn),ispn=1,2)
          else !legacy format
            read(300,*)nfrm,nidt
          endif
          write(302,*)nidt
          do ispn=1,2
            do ifrm=1,nfrm(ispn)
              if(ifrm.le.nidt(ispn))then
                read(300,*)(fodidt(j,ifrm,ispn),j=1,3),
     &                      nwht(ifrm,ispn),
     &          (nequiv(i,ifrm,ispn),i=1,nwht(ifrm,ispn))
              else
                read(300,*)(fodidt(j,ifrm,ispn),j=1,3) 
              end if
              read(301,*)(fordrv(j,ifrm,ispn),j=1,3)
            end do
            do ifrm=1,nfrm(ispn)
              print *,ifrm,(fodidt(j,ifrm,ispn),j=1,3)
            end do
            do ifrm=1,nfrm(ispn)
              print *,ifrm,(fordrv(j,ifrm,ispn),j=1,3)
            end do
c   symmetrize the force:
            do idnt=1,nidt(ispn)
              do j=1,3
                h(j)=0.0d0
              end do
              fnorm=fodidt(1,idnt,ispn)**2+
     &              fodidt(2,idnt,ispn)**2+
     &              fodidt(3,idnt,ispn)**2 
              fnorm=sqrt(fnorm)
              nreplica=0
              do iwht=1,nwht(idnt,ispn)         
                kwht=nequiv(iwht,idnt,ispn)
                do igp=1,ngp
                  err=0.0d0
                  do j=1,3
                    r(j)=0.0d0
                    f(j)=0.0d0
                    do k=1,3
                      r(j)=r(j)+g(k,j,igp)*fodidt(k,kwht,ispn)
                      f(j)=f(j)+g(k,j,igp)*fordrv(k,kwht,ispn)
                    end do
                    err=err+abs(r(j)-fodidt(j,idnt,ispn))
                  end do
                  if(err.le.0.001*fnorm)then
                    nreplica=nreplica+1
                    do j=1,3
                      h(j)=h(j)+f(j)
                    end do
                    print 30,h,r,f
 30                 format(9f12.6)
                  end if
                end do
              end do
              print*,idnt,nreplica,h/nreplica
              do j=1,3
                fordrv(j,idnt,ispn)=h(j)/nreplica
              end do
              do j=1,3
                fordrv(j,idnt,ispn)=fordrv(j,idnt,ispn)
     &                             *nwht(idnt,ispn)
              end do
              write(302,101)(fodidt(j,idnt,ispn),j=1,3),
     &                       nwht(idnt,ispn),  
     &         (nequiv(i,idnt,ispn),i=1,nwht(idnt,ispn))
              write(302,100)(fordrv(j,idnt,ispn),j=1,3)
            end do
            write(302,*)
          end do
          close(300)
          close(301)
          close(302)
        end if
 20     format(' ',2I5,'    ',2I5)
        return
        end
                     
                     

                    
  

         
       
 
        
