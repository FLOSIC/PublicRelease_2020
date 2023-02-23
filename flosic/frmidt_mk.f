         subroutine shell_analysis
         implicit real*8 (a-h,o-z)
         parameter (lm=4)
         parameter (mxd=(lm+1)**2)
         dimension g(3,3,121)
         dimension shell(6,1000),r(3),v(3),nshell(300)
         dimension angle(3,300),domega(300),rank(300)
         dimension nsoln(100,100)
         dimension nmem(100,2)
         dimension ylm(121,mxd)
         character*100 line
c        do LMAX=1,5
c        CALL BFSA(300,LMAX,NPTS,ANGLE,DOMEGA)
c        PRINT*,LMAX,NPTS
c             DO IPTS=1,NPTS
c             PRINT 30,(ANGLE(J,IPTS),J=1,3)
c             END DO  
c        END DO
         open(400,file='SHELL_ANALYSIS')
         call atomic_radii
         open(90,file='GRPMAT')
         read(90,*)ngp
         read(90,*)(((g(j,i,igp),j=1,3),i=1,3),igp=1,ngp)
         close(90)
         best=0.0d0
         print*,'Group Generators:' 
         do igp=1,ngp
              write(400,30)((g(j,i,igp),j=1,3),i=1,3)
 30           format(9F12.6)
         end do
         write(400,*)
             do itry=1,10
                mgp=ngp
                do igp=1,mgp
                do jgp=1,mgp
                  do m=1,3
                  do n=1,3
                  g(m,n,ngp+1)=0.0d0
                     do k=1,3
                     g(m,n,ngp+1)=g(m,n,ngp+1)+
     &                   g(m,k,igp)*g(k,n,jgp)
                     end do
                  end do
                  end do
                  err_min=1.0D30
                  do kgp=1,ngp
                  err=0.0d0
                     do m=1,3
                     do n=1,3
                     err=abs(g(m,n,kgp)-g(m,n,ngp+1))+err
                     end do
                     end do
                  err_min=min(err_min,err)
                  end do
                  if(err_min.gt.0.001)ngp=ngp+1
                end do
                end do
             end do
             write(400,*)'Order of Group:',mgp,ngp
         ns=0
         do itry=1,2
         do lmax=1,13
         if(itry.eq.1)then
         call angmsh(300,lmax,nang,angle,domega)
         else
c        call bfsa  (300,lmax,nang,angle,domega)
         end if
           do iang=1,nang
                 ns=ns+1        
                 shell(1,ns)=angle(1,iang)
                 shell(2,ns)=angle(2,iang)
                 shell(3,ns)=angle(3,iang)
                 shell(4,ns)=1.0d0
           idisp=0
           do js=1,ns-1
            do igp=1,ngp
              do k=1,3
              v(k)=0.0d0
               do l=1,3
               v(k)=v(k)+g(k,l,igp)*shell(l,js)
               end do
              end do
              err=abs(v(1)-shell(1,ns))+
     &            abs(v(2)-shell(2,ns))+
     &            abs(v(3)-shell(3,ns)) 
              if(err.lt.0.001)idisp=1
            end do
           end do
           if(idisp.eq.1)ns=ns-1
           end do
         end do
         end do
         do is=1,ns
         print*,is,(shell(j,is),j=1,3)
         end do
         do isx=6,-6,-1 
         do isy=6,-6,-1   
         do isz=6,-6,-1   
         ns=ns+1
         shell (1,ns)=float(isx)
         shell (2,ns)=float(isy)
         shell (3,ns)=float(isz)
         shell (4,ns)=1.0
         ss=shell(1,ns)**2+shell(2,ns)**2+shell(3,ns)**2
                      if(ss.gt.0.001)then
                      do j=1,3
                      shell(j,ns)=shell(j,ns)/sqrt(ss)
                      end do
                      end if
           idisp=0
           do js=1,ns-1
            do igp=1,ngp
              do k=1,3
              v(k)=0.0d0
               do l=1,3
               v(k)=v(k)+g(k,l,igp)*shell(l,js)
               end do
              end do
              err=abs(v(1)-shell(1,ns))+
     &            abs(v(2)-shell(2,ns))+
     &            abs(v(3)-shell(3,ns)) 
              if(err.lt.0.001)idisp=1
            end do
           end do
         if(idisp.eq.1)ns=ns-1
         end do
         end do         
         end do         
         print*,9**3,ns
         ms=ns
         do is=1,ns
            jgp=0
            do igp=1,ngp
            jgp=jgp+1
              do k=1,3
              angle(k,jgp)=0.0d0
               do l=1,3
               angle(k,jgp)=angle(k,jgp)+g(k,l,igp)*shell(l,is)
               end do
              end do
              idisp=0
              do kgp=1,jgp-1
              err=abs(angle(1,jgp)-angle(1,kgp))+
     &            abs(angle(2,jgp)-angle(2,kgp))+
     &            abs(angle(3,jgp)-angle(3,kgp)) 
              if(err.lt.0.001)idisp=1
              end do
              if(idisp.eq.1)jgp=jgp-1
              shell(4,is)=float(jgp)
            end do
         end do
         do is=1,ms 
c find total number of points for this shell:
                
         do js=is+1,ms
             if(shell(4,js).lt.shell(4,is))then
                   do k=1,4
                   sh=shell(k,js)
                   shell(k,js)=shell(k,is)
                   shell(k,is)=sh
                   end do
              end if
         end do
         print *,(shell(j,is),j=1,4)
         end do
         ns=ms
c        do is=1,ns
c        if(shell(4,is).ge.0.0d0)then
c          do js=is+1,ns
c          err=abs(shell(1,js)-shell(1,is))+
c    &         abs(shell(2,js)-shell(2,is))+
c    &         abs(shell(3,js)-shell(3,is)) 
c          if(err.lt.0.0001)shell(4,js)=-1.0d0
c          end do
c        end if
c        end do
c        ms=0
c        do is=1,ns
c        if(shell(4,is).ge.0.0)then
c        ms=ms+1
c         do j=1,4
c         shell(j,ms)=shell(j,is)
c         end do
c        end if
c        end do
         mtyp=0
         ns=ms
         detm=0.0d0
         do is=1,ns
                shell(5,is)=0.0D0
                do i=1  ,nint(shell(4,is))
                do j=i+1,nint(shell(4,is))
                   dist=(shell(1,i)-shell(1,j))**2+
     &                 +(shell(2,i)-shell(2,j))**2
     &                 +(shell(3,i)-shell(3,j))**2
                   dist=sqrt(dist)
                shell(5,is)=shell(5,is)+1.0d0/dist
                end do
                end do
           call  all_point(ngp,mgp,shell(1,is),g,angle)
c generate all equivalent sites:
c                mgp=0
c                do igp=1,ngp
c                mgp=mgp+1
c                 do i=1,3
c                   angle(i,mgp)=0.0d0
c                   do k=1,3
c                   angle(i,mgp)=angle(i,mgp)+g(i,k,igp)*shell(k,is)
c                   end do
c                 end do
c                 errmin=1.0D30
c                 do jgp=1,mgp-1
c                 err=abs(angle(1,mgp)-angle(1,jgp))
c    &               +abs(angle(2,mgp)-angle(2,jgp))
c    &               +abs(angle(3,mgp)-angle(3,jgp))
c                 errmin=min(err,errmin)
c                 end do
c                 if(errmin.lt.0.0001)mgp=mgp-1
c                end do
          CALL HARMONICS(121,mgp,LM,ANGLE,YLM,NPOL)
          det=1.0d0
          call find_rank(det,mxd,mgp,ylm,nrank)
          best=max(abs(det),best)
          if(abs(det).gt.0.01)then
c         write(400,*) -1000, ' Successful Shell',mgp,best
          shell(6,is)=abs(det)
c         write(400,35)is,mgp,abs(det),(shell(j,is),j=1,6),best
 35       format(2I4,f12.4,8f15.6)
c         do jgp=1,mgp
c         write(400,40)(ylm(jgp,ipol),ipol=1,16)
c40       format(f6.2,4x,3f6.2,4x,5f6.2,4x,7f6.2)
c         end do
          end if
         end do
         write(400,*)'ordered'
         ms=0
            do jgp=1,ngp
               do is=1,ns
                  if(nint(shell(4,is)).eq.jgp)iend=is
               end do
               do is=ns,1,-1
                  if(nint(shell(4,is)).eq.jgp)ibeg=is
               end do
               do is=ibeg,iend
                do js=is+1,iend
                   if(shell(6,js).gt.shell(6,is))then
                      do k=1,6
                      sss=shell(k,js)
                      shell(k,js)=shell(k,is)
                      shell(k,is)=sss
                      end do
                   end if
                end do
               end do
            end do
          ms=0
          do is=1,ns
          if(shell(6,is).gt.0.01)then
          ms=ms+1
              do j=1,6
              shell(j,ms)=shell(j,is)
              end do
          write(400,35)ms,mgp,(shell(j,ms),j=1,6)
          end if
          end do
            ns=ms
c look for mixed meshes:
          det_max=0.0d0
          ms=ns
          do is=1,ns
          do js=is+1,ns
           call  all_point(ngp,mgp1,shell(1,is),g,angle(1,1))
           call  all_point(ngp,mgp2,shell(1,js),g,angle(1,1+mgp1))
           mgp=mgp1+mgp2
          CALL HARMONICS(121,mgp,LM,ANGLE,YLM,NPOL)
          det=1.0d0
          call find_rank(det,mxd,mgp,ylm,nrank)
          if(abs(det).gt.0.01)then
          ms=ms+1
          det_max=abs(det)
          write(400,*)is,js,mgp,det, 'mixed shell'
          print     *,is,js,mgp,det, 'mixed shell'
                   do ks=1,mgp
                   write(400,35)ms,mgp,(angle(j,ks),j=1,3)
                   end do
          end if
          end do
          end do
         do is=1,ns
         errmin=1.0D30
          do ityp=1,mtyp
          errmin=min(abs(shell(4,is)-nmem(ityp,1)),errmin)
          end do
          if(errmin.gt.0.001)then
          mtyp=mtyp+1
          nmem(mtyp,1)=nint(shell(4,is))
          write(400,*)mtyp,nmem(mtyp,1)
          end if
         end do
         nelec=43          
         do ityp=1,mtyp
         nmem(ityp,2)=nelec/nmem(ityp,1)
         write(400,*)ityp,nmem(ityp,1),nmem(ityp,2)
         end do
         call srand(12345)
         isoln=1
         do jsoln=1,30
         n_elec=nelec
         do ityp=mtyp,2,-1
            ntme=n_elec/nmem(ityp,1)
            ncoin=int(rand()*ntme+0.5)
            nsoln(ityp,isoln)=ncoin  
            m_elec=n_elec
                  m_elec=m_elec-nsoln(ityp,isoln)*nmem(ityp,1)
c           write(400,55)ityp,nmem(ityp,1),ncoin,m_elec,
c    &(nsoln(jtyp,isoln),jtyp=1,mtyp)
            n_elec=m_elec
         end do
         nsoln(1,isoln)=m_elec
                 m_elec=0
                 do jtyp=1,mtyp
                 m_elec=m_elec+nsoln(jtyp,1)*nmem(jtyp,1)
                 end do
            nerr=1000
            do ksoln=1,isoln-1
            ierr=0
                do jtyp=1,mtyp
                ierr=ierr+abs(nsoln(jtyp,isoln)-nsoln(jtyp,ksoln))    
                end do
            nerr=min(ierr,nerr)
            end do
         if(nerr.ne.0)then
C00000000000000000000000000000000000000000000000000000000000000000000000
      write(line,60)m_elec,(nsoln(jtyp,isoln),nmem(jtyp,1),jtyp=1,mtyp) 
         do l=1,98
         if(line(l:l+1).eq.'00')line(l:l+4)='     '
         end do
         write(400, *)line
         isoln=isoln+1
         end if
         end do
         msoln=isoln
             do isoln=1      ,msoln
             do jsoln=isoln+1,msoln
                irank=0
                jrank=0
                do ityp=1,mtyp
                irank=irank+nsoln(ityp,isoln)
                jrank=jrank+nsoln(ityp,jsoln)
                end do
                if(jrank.lt.irank)then
                      do ityp=1,mtyp
                      isv=nsoln(ityp,isoln)
                      nsoln(ityp,isoln)=nsoln(ityp,jsoln)
                      nsoln(ityp,jsoln)=isv
                      end do
                end if 
             end do
             end do
         write(400,*)' '
         do isoln=1,msoln
      write(line,60)m_elec,(nsoln(jtyp,isoln),nmem(jtyp,1),jtyp=1,mtyp) 
         do l=1,98
         if(line(l:l+1).eq.'00')line(l:l+4)='     '
         end do
         write(400,*)line
         end do
         call flush(400)
         close(400)
 50      format(i5,8f12.4)
 55      format(4I4,'   ',25I4)
 60      format(I5,'    ',25(I2.2,'*',I2.2,' '))
         end
          subroutine find_rank(det,mxd,mgp,ylm,nrank)
          implicit real*8 (a-h,o-z)
          dimension ylm(121,mxd),squ(121,242)
          nret=0
           det=0.0d0
          if(mgp.eq. 4)nret=1
          if(mgp.eq. 6)nret=1
          if(mgp.eq. 8)nret=1
          if(mgp.eq.12)nret=1
          if(mgp.eq.16)nret=1
          if(nret.eq.0)return
          print*,'mgp:',mgp
c         read*,dummy
                if(mgp.eq. 4)n4=0   !mix L=0 and L=1
                if(mgp.eq. 6)n4=1   !mix L=1 and L=2
                if(mgp.eq. 8)n4=1   !mix L=1 and L=2
                if(mgp.eq.12)n4=4   !mix L=2 and L=3
                if(mgp.eq.16)n4=0   !mix L=0, 1, 2, 3
C         if(mgp.eq.12)then
          squ=0.0d0
          do igp=1,mgp
            do jgp=1,mgp
            squ(igp,jgp)=ylm(igp,jgp+n4)
            end do
          squ(igp,igp+mgp)=1.0d0 
          print 75,(squ(igp,jgp),jgp=1,mgp)
          end do
          print*,n4,mgp
c         read*,dummy
          det=1.0d0
          do igp=1,mgp
             sqm=0.0d0
             do jgp=igp,mgp
               if(abs(squ(jgp,igp)).gt.abs(sqm))then
               kgp=jgp
               sqm=squ(jgp,igp)
               end if
             end do
             det=det*sqm
             if(abs(det).le.1.0d-6)then
             print*,'sqm:',sqm,' det:',det
                   return
             end if
             do jgp=igp,mgp!*2
             sqs=squ(kgp,jgp)
             squ(kgp,jgp)=squ(igp,jgp)
             squ(igp,jgp)=sqs
             end do
             do kgp=igp,mgp!*2
             squ(igp,kgp)=squ(igp,kgp)/sqm
             end do
             do jgp=1,mgp
             sqs=squ(jgp,igp)
             if(jgp.ne.igp)then
             do kgp=igp,mgp!*2
             squ(jgp,kgp)=squ(jgp,kgp)-squ(igp,kgp)*sqs
             end do
             end if
             end do
          print 76,(squ(igp,kgp),kgp=1,mgp)
          end do 
          print*,'det:',det
c         read*,dummy
               if(abs(det).gt.detm)detm=det
 76       format(48f12.6)
c         end if 
 75      format(32f6.2)
         end
          subroutine atomic_radii
          dimension alp(100),con(100,10,10)
          dimension ncon(10,2),gint(0:30)
         dimension radii(3,100)
          character*7 name
          character*10 number
          character*80 line
          data number/'0123456789'/
          open(90,file='ISYMGEN')
          read(90,*)ntyp
          do ityp=1,ntyp
          ncon=0
          read(90,*)nz,mz
          read(90,*)
          read(90,*)natm
                   do iatm=1,natm
                   read(90,'(A7)')name
                   end do
          read(90,*)
          read(90,*)nalp
          lmx=0
          do i=1,2
          read(90,'(A)')line
                do j=1,80         
                do k=1,10
                if(line(j:j).eq.number(k:k))lst=j
                end do
                end do
          read(line,*,err=10,end=10)(ncon(j,i),j=1,10)
 10       continue
          do j=1,10
          if(ncon(j,i).gt.0)lmx=max(j,lmx)
          end do
c         print*,nalp,natm,line(1:lst),'!',(ncon(j,i),j=1,lmx)
          end do
          read(90,*)(alp(j),j=1,nalp)
           do l=1,lmx
           read(90,*)((con(j,i,l),j=1,nalp),i=1,ncon(l,1)+ncon(l,2))
           mcon=0
           do i=1,ncon(l,1)+ncon(l,2)
              nzero=0
              do j=1,nalp
              if(abs(con(j,i,l)).gt.0.000001)nzero=nzero+1
              end do
            if(nzero.gt.3)then
            mcon=mcon+1
            do j=1,nalp
            con(j,mcon,l)=con(j,i,l)
            end do
            end if
           end do
           ncon(l,1)=mcon
           ncon(l,2)=0
           end do
                   lmxx=lmx
                   do l=1,lmx
                   if(ncon(l,1).gt.0)lmxx=lmx
                   end do
                   lmx=lmxx
          print*,nz,name,(ncon(l,1),l=1,lmx)
          do l=1,lmx
          do i=1,ncon(l,1)
           print'(20G12.3)',(con(j,i,l),j=1,nalp)
          end do
          end do
c calculate radii:
               write(400,*)(ncon(l,1),l=1,lmx),nz,name
               nrad=0
               do l=1,lmx
               do i=1,ncon(l,1)
               chg=0.0d0
               rsq=0.0d0
                 do m=1,nalp
                 do n=1,nalp
                 gint(0)=sqrt(atan(1.0d0)/(alp(m)+alp(n)))
                       do k=1,lmx+4
                       gint(k)=(2*k-1)*gint(k-1)/(2.0*(alp(m)+alp(n)))
                       end do
                 chg=chg+con(m,i,l)*con(n,i,l)*gint(l)
                 rsq=rsq+con(m,i,l)*con(n,i,l)*gint(l+1)
                 end do
                 end do
                 print*,l,i,chg*16.0*atan(1.0d0),sqrt(rsq/chg)
               nrad=nrad+1
               radii(1,nrad)=i
               radii(2,nrad)=l-1
               radii(3,nrad)=sqrt(rsq/chg)
               end do
               end do
               radl=0.0d0
               do irad=1,nrad
               do jrad=irad+1,nrad
                   if(radii(3,jrad).lt.radii(3,irad))then
                      do m=1,3
                      rs=radii(m,irad)
                      radii(m,irad)=radii(m,jrad)
                      radii(m,jrad)=rs
                      end do
                   end if
               end do
               if(radii(3,irad).gt.1.5*radl)then
                    avg=0.0d0
                    knt=0
                    write(400,*)' '
               end if
               radl=radii(3,irad)
                  avg=avg+radl
                  knt=knt+1
               write(400,40)nint(radii(1,irad)+radii(2,irad)),
     &         nint(radii(2,irad)),2*nint(radii(2,irad))+1,
     &         avg/knt,radii(3,irad),name
               end do
          end do
 40            format(3I5,2F12.4,' ',a7)
          close(90)
          end
          subroutine  all_point(ngp,mgp,shell,g,angle)
          implicit real*8 (a-h,o-z)
          dimension g(3,3,*),angle(3,*),shell(6)
c generate all equivalent sites:
                 mgp=0
                 do igp=1,ngp
                 mgp=mgp+1
                  do i=1,3
                    angle(i,mgp)=0.0d0
                    do k=1,3
                    angle(i,mgp)=angle(i,mgp)+g(i,k,igp)*shell(k)
                    end do
                  end do
                  errmin=1.0D30
                  do jgp=1,mgp-1
                  err=abs(angle(1,mgp)-angle(1,jgp))
     &               +abs(angle(2,mgp)-angle(2,jgp))
     &               +abs(angle(3,mgp)-angle(3,jgp))
                  errmin=min(err,errmin)
                  end do
                  if(errmin.lt.0.0001)mgp=mgp-1
                 end do
            return
            end
