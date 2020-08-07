C UTEP Electronic Structure Lab (2020)
        subroutine r_expect_val(max_pts,mxspn,nmsh,ipts,mx_grp,ngrp,
     &   Rmat,Rhog,Wmsh,Rmsh)
!----6------------------------------------------------------------------
!   Inputs: ipts, ngrt, rhog, wmsh
! Routine compute various expectations of electron density: <r^n> (n=-1,0,1,2,3)
! <r^2> moment is checked  against Gaussian (what one gets from Gaussian code).
! Quadrupole moments are also computed but they do not agree with Gaussian results.
! Atomic units are used everywhere.
! Shanon information entropy Int_ rho*log(rho) d^3r
!               Raja, Sept. 25, 2005

        implicit none
         integer max_pts, mxspn, nmsh,ipts,i,j,igrp, mx_grp,ngrp,krhog
        real*8   totqnum,rsq,rcube,rav,rinv,psave,aux,beriz,rdist
        PARAMETER(KRHOG=10)
        !real*8  RHOG(MAX_PTS,KRHOG,MXSPN),WMSH(NMSH)
        real*8  RHOG,wmsh,dipole,rmat,rmsh
        dimension RHOG(MAX_PTS,KRHOG,MXSPN),WMSH(NMSH),DIPOLE(3)
        dimension RMAT(3,3,MX_GRP),Rmsh(3,max_pts)
        real*8  rsq_aux,rrot,q_xx,qxy,q_yy,q_yz,q_zz,qxz,a_xx,a_yy,a_zz
       


         TOTQNUM=0.0D0
       DO I=1,3
        DIPOLE(I)=0.0D0
       END DO

        rsq = 0.d0
        rcube = 0.d0
        rav = 0.d0
        rinv = 0.d0
        q_xx = 0.d0
        q_yy = 0.d0
        q_zz = 0.d0
       DO 120 IPTS=1,NMSH
         PSAVE=RHOG(IPTS,1,1)*WMSH(IPTS) 
         TOTQNUM=TOTQNUM+PSAVE
        DO I=1,3
          RROT=0.0D0
          rsq_aux = 0.0d0
          rdist=0.d0
          a_xx = 0.d0
          a_yy = 0.d0
          a_zz = 0.d0
         DO IGRP=1,NGRP
          DO J=1,3
           RROT=RROT+RMAT(J,I,IGRP)*RMSH(J,IPTS)
           rdist =  rdist + (Rmat(j,i,igrp)*Rmsh(j,ipts))**2  ! Compute r**2
          END DO
           a_xx =  a_xx + (Rmat(1,i,igrp)*Rmsh(1,ipts))**2  ! r_x**2 (Xcomponent only)
           a_yy =  a_yy + (Rmat(2,i,igrp)*Rmsh(2,ipts))**2  ! r_y**2 (Xcomponent only)
           a_zz =  a_zz + (Rmat(3,i,igrp)*Rmsh(3,ipts))**2  ! r_z**2 (Xcomponent only)
         END DO
         DIPOLE(I)=DIPOLE(I)+RROT*PSAVE
           q_xx = q_xx + a_xx*PSAVE
           q_yy = q_yy + a_yy*PSAVE
           q_zz = q_zz + a_zz*PSAVE
           rdist = a_xx + a_yy + a_zz
           rdist = dsqrt(rdist)
           rsq = rsq + (rdist**2)*PSAVE                    ! <r**2>
           rcube = rcube + ((rdist)**3)*PSAVE              ! <r**3>
           rav = rav + (rdist)*PSAVE                       ! <r>
           if (rdist .gt. 1.d-10) then
            rinv = rinv + PSAVE/rdist                      ! <1/r>
           endif
        END DO
  120  CONTINUE 
           q_xx = (3.d0*q_xx - rsq)*0.5d0
           q_yy = (3.d0*q_yy - rsq)*0.5d0
           q_zz = (3.d0*q_zz - rsq)*0.5d0

        WRITE(6,*)'SUB: Expectation value: <r>    =',rav
        WRITE(6,*)'TOTQNUM =',TOTQNUM
        WRITE(6,*)'Qxx value:                =',q_xx
        WRITE(6,*)'Qyy value:                =',q_yy
        WRITE(6,*)'Qzz value:                =',q_zz
        WRITE(6,*)'Expectation value: <r>    =',rav
        WRITE(6,*)'Expectation value: <r**2> =',rsq
        WRITE(6,*)'Expectation value: <r**3> =',rcube
        WRITE(6,*)'Expectation value: <1/r>  =',rinv
       DO I=1,3
        DIPOLE(I)=DIPOLE(I)/NGRP
       END DO

          return 
          end
