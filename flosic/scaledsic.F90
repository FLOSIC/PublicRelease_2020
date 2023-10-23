
!read scaling factor from file. 
!Assign them to TAUW_TAU or BETA in the module
!MPI bcast 

! input file tauw_tau1.dat or beta1.dat
!       module scaledpzsic
! output module TAUW_TAU or BETA array

subroutine readrsicfactor(IORBX,LSPX)

 use mesh1,only    : WMSH,RMSH,NMSH
 use pot_dens,only : COULOMB,RHOG
!MGGA module
 USE XTMP2A,ONLY : TAUW_TAU,BETA
 use scaledpzsic,only : LSICON,ORBSICON,GSICON,AVGSICON,SDSICON, &
                        AVGSIC,SICEXC,SDSICSAV,oneshot,scaledsic
 INCLUDE  'PARAMA2'
 integer,intent(in) :: IORBX,LSPX
 integer :: IMSH
 real(8) :: FPARAMA,FPARAMB,FPARAMC,CHI,POW,SDSCALE,SDSCALE1, &
            SDSCALE2


!Scaling
!Read RSIC factor here 
 IF(LSICON.or.ORBSICON.or.SDSICON) THEN
  IF(LSPX .EQ. 1) THEN
   OPEN(214,FILE='tauw_tau1.dat',status='OLD')
  ELSE
   OPEN(214,FILE='tauw_tau2.dat',status='OLD')
  END IF
  REWIND(214)
  allocate(TAUW_TAU(NMSH))
  READ(214,*) (TAUW_TAU(IMSH),IMSH=1,NMSH)
  !TAUW_TAU(:)=1.0d0 !Debugging - restore pure PZSIC
  CLOSE(214)
! Scuseria's Orbital-Scaling start
  if(ORBSICON) then
   !open(888,file='orbsicpower',form='formatted',status='old')
   !rewind(888)
   !read(888,*) POW
   !close(888)
   POW=1.0d0
   CHI=0.0d0
   DO IMSH=1,NMSH
    CHI = CHI + (TAUW_TAU(IMSH)**POW)*RHOG(IMSH,1,1)*WMSH(IMSH)
   END DO
   print *,"OrbSIC scaling factor (Spin,Orb,Chi)",LSPX,IORBX,CHI
   !Overwrite TAUW_TAU array with the same value (for simplicity)
   DO IMSH=1,NMSH
    TAUW_TAU(IMSH) = CHI
   END DO
  end if
! Orbital-Scaling end
  CALL SENDDATA(216) !Broadcast TAUW_TAU
  !DO IMSH=1,NMSH
  ! COULOMB(IMSH)=(TAUW_TAU(IMSH)**1)*COULOMB(IMSH)
  !END DO
 END IF
!John's GSIC modSCAN
 IF(GSICON) THEN
!In case of modSCANx, read the fitting params.
! open(888,file='modscanfitparam',form='formatted',status='old')
! rewind(888)
! read(888,*) FPARAMA,FPARAMB,FPARAMC
! close(888)

  FPARAMA = 0.349741D-01
  FPARAMB = -0.281332D00    
  FPARAMC = 0.127402D00 
  
  IF(LSPX .EQ. 1) THEN
   OPEN(214,FILE='beta1.dat',STATUS='OLD')
  ELSE
   OPEN(214,FILE='beta2.dat',STATUS='OLD')
  END IF
  REWIND(214)
  allocate(BETA(NMSH))
  READ(214,*) (BETA(IMSH),IMSH=1,NMSH)
  CLOSE(214)
  DO IMSH=1,NMSH
   BETA(IMSH) = 1.0d0+FPARAMA*BETA(IMSH) &
                     +FPARAMB*BETA(IMSH)*BETA(IMSH) &
                     +FPARAMC*BETA(IMSH)*BETA(IMSH)*BETA(IMSH)
  END DO
  CALL SENDDATA(214) !Broadcast beta
 END IF

end subroutine readrsicfactor

!###############################


! inout COULOMB,POT
! in    module TAUW_TAU,BETA,SICEXC
! out   module SDSICSAV 
subroutine scalepotential(IORBX,LSPX)

 use common2,only  : ERGXL,ERGXN,ERGCL,ERGCN
 use pot_dens,only : COULOMB
 use mixpot1,only  : POT=>POTOUT
 use mesh1,only    : WMSH,RMSH,NMSH
 USE XTMP2A,ONLY : TAUW_TAU,BETA
 use scaledpzsic,only : LSICON,ORBSICON,GSICON,AVGSICON,SDSICON, &
                        AVGSIC,SICEXC,SDSICSAV,oneshot,scaledsic
 implicit none
 integer,intent(in) :: IORBX,LSPX
 integer :: IMSH
 real(8) :: POW,SDSCALE,SDSCALE1,SDSCALE2        

 if(LSICON.OR.ORBSICON) then 
  !Scaling COULOMB and POT 
  DO IMSH=1,NMSH
   COULOMB(IMSH)=(TAUW_TAU(IMSH))*COULOMB(IMSH)
   !POT(IMSH)=(TAUW_TAU(IMSH))*POT(IMSH)
   !Scaling of Vxc is done inside subvlxc
  END DO
 end if

 if(SDSICON) then
!obtain exterior scaling
  SDSCALE1=0.0d0
  SDSCALE2=0.0d0
  POW=3.0d0
  DO IMSH=1,NMSH
   SDSCALE1=SDSCALE1+SICEXC(IMSH)*WMSH(IMSH) &
             *(POW*TAUW_TAU(IMSH)**POW       &
      -(POW-1.0d0)*TAUW_TAU(IMSH)**(POW+1.0d0))
   SDSCALE2=SDSCALE2+SICEXC(IMSH)*WMSH(IMSH)
  END DO
  SDSCALE=SDSCALE1/SDSCALE2
  print *,"SDSIC",IORBX,POW,SDSCALE
!SD scaling -multiple energy
  ERGXL = ERGXL*SDSCALE
  ERGXN = ERGXN*SDSCALE
  ERGCL = ERGCL*SDSCALE
  ERGCN = ERGCN*SDSCALE
!SD scaling - multiply potential here
!Coulomb energy is scaled with this
  DO IMSH=1,NMSH
   POT(IMSH)=POT(IMSH)*SDSCALE
   COULOMB(IMSH)=COULOMB(IMSH)*SDSCALE
  END DO
  SDSICSAV(IORBX,LSPX)=SDSCALE
 end if

end subroutine scalepotential

!#################################

subroutine deallocatersicfactor
 USE XTMP2A,ONLY : TAUW_TAU,BETA
 use scaledpzsic,only : LSICON,ORBSICON,GSICON,AVGSICON,SDSICON
 implicit none

  IF(LSICON.or.ORBSICON.or.SDSICON) THEN
#ifndef MPI
   deallocate(TAUW_TAU) !Serial
#else
   CALL SENDDATA(217)    !MPI
#endif
  END IF
  IF(GSICON) THEN
#ifndef MPI
   deallocate(BETA)
#else
   CALL SENDDATA(215)   !Deallocate BETA
#endif
  END IF
end subroutine deallocatersicfactor
