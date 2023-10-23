C UTEP Electronic Structure Lab (2020)
       subroutine allocate_sic(mode)
       use pot_dens,only : COULOMB,RHOG
       use mesh1,only : NMSH
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:46 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IERR
       logical, intent(in) :: mode

       if(mode) then
#ifdef MPI
        CALL SENDDATA(108)
#else
        ALLOCATE(COULOMB(NMSH*MXSPN),STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*)'ERROR ALLOCATING COULOMB'
        ALLOCATE(RHOG(NMSH,KRHOG,MXSPN),STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*)'ERROR ALLOCATING RHOG'
#endif
       else
#ifdef MPI
        CALL SENDDATA(109)
#else
        DEALLOCATE(RHOG)
        DEALLOCATE(COULOMB)
#endif
       end if
       return
       end
