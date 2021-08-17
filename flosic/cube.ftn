! Simple fortran program to evaluate the density, gradients, and tau 
! on an array of points. This is meant to be used as post-processing.
! the integer number of points nxyz and the xyz coordinates in xyzmsh 
!  are read from the input stream. 
!
!  juan e. peralta and duyen b. nguyen, CMU 2021
!   
!

      program runcube
      implicit none 
      real*8, dimension (:,:), allocatable :: xyzmsh  
      real*8, dimension (:,:,:), allocatable :: xyzg 

      integer :: nxyz,i

! read the number of point to be read
      read*, nxyz 

      allocate (xyzmsh(3,nxyz))


! read the points xyz
      Do i = 1, nxyz
         read*, xyzmsh(1,i), xyzmsh(2,i), xyzmsh(3,i)
         enddo

      allocate (xyzg(nxyz,11,2))



! do the work 
      call cube(nxyz, xyzmsh,xyzg) 


! print the density and its gradients of two spin channels
      do i =1, nxyz
         write(6,*) 'dens= ',xyzmsh(1,i), xyzmsh(2,i),xyzmsh(3,i), 
     &    xyzg(i,1,1), xyzg(i,2, 1), xyzg(i,3, 1), xyzg(i,4,1),
     &    xyzg(i,5,1), xyzg(i,6, 1), xyzg(i,7, 1), xyzg(i,8,1),
     &    xyzg(i,9,1), xyzg(i,10,1), xyzg(i,11,1),
     &    xyzg(i,1,2), xyzg(i,2, 2), xyzg(i,3, 2), xyzg(i,4,2),
     &    xyzg(i,5,2), xyzg(i,6, 2), xyzg(i,7, 2), xyzg(i,8,2),
     &    xyzg(i,9,2), xyzg(i,10,2), xyzg(i,11,2)
         enddo





      deallocate (xyzmsh)  
      deallocate (xyzg)  

      return
      end 


! do the work here
 
      subroutine cube(nxyz, xyzmsh, xyzg)
C     ------------------------------------------------------------------
      use common2,only : IGGA

       INCLUDE  'PARAMA2'  
       INTEGER :: nxyz
       real*8 xyzmsh,xyzg,TIME1,RDUMMY,RVEC(1)

      SAVE

      LOGICAL EXIST,FAILED
      dimension xyzmsh(3,nxyz), xyzg(nxyz,11,2)
C*****************************************************************

      IGGA(1) = 3
      IGGA(2) = 3
      CALL FGMAT   !READS THE GRPMAT FILE AND CHECKS FOR A LARGER GROUP
      CALL CREPMAT !DETERMINING SYMMETRY GROUP
      CALL READINP !READ ATOM DATA FROM INPUT & ELECTRIC FIELD FROM EFIELD
      CALL TESTBAS !TESTING BASIS SET
C
C
      CALL READWF(FAILED) !READS IN WAVEFUNCTIONS TO ALLOW RESTARTING A
                          !CALCULATION AT A SELF-CONSISTENT END POINT
                          !FOR FINE-TUNING A SET OF FIT GAUSSIANS AND
                          !FOR GETTING GOOD NEW STARTING POINTS

! this is the main worker 
      Call  DENS_CUBE(TIME1,xyzmsh,nxyz,xyzg) ! GETTING THE DENSITY


      RETURN
      END



C*******************************************************************
C*******************************************************************

