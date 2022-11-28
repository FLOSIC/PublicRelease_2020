! UTEP Electronic Structure Lab (2020)
!A subroutine to clean out ZPOTXXX files and MGGAXXX files
!Those files are created during a calculation and recommended to 
!be cleaned out upon starting a new calculation.
subroutine delete_oldfiles

integer :: INDEX
logical :: EXIST
character*12 :: ZPOTSTR, MGGASTR
!
! REMOVE OLD EVALUE FILES AND CREATE DUMMY FILE FOR FIRST ITERATION
!
INDEX=0
300   INDEX=INDEX+1
! WRITE(ZPOTSTR,'(A,I3.3)')'ZPOT',INDEX
!<LA: updated here to work correctly
 WRITE(ZPOTSTR,'(A,I4.4)')'ZPOT',INDEX
 INQUIRE(FILE=ZPOTSTR,EXIST=EXIST)
 IF (EXIST) THEN
  OPEN(198,FILE=ZPOTSTR,FORM='UNFORMATTED',STATUS='OLD')
  CLOSE(198,STATUS='DELETE')
  GOTO 300
 END IF
CONTINUE

INDEX=0
400   INDEX=INDEX+1
! WRITE(MGGASTR,'(A,I3.3)')'MGGA',INDEX
 WRITE(MGGASTR,'(A,I4.4)')'MGGA',INDEX
 INQUIRE(FILE=MGGASTR,EXIST=EXIST)
 IF (EXIST) THEN
  OPEN(298,FILE=MGGASTR,FORM='UNFORMATTED',STATUS='OLD')
  CLOSE(298,STATUS='DELETE')
  GOTO 400
 END IF
CONTINUE


end subroutine delete_oldfiles
