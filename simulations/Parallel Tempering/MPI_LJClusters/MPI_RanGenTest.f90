Program MPIRanGenTest

!use RanGen
include 'mpif.h'

character charMyId*2
integer i
real*8 ranNumber

integer, parameter :: seed = 0
real*8 DUMMY

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myId, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcesses, ierr)

write(charMyId,'(I2.2)') myId
!call ZBQLINI(seed)
call rmaset(-6,10,1,myId,'nexiste.pa')

open(6,file="ranGenTest"//charMyId//".dat", form="formatted", status="unknown")
write(6,'("myId = ",I2.2)') myId
close(6)

do i=1, 10
  call ranmar(ranNumber)
  open(6,file="ranGenTest"//charMyId//".dat", form="formatted", status="old", access="append")
  write(6,'(F10.5)') ranNumber
  close(6)
end do

call MPI_FINALIZE(ierr)

end

include 'ranmar.f90'
include 'rmaset.f90'