program testMT

use mt_stream
use mpi

Type(mt_state) mts, myMts
integer i, seed, ierr
integer myId, numProcesses

call set_mt19937
call new(mts)
seed = (113 + 30491*myId)
call init(mts,seed)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myId, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcesses, ierr)

call create_stream(mts, myMts, myId)

do i=1, 10
  write(*,'("ID ",I2," : "F10.7)') myId, genrand_double1(myMts)
end do

call MPI_FINALIZE(ierr)

end program