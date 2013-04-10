! 7.4.2013 Juraj Hasik

Program MPIRanGenTest

use rngmar
use LJSystem

include 'mpif.h'

!===================================================================================
character charMyId*2
character charNumberOfAtoms*3

integer i

REAL*8 DRMAX

Type(LjEnsamble) LjEns
Real*8 initialBeta
Real*8 maximalDisplacement

! ** Diagnostic variables **



!       John A. White "Lennard-Jones as a model for argon and test of extended renormalization group calculations", 
!	Journal of Chemical Physics 111 pp. 9352-9356 (1999)
!	Parameter Set #4
!	epsilon/K_b = 125.7 K
Real*8, parameter :: sigma = 3.4345 !Ang

! Minimum temperature reached during simulation in LJ units
REAL*8, Parameter :: minTemp = 0.039777
REAL*8, Parameter :: maxTemp = 0.318218
Integer, parameter :: numberOfAtoms = 55
Integer, parameter :: initialEquilibriation = 5000 

REAL*8, PARAMETER :: PI = 3.1415927
!===================================================================================

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myId, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcesses, ierr)

! ** Log info about current run **
if(myId.eq.0) then
  write(charNumberOfAtoms,'(I3.3)') numberOfAtoms
  open(7,file="ptMpi"//charNumberOfAtoms//".log", form="formatted", status="unknown")
  write(7,'("Parallel Tampering with Mpi")')
  write(7,'("Number of processes: ",I2)') numProcesses
  write(7,'("Number of Atoms: ",I3)') numberOfAtoms
  write(7,'("Temperature range: ",F10.6," - ",F10.6," [in LJ units]")') minTemp, maxTemp
end if

! ** Initialize PRNG **
call rmaset(-6,10,1,myId,'nexiste.pa')

write(charMyId,'(I2.2)') myId
open(6,file="ranGenTest"//charMyId//".dat", form="formatted", status="unknown")
write(6,'("LjEnsamble - governed by process with Id = ",I2)') myId
close(6)

! ** Initialze LJ Ensamble **
initialBeta = 1/(minTemp+dble(myId)*(maxTemp-minTemp)/dble(numProcesses))
call InitLjSystemBeta( LjEns, initialBeta, (myId-1), myId, (myId+1))
call InitLjSystemCoords( numberOfAtoms, LjEns )

! ** Log initial values to separate coressponding file **
open(6,file="ranGenTest"//charMyId//".dat", form="formatted", status="old", access="append")
write(6,'("Beta temperature = ",F10.5)') LjEns%Beta
write(6,'("Initial indexes by Beta Temperature: ",I2,"  ",I2,"  ",I2," [colder - my index - hotter]")') LjEns%higherBetaId, LjEns%myBetaId, LjEns%lowerBetaId
write(6,'("Initial Coordinates for N = ",I3," atoms")') numberOfAtoms
  do i=1, numberOfAtoms
    write(6,'(3F10.5)') LjEns%X(i), LjEns%Y(i), LjEns%Z(i)
  end do
close(6)

! ** equilibriate after initialization **
maximalDisplacement = sigma
do i=1, initialEquilibriation
  call sweepOverReplica( numberOfAtoms, sigma, maximalDisplacement, LjEns )
end do

call energyOfSystem( sigma, numberOfAtoms, LjEns)
! ** log data after initial equilibriation **
open(6,file="ranGenTest"//charMyId//".dat", form="formatted", status="old", access="append")
write(6,'("Energy after initial equilibriation = ",F10.5)') LjEns%V
write(6,'("Maximal displacement after initial equilibriation = ",F10.5)') maximalDisplacement
close(6)

! ** 
call tryReplicaSwap( numProcesses, LjEns, myId )
open(6,file="ranGenTest"//charMyId//".dat", form="formatted", status="old", access="append")
write(6,'("Extended Neighbour list: ",I2," ",I2," ",I2)') LjEns%neighbourList(1), LjEns%neighbourList(2), LjEns%neighbourList(3) 
close(6)

call MPI_FINALIZE(ierr)

end