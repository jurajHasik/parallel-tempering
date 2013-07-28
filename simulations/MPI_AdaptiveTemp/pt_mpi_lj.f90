! 7.4.2013 Juraj Hasik

Program LjClusteringWithMPI

use mt_stream
use LJSystem
use mpi

!===================================================================================
! ** Output files spcifieres **
character charMyId*2
character charMyBetaId*2
character charNumberOfAtoms*3

integer i, j, k, l, seed
integer itag, itagReal, tick

! ** keeps state with minimal energy during whole simulation **
Real*8 minEnergy
Real*8, allocatable, dimension(:) :: minEx, minEy, minEz
Real*8, allocatable, dimension(:) :: minEnergyOfReplicas

! **  
Type(LjEnsamble) LjEns
Real*8, allocatable, dimension(:) :: betaArray, energyOfReplicas

integer, allocatable, dimension(:) :: currentBetaIndexes, acceptanceRates, IPI_B
Real*8, allocatable, dimension(:) :: BASUM

! ** Rng **
Type(mt_state) :: mts, myMts


!       John A. White "Lennard-Jones as a model for argon and test of extended renormalization group calculations", 
!	Journal of Chemical Physics 111 pp. 9352-9356 (1999)
!	Parameter Set #4
!	epsilon/K_b = 125.7 K
!	sigma = 3.4345 Ang
Real*8, parameter :: sigma = 1.0 !3.4345D-10 m 
Real*8, parameter :: LJepsilon = 1.0 !1.735478D-21 J

! Minimum temperature reached during simulation in LJ units
REAL*8, Parameter :: minTemp = 0.007955 !0.039777 5 K
REAL*8, Parameter :: maxTemp = 0.397772 !0.397772 50K 0.357995 45K 0.318218 40K 0.278441 35 K 0.238663 - 30 K
Integer, parameter :: numberOfAtoms = 75

Integer, parameter :: initialEquilibriation = 100000
Integer, parameter :: numberOfSwapTriesPriorToRecursion = 300
Integer, parameter :: numberOfBetaRecursions = 20
Integer, parameter :: betaDistributionEq = 100
Integer, parameter :: productionSweeps = 10000000
Integer, parameter :: productionSwapFrequency = 100
Integer, parameter :: loggingFrequency = 10000

!===================================================================================

! ** Initialize base Marsenne Twister stream **
call set_mt19937
call new(mts)
seed = 1
call init(mts,seed)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myId, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcesses, ierr)

! ** assign part of rng stream to each process **
call create_stream(mts, myMts, myId)

! ** Log info about current run **
if(myId.eq.0) then
  write(charNumberOfAtoms,'(I3.3)') numberOfAtoms
  open(7,file="ptMpi"//charNumberOfAtoms//".log", form="formatted", status="unknown")
  write(7,'("Lennard Jones parameters: ")')
  write(7,'("epsilon = ",F10.5," kg.Ang2.s-2")') LJepsilon
  write(7,'("sigma = ",F10.5," Ang")') sigma
  write(7,'("Parallel Tampering with Mpi")')
  write(7,'("Number of processes: ",I2)') numProcesses
  write(7,'("Number of Atoms: ",I3)') numberOfAtoms
  write(7,'("Temperature range: ",F10.6," - ",F10.6," [in LJ units]")') minTemp, maxTemp
  write(7,'("Initial Equilibriation - number of sweeps: ",I10)') initialEquilibriation
  write(7,'("Number Of Swap Tries Prior To Recursion: ",I10)') numberOfSwapTriesPriorToRecursion
  write(7,'("Number Of Beta Recursions: ",I10)') numberOfBetaRecursions
  write(7,'("Equilibriation after Beta recursion - number of sweeps: ",I10)') betaDistributionEq
  write(7,'("Number of Production sweeps: ",I12)') productionSweeps
  write(7,'("Production Swap attempt Frequency: [every * sweeps] ",I10)') productionSwapFrequency
  write(7,'("seed = "I10)') seed
  close(7)
end if

write(charMyId,'(I2.2)') myId
open(6,file="ptArData"//charMyId//".dat", form="formatted", status="unknown")
write(6,'("LjEnsamble - governed by process with Id = ",I2)') myId
close(6)

! ** Initialze LJ Ensamble **
acceptance = 0
if(myId.eq.0) write(* ,'("Initialization ensambles [Beta temperature and coordinates]")')
allocate(betaArray(0:(numProcesses-1)))
do i=0, numProcesses-1
  !betaArray(i) = 1/(minTemp+dble(numProcesses-1-i)*(maxTemp-minTemp)/dble(numProcesses-1))
  betaArray(i)=1/(maxTemp*((minTemp/maxTemp)**(dble(i)/(numProcesses-1))))
end do
call InitLjSystemBeta( LjEns, betaArray, numProcesses, myId, sigma, sigma, LJepsilon)
call InitLjSystemCoords( numberOfAtoms, LjEns, myMts )

write(charMyBetaId,'(I2.2)') LjEns%myBetaId
open(8,file="ptAr"//charMyBetaId//".xyz", form="formatted", status="unknown")
close(8)

if(myId.eq.0) then
  open(7,file="ptMpi"//charNumberOfAtoms//".log", form="formatted", status="old", access="append")
  write(7,'("Initialized Beta Temperatures for Processes as: [hotter -> colder]")')
  do i=0, numProcesses-2
    write(7,'(F10.4)', ADVANCE='NO') betaArray(i)
  end do
  write(7,'(F10.4)') betaArray(numProcesses-1)
  close(7)
end if

! ** Log initial values to separate coresponding file **
open(6,file="ptArData"//charMyId//".dat", form="formatted", status="old", access="append")
write(6,'("Initial Beta temperature = ",F10.5)') LjEns%staticBetaList(LjEns%myBetaId)
write(6,'("Initial indexes for Neighbours [in Beta Temperature]: ",I2,"  ",I2," [hotter - colder]")') LjEns%betaNeighbours(1), LjEns%betaNeighbours(2)
write(6,'("Initial Coordinates for N = ",I3," atoms")') numberOfAtoms
  do i=1, numberOfAtoms
    write(6,'(3F10.5)') LjEns%X(i), LjEns%Y(i), LjEns%Z(i)
  end do
close(6)

! ** equilibriate after initialization and set optimal maximalDisplacement for trial move**
if(myId.eq.0) write(* ,'("Starting initialization equilibriation and setting optimal maximalDisplacement for trial move")')
do i=1, initialEquilibriation
  call sweepOverReplica( numberOfAtoms, LjEns, myMts )
end do
if(myId.eq.0) write(* ,'("Initialization DONE")')

call energyOfSystem( numberOfAtoms, LjEns)
! ** log data after initial equilibriation **
open(6,file="ptArData"//charMyId//".dat", form="formatted", status="old", access="append")
write(6,'("Energy after initial equilibriation = ",F12.5)') LjEns%V
write(6,'("Maximal displacement after initial equilibriation = ",F10.5)') LjEns%maxDisplacement
close(6)

open(6,file="ptArData"//charMyId//".dat", form="formatted", status="old", access="append")
write(6,'("Starting indexes for Neighbours [in Beta Temperature]: ",I2," ",I2," ",I2," [hotter - myBetaId - colder]")') LjEns%betaNeighbours(1), LjEns%myBetaId, LjEns%betaNeighbours(2)
close(6)

itag = 1
itagReal = 10001
allocate( currentBetaIndexes(0:(numProcesses-1)), acceptanceRates(0:(numProcesses-1)), IPI_B(0:(numProcesses-1)) )
allocate( BASUM(0:(numProcesses-1)) )



! ** beta temperature rescaling **
if(myId.eq.0) write(* ,'("Starting recursion to rescale Beta temperature intervals")')
WEIGHt = 0.0
ICALL = 0
IADD = 0
do i=1, numberOfBetaRecursions 
  LjEns%acceptance = 0
  do j=1, numberOfSwapTriesPriorToRecursion
    do k=1, betaDistributionEq
      call sweepOverReplica( numberOfAtoms, LjEns, myMts )
    end do
    call tryReplicaSwap( LjEns, numProcesses, mod(j,3), itag, itagReal, myMts )
    open(6,file="ptArData"//charMyId//".dat", form="formatted", status="old", access="append")
    write(6,'("New indexes for Neighbours [in Beta Temperature]: ",I2," ",I2," ",I2," [hotter - myBetaId - colder]")') LjEns%betaNeighbours(1), LjEns%myBetaId, LjEns%betaNeighbours(2)
    write(6,'("New Beta = ",F10.5)') LjEns%staticBetaList(LjEns%myBetaId)    
    close(6)
  end do
  ! ** beta recursion **
  call MPI_ALLGATHER(LjEns%myBetaId, 1, MPI_INTEGER, currentBetaIndexes(0), 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  call MPI_ALLGATHER(LjEns%acceptance, 1, MPI_INTEGER, acceptanceRates(0), 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  
  if(myId.eq.0) then
    open(7,file="ptMpi"//charNumberOfAtoms//".log", form="formatted", status="old", access="append")
    write(7,'("Acceptance Rates after ",(I3)," Swap attempts for Processes: ")') numberOfSwapTriesPriorToRecursion
    do l=0, numProcesses-2
      write(7,'(I4)', ADVANCE='NO') acceptanceRates(l)
    end do
    write(7,'(I4)') acceptanceRates(numProcesses-1)
    close(7)
  end if
  
  call betaRecursion1(LjEns%staticBetaList, BASUM, acceptanceRates, currentBetaIndexes, IPI_B, numProcesses, numberOfBetaRecursions,weight,icall,iadd)
  
  if(myId.eq.0) then
    open(7,file="ptMpi"//charNumberOfAtoms//".log", form="formatted", status="old", access="append")
    write(7,'("Beta Temperatures after recursion: [hotter -> colder]")')
    do l=0, numProcesses-2
      write(7,'(F10.4)', ADVANCE='NO') LjEns%staticBetaList(l)
    end do
    write(7,'(F10.4)') LjEns%staticBetaList(numProcesses-1)
    close(7)
  end if
    
end do
if(myId.eq.0) write(* ,'("Beta recursion DONE")')

call energyOfSystem( numberOfAtoms, LjEns)
! ** log data after beta recursion **
open(6,file="ptArData"//charMyId//".dat", form="formatted", status="old", access="append")
write(6,'("Energy after Beta Recursion = ",F12.5)') LjEns%V
write(6,'("Maximal displacement after Beta Recursion = ",F10.5)') LjEns%maxDisplacement
close(6)

! ** Production **
allocate(minEx(numberOfAtoms),minEy(numberOfAtoms),minEz(numberOfAtoms))
tick = productionSweeps / 100
itag = 1
itagReal = 10001
minEnergy = 0.0
j=0
if(myId.eq.0) write(* ,'("Entering production run ... ")')
do i=0, productionSweeps
 
    call sweepOverReplica( numberOfAtoms, LjEns, myMts )
  
    if (minEnergy .gt. LjEns%V) then
      minEnergy = min(minEnergy, LjEns%V)
      minEx = LjEns%x
      minEy = LjEns%y
      minEz = LjEns%z
    end if
  
    ! ** Data logging **
    if(mod(i,loggingFrequency) .eq. 0) then
      write(charMyBetaId,'(I2.2)') LjEns%myBetaId
      open(8,file="ptAr"//charMyBetaId//".xyz", form="formatted", status="old", access="append")
      write(8,'(I3)') numberOfAtoms
      write(8,'("Energy = ",F12.5," || minimal Energy =",F12.5," at Beta: ",F10.5)') LjEns%V, minEnergy, betaArray(LjEns%myBetaId)
      do k=1, numberOfAtoms
	write(8,'("Ar ",3F10.6)') LjEns%x(k), LjEns%y(k), LjEns%z(k)
      end do
      close(8)
    end if
        
    if(mod(i,productionSwapFrequency).eq.0) then
      call tryReplicaSwap(LjEns, numProcesses, mod(j,3), itag, itagReal, myMts )
      j=j+1
    end if
    
    if(mod(i,tick) .eq. 0) then
      if(myId.eq.0) write(*,'("-")', ADVANCE='NO')
    end if
    
end do
if(myId.eq.0) write(* ,'("Finished ! ")')

call energyOfSystem( numberOfAtoms, LjEns)
open(6,file="ptArData"//charMyId//".dat", form="formatted", status="old", access="append")
write(6,'("Beta = ",F10.5)') LjEns%staticBetaList(LjEns%myBetaId)
write(6,'("Final energy = ",F12.5)') LjEns%V
write(6,'("Lowest energy configuration: Energy = ",F10.5)') minEnergy
do i=1, numberOfAtoms
  write(6,'("Ar ",3F10.6)') minEx(i), minEy(i), minEz(i)
end do
close(6)

allocate(energyOfReplicas(0:(numProcesses-1)))
call MPI_GATHER(LjEns%V, 1, MPI_REAL8, energyOfReplicas(0), 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
allocate(minEnergyOfReplicas(0:(numProcesses-1)))
call MPI_GATHER(minEnergy, 1, MPI_REAL8, minEnergyOfReplicas(0), 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

if(myId.eq.0) then
  open(7,file="ptMpi"//charNumberOfAtoms//".log", form="formatted", status="old", access="append")
  write(7,'("Final Energy of Replicas")')
  do i=0, numProcesses-1
    write(7, '(2F12.5)') energyOfReplicas(i), minEnergyOfReplicas(i)
  end do
  close(7)
endif

call MPI_FINALIZE(ierr)

end

SUBROUTINE tryReplicaSwap ( LjEns, numProcesses, offset, itag, itagReal, mts )
        
!    *******************************************************************
!    ** Attepmts an exchange of Beta temperatures                     **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    **								      **
!    *******************************************************************

! Copyright, Bernd Berg, December 20, 2001.
! Parallel tempering for the Potts models heat bath code.
! MPI implementation of the beta exchange.
      
      use mt_stream
      use LJSystem
      use mpi
      
      Type(LjEnsamble) LjEns
      integer numProcesses, offset
      Type(mt_state) mts
      
      integer itag, itagReal
      integer status(MPI_STATUS_SIZE)
      
      logical mcCriterionForReplicaSwap
      integer my_ind1, my_ind2, my_ind3
      integer ndest3l, nsend3l, ndest3r, nsend3r, ndest1, NRECV3R, NRECV3L
            
!    *******************************************************************
      call MPI_COMM_RANK(MPI_COMM_WORLD, myId, ierr)
      
      itag=mod(itag+1,9000)
      itagReal=10000+mod(itagReal+1,10000)

      MY_IND1=MOD(LjEns%myBetaId,3)
      MY_IND2=MOD(LjEns%myBetaId+2,3) ! MY_B-1+3
      MY_IND3=MOD(LjEns%myBetaId+1,3) ! MY_B-2+3

      IF(MY_IND1.EQ.offset) THEN ! Processes MY_IND1.
        NDEST3L=LjEns%betaNeighbours(1)
        NSEND3L=-2
        
        IF(LjEns%betaNeighbours(2).LT.numProcesses) THEN
          LjEns%toSend(1)=LjEns%betaNeighbours(1)
          LjEns%toSend(2)=0 ! iact *former energy*
	  LjEns%toSend(3)=LjEns%acceptance
          LjEns%toSendReals(1) = LjEns%V
          
          CALL MPI_SEND(LjEns%toSend,3,MPI_INTEGER,LjEns%betaNeighbours(2),itag,MPI_COMM_WORLD,IERR)
          CALL MPI_SEND(LjEns%toSendReals,1,MPI_REAL8,LjEns%betaNeighbours(2),itagReal,MPI_COMM_WORLD,IERR)
          CALL MPI_RECV(LjEns%toRecieve,2,MPI_INTEGER,LjEns%betaNeighbours(2),itag,MPI_COMM_WORLD,STATUS,IERR)
          
          IF(LjEns%toRecieve(1).NE.-2) THEN
            NSEND3L=LjEns%betaNeighbours(2)
            LjEns%betaNeighbours(1)=LjEns%betaNeighbours(2)
            LjEns%betaNeighbours(2)=LjEns%toRecieve(1)
            LjEns%acceptance=LjEns%toRecieve(2) !**acceptance rate
            LjEns%myBetaId=LjEns%myBetaId+1
          END IF
        
        END IF
        IF(NDEST3L.GE.0) then
	  CALL MPI_SEND(NSEND3L,1,MPI_INTEGER,NDEST3L, itag,MPI_COMM_WORLD,IERR)
	end if
      END IF

      IF(MY_IND2.EQ.offset) THEN ! Processes MY_IND2.
        NDEST3R=LjEns%betaNeighbours(2)
        NSEND3R=-2
        
        IF(LjEns%betaNeighbours(1).GE.0) THEN
          CALL MPI_RECV(LjEns%toRecieve,3,MPI_INTEGER,LjEns%betaNeighbours(1),itag, MPI_COMM_WORLD,STATUS,IERR)
          CALL MPI_RECV(LjEns%toRecieveReals,1,MPI_REAL8,LjEns%betaNeighbours(1),itagReal, MPI_COMM_WORLD,STATUS,IERR)
          
          NDEST1=LjEns%betaNeighbours(1)
          
          IF(mcCriterionForReplicaSwap(LjEns%staticBetaList(LjEns%myBetaId),LjEns%staticBetaList(LjEns%myBetaId-1),LjEns%V,LjEns%toRecieveReals(1), mts)) THEN
            LjEns%toSend(1)=LjEns%betaNeighbours(2)
            LjEns%toSend(2)=LjEns%acceptance+1
            LjEns%acceptance=LjEns%toRecieve(3)
            NSEND3R=LjEns%betaNeighbours(1)
            LjEns%betaNeighbours(2)=LjEns%betaNeighbours(1)
            LjEns%betaNeighbours(1)=LjEns%toRecieve(1)
            LjEns%myBetaId=LjEns%myBetaId-1
          ELSE
            LjEns%toSend(1)=-2
          END IF
          
          CALL MPI_SEND(LjEns%toSend,2,MPI_INTEGER,NDEST1,itag, MPI_COMM_WORLD,IERR)
        END IF
        IF(NDEST3R.LT.numProcesses) then
	  CALL MPI_SEND(NSEND3R,1,MPI_INTEGER, NDEST3R,itag,MPI_COMM_WORLD,IERR)
	end if
      END IF

      IF(MY_IND3.EQ.offset) THEN ! Processes MY_IND3.
        IF(LjEns%betaNeighbours(1).GE.0) THEN
          CALL MPI_RECV(NRECV3R,1,MPI_INTEGER,LjEns%betaNeighbours(1), itag,MPI_COMM_WORLD,STATUS,IERR)
          IF(NRECV3R.NE.-2) LjEns%betaNeighbours(1)=NRECV3R
        END IF
        IF(LjEns%betaNeighbours(2).LT.numProcesses) THEN
          CALL MPI_RECV(NRECV3L,1,MPI_INTEGER,LjEns%betaNeighbours(2), itag,MPI_COMM_WORLD,STATUS,IERR)
          IF(NRECV3L.NE.-2) LjEns%betaNeighbours(2)=NRECV3L
        END IF
      END IF
      
      RETURN
END subroutine tryReplicaSwap

FUNCTION mcCriterionForReplicaSwap( higherBeta, lowerBeta, lowerEnergy, higherEnergy, mts ) 

!    *******************************************************************
!    ** Metropolis criterion for accepting the replica swap move      **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    **								      **
!    *******************************************************************

      use mt_stream
      
      real*8 higherBeta, lowerBeta, lowerEnergy, higherEnergy
      real*8 delta
      Type(mt_state) mts
      
      mcCriterionForReplicaSwap = .false.
      delta = (lowerBeta - higherBeta)*(higherEnergy - lowerEnergy)
      if( genrand_double1(mts) .lt. exp(-delta)) mcCriterionForReplicaSwap = .true.
      RETURN
end FUNCTION mcCriterionForReplicaSwap

SUBROUTINE betaRecursion1(betaArray, BASUM, acceptanceRates, currentBetaIndexes, IPI_B, numProcesses, NCALL, weight, icall, iadd)

!    *******************************************************************
!    ** Recursion for beta temperature distribution		      **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    *******************************************************************

! Copyright Bernd Berg, Jan 8 2002.
      
      !include 'implicit.sta'
      
      real*8, dimension(0:(numProcesses-1)) :: betaArray, BASUM
      integer, dimension(0:(numProcesses-1)) :: acceptanceRates, currentBetaIndexes, IPI_B
      
      integer I, icall, iadd
      integer NACPT_MIN, NACPT_MAX, NCALL
      real*8 XLA_D, XLA, BIM1_NEW, BI_NEW
      
      real*8 weight
     
      IF(ICALL.EQ.0) THEN
        DO I=0, (numProcesses-1)
          BASUM(I)= 0.0
        END DO
      END IF
      ICALL=ICALL+1
      NACPT_MAX=0
      DO I=0,(numProcesses-1)
        IPI_B(currentBetaIndexes(I))=I
        acceptanceRates(I)=2*acceptanceRates(I)
        NACPT_MAX=MAX(NACPT_MAX,acceptanceRates(I))
      END DO
      !IF(NACPT_MAX.EQ.0) CALL STOP_MPI(IUD,MY_ID,'PT_REC1: NACPT_MAX=0')
      NACPT_MIN=NACPT_MAX
      DO I=1,(numProcesses-1)
        NACPT_MIN=MIN(NACPT_MIN,acceptanceRates(IPI_B(I)))
      END DO
      IF(NACPT_MIN.GT.0) IADD=IADD+1
      WEIGHT=WEIGHT+dble(NACPT_MIN)
      DO I=0,(numProcesses-1)
        BASUM(I)=BASUM(I)+dble(NACPT_MIN*betaArray(I)) 
      END DO
      DO I=0,(numProcesses-1)
        IF(acceptanceRates(I).EQ.0) acceptanceRates(I)=1
      END DO
      XLA_D=0.0
      DO I=1,(numProcesses-1)
        XLA_D=XLA_D+dble(acceptanceRates(IPI_B(I)))*(betaArray(I)-betaArray(I-1))
      END DO
      XLA=( betaArray(numProcesses-1) - betaArray(0) )/XLA_D
      BIM1_NEW=betaArray(0)
      DO I=1, (numProcesses-1)
        BI_NEW=BIM1_NEW+(XLA*dble(acceptanceRates(IPI_B(I))))*(betaArray(I)-betaArray(I-1))
        betaArray(I-1)=BIM1_NEW
        BIM1_NEW=BI_NEW
      END DO
      IF(ICALL.EQ.NCALL) THEN
        !IF(IADD.EQ.0) CALL STOP_MPI(IUD,MY_ID,'PT_REC1: IADD=0.')
        DO I=0, (numProcesses-1)
          betaArray(I) = BASUM(I)/WEIGHT
        END DO
      END IF
      RETURN
END subroutine betaRecursion1
