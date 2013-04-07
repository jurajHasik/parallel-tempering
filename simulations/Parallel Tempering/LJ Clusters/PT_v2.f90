! -*-fortran-*-

! *****************************************************************************
! ** FICHE F.11.  CONSTANT-NVT MONTE CARLO FOR LENNARD JONES ATOMS           **
! ** This FORTRAN code is intended to illustrate points made in the text.    **
! ** To our knowledge it works correctly.  However it is the responsibility  **
! ** of the user to test it, if it is to be used in a research application.  **
! *****************************************************************************

! Modified by J.Hasik, Dec 2012
! Continous modifications 8.1.2013

! Optimization by Simulated Annealing, S Kirckpatrick, C. D. Gellaat; M.P. Vecchi
! Science, New Series, Vol. 220, No. 4598. (May 13, 1983)

! Parallel Tempering with simulated Annealing 20.1.2013
  

        PROGRAM PT
	
	use RanGen
	use LjSystem
	
!    *******************************************************************
!    ** MONTE CARLO SIMULATION PROGRAM IN THE CONSTANT-NVT ENSEMBLE.  **
!    **                                                               **
!    ** THIS PROGRAM TAKES A CONFIGURATION OF LENNARD JONES ATOMS     **
!    ** AND PERFORMS A CONVENTIONAL NVT MC SIMULATION. THE BOX IS OF  **
!    ** UNIT LENGTH, -0.5 TO +0.5 AND THERE ARE NO LOOKUP TABLES.     **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                   NUMBER OF MOLECULES               **
!    ** INTEGER NSTEP               MAXIMUM NUMBER OF CYCLES /unused/ **
!    ** INTEGER NEQU                NUMBER OF EQUILIBRATION CYCLES    **
!    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
!    ** REAL    TEMP                REDUCED TEMPERATURE               **
!    ** REAL    SIGMA               REDUCED LJ DIAMETER               **
!    ** REAL    RMIN                MINIMUM REDUCED PAIR SEPARATION   **
!    ** REAL    RCUT                REDUCED CUTOFF DISTANCE           **
!    ** REAL    DRMAX               REDUCED MAXIMUM DISPLACEMENT      **
!    ** REAL    V                   THE POTENTIAL ENERGY              **
!    ** REAL    W                   THE VIRIAL                        **
!    ** REAL    PRES                THE PRESSURE                      **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THE PROGRAM TAKES IN A CONFIGURATION OF ATOMS                 **
!    ** AND RUNS A MONTE CARLO SIMULATION AT THE GIVEN TEMPERATURE    **
!    ** FOR THE SPECIFIED NUMBER OF CYCLES.                           **
!    **                                                               **
!    ** UNITS:                                                        **
!    **                                                               **
!    ** THE PROGRAM USES LENNARD-JONES UNITS FOR USER INPUT AND       **
!    ** OUTPUT BUT CONDUCTS THE SIMULATION IN A BOX OF UNIT LENGTH.   **
!    ** FOR EXAMPLE, FOR A BOXLENGTH L, AND LENNARD-JONES PARAMETERS  **
!    ** EPSILON AND SIGMA, THE UNITS ARE:                             **
!    **                                                               **
!    **     PROPERTY       LJ  UNITS            PROGRAM UNITS         **
!    **                                                               **
!    **     TEMP           EPSILON/K            EPSILON/K             **
!    **     PRES           EPSILON/SIGMA**3     EPSILON/L**3          **
!    **     V              EPSILON              EPSILON               **
!    **     DENS           1/SIGMA**3           1/L**3                **
!    **                                                               **
!    ** ROUTINES REFERENCED (AND INCLUDED IN THIS FILE):              **
!    **                                                               **
!    ** SUBROUTINE SUMUP ( RCUT, RMIN, SIGMA, OVRLAP, V, W )          **
!    **    CALCULATES THE TOTAL POTENTIAL ENERGY FOR A CONFIGURATION  **
!    ** SUBROUTINE ENERGY ( RXI, RYI, RZI, I, RCUT, SIGMA, V, W )     **
!    **    CALCULATES THE POTENTIAL ENERGY OF ATOM I WITH ALL THE     **
!    **    OTHER ATOMS IN THE LIQUID                                  **
!    ** SUBROUTINE READCN (CNFILE )                                   **
!    **    READS IN A CONFIGURATION                                   **
!    ** SUBROUTINE WRITCN ( CNFILE )                                  **
!    **    WRITES OUT A CONFIGURATION                                 **
!    ** REAL FUNCTION RANF ( SEED )                                   **
!    **    RETURNS A UNIFORM RANDOM NUMBER BETWEEN ZERO AND ONE       **
!    ** SUBROUTINE ORDER ( KX, KY, KZ, RHO )                          **
!    **    RETURNS THE ORDER PARAMETER RHO FOR K-VECTOR (KX,KY,KZ)    **
!    **                                                               **
!    *******************************************************************

        REAL*8        DRMAX, DENS, TEMP, SIGMA
        REAL*8        V, VNEW, VOLD, VEND, VN
        REAL*8        W, WEND, WNEW, WOLD, PRES
        INTEGER       SwitchAttempts, numSwitch
        INTEGER       SEED
        REAL*8        TSTEP
        Logical	      swapSuccess
        
!       Minimum temperature reached during simulation in LJ units
        REAL*8, Parameter :: finalTemp = 0.039777
!       Attempt ensamble switch every -this- moves        
        INTEGER, Parameter :: PTSwitch = 100
        
	INTEGER     NEQU, ANNEAL
	
	INTEGER     N
	REAL*8, ALLOCATABLE, dimension(:) :: RX, RY, RZ
	
	INTEGER     ENSMAX, ENS
	PARAMETER   ( ENSMAX = 20 ) ! Maximal number of ensambles 
	TYPE(LJENSAMBLE), dimension( ENSMAX ) :: LJENSAMBLES
	
        INTEGER     STEP, j, k, IPRINT, ISAVE
        CHARACTER   TITLE*80, CNFILE*30
        CHARACTER   SAVEFILE*30
        CHARACTER   energyData*30
        
        REAL*8, PARAMETER :: PI = 3.1415927
	REAL*8, PARAMETER :: DUMM = 1.0 

!       ****************************************************************

!    ** READ INPUT DATA **

        WRITE(*,'(//,'' **** PROGRAM MCLJ ****                   ''/)')
        WRITE(*,'('' CONSTANT-NVT MONTE CARLO PROGRAM            '' )')
        WRITE(*,'('' FOR LENNARD JONES ATOMS                      '')')

        WRITE(*,'('' ENTER THE RUN TITLE                          '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER Number of Atoms '')')
	READ (*,*) N
	WRITE(*,'('' ENTER Starting Temperature in Lj units '')')
	READ (*,*) Temp
        WRITE(*,'('' ENTER NUMBER OF EQUILIBRATION CYCLES         '')')
        READ (*,*) NEQU
        WRITE(*,'('' ENTER NUMBER OF Parallel Tempering CYCLES    '')')
        READ (*,*) ANNEAL
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN OUTPUT LINES   '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER NUMBER OF STEPS BETWEEN DATA SAVES     '')')
        READ (*,*) ISAVE
        WRITE(*,'('' ENTER THE CONFIGURATION FILE NAME            '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER THE Configuration SAVE FILE NAME       '')')
        READ (*,'(A)') SAVEFILE
        WRITE(*,'('' ENTER THE Energy data Save FILE NAME         '')')
        READ (*,'(A)') energyData
        WRITE(*,'('' RANDOM NUMBER GENERATOR SEED                 '')')
        READ (*,*) SEED

!    ** WRITE INPUT DATA **

        WRITE(*,'(       //1X                    ,A     )') TITLE
        WRITE(*,'('' OUTPUT FREQUENCY          '',I10   )') IPRINT
        WRITE(*,'('' SAVE FREQUENCY            '',I10   )') ISAVE
        WRITE(*,'('' RANDOM NUMBER GEN. SEED   '',I10   )') SEED
        WRITE(*,'('' CONFIGURATION FILE  NAME  '',A     )') CNFILE
        WRITE(*,'('' TEMPERATURE               '',F10.4 )') TEMP
        WRITE(*,'('' PT Cycles               '',F10.4 )') ANNEAL

!    ** READ INITIAL CONFIGURATION **
	ALLOCATE(RX(N),RY(N),RZ(N))
        CALL READCN ( CNFILE, N, RX, RY, RZ )
        
        WRITE(*,'('' Number of molecules     '',F10.4 )') N
	Do j=1, N
	  write(*,'(" 1 ",F10.4)') RX(j), Ry(j), Rz(j)
	enddo
        
!    ** Copy starting configurations to each ensamble **
        DO j=1, ENSMAX
	  ALLOCATE( LJENSAMBLES(j)%X(N), LJENSAMBLES(j)%Y(N), LJENSAMBLES(j)%Z(N) )
	  LJENSAMBLES(j)%X = RX
	  LJENSAMBLES(j)%Y = RY
	  LJENSAMBLES(j)%Z = RZ
        ENDDO
        
!       John A. White "Lennard-Jones as a model for argon and test of extended renormalization group calculations", 
!	Journal of Chemical Physics 111 pp. 9352-9356 (1999)
!	Parameter Set #4
!	epsilon/K_b = 125.7 K
!	sigma = 3,4345 Ang


        SIGMA  = 3.4345
        DRMAX  = 0.15 * SIGMA

!    ** WRITE OUT SOME USEFUL INFORMATION **

        WRITE(*,'('' SIGMA              =  '',F10.4)')  SIGMA

!    ** CALCULATE INITIAL ENERGY. Set V, VN, W, Beta for each ensamble **        
	DO j=1, ENSMAX
	   CALL SUMUP ( SIGMA, N, LJENSAMBLES(j) )
	   LJENSAMBLES(j)%VN = LJENSAMBLES(j)%V / REAL ( N )
        ENDDO

!    ** Initialize temperatures across ensambles ** 
        TSTEP = (finalTemp/TEMP)**(1.0/(REAL(ENSMAX)-1.0))
        
        Write(*, '('' Geometric series Temperature coefficient  '', F10.6)') TSTEP
        
        DO j=1, ENSMAX
	   LJENSAMBLES(j)%BETA = 1.0 / (TEMP*(TSTEP**(REAL(j)-1.0)))	   
	   LJENSAMBLEs(j)%lowerBetaEns = j-1
	   LJENSAMBLEs(j)%higherBetaEns = j+1
        ENDDO
        
        WRITE(*,'('' INITIAL VN for each ensamble =  '', F10.4 )' ) LJENSAMBLES(1)%VN
        
        WRITE(*,'(//'' START OF PARALLEL TEMPERING          ''//)')
        WRITE(*,'(''    STEP      V/N       P       TEMP ''/)')
        
!    ** INITIALIZE RANDOM NUMBER GENERATOR **        
	CALL ZBQLINI ( SEED )

!    *******************************************************************
!    ** LOOPS OVER ALL ENSAMBLES                                      **
!    *******************************************************************
        DO STEP = 1, NEQU
	  DO j = 1, ENSMAX
              CALL DOMOVES ( N, SIGMA, DRMAX, LjEnsambles(j) )
	  EndDo
	  
!	  ** PERFORM PERIODIC OPERATIONS  **
	      IF ( MOD ( STEP, IPRINT ) .EQ. 0 ) THEN
!         ** WRITE OUT RUNTIME INFORMATION **
		WRITE(*,'(1I8 , 3F10.4)') STEP, LJensambles(20)%Vn, Pres, Temp
	      ENDIF
	      IF ( MOD ( STEP, ISAVE ) .EQ. 0 ) THEN
!         ** WRITE OUT THE CONFIGURATION AT INTERVALS **
		CALL WRITEenergy ( energyData, ENSMAX, LJENSAMBLES )
	      ENDIF
        ENDDO
           
        WRITE(6,'(20X,''EQUILIBRATION FINISHED '',I10)')
        
        Write(*, '('' Ensamble    Temp    energy    lowerBeta     higherBeta'')')
        DO j=1, ENSMAX
	   CALL SUMUP ( SIGMA, N, LJENSAMBLES(j) )
	   LJENSAMBLES(j)%VN = LJENSAMBLES(j)%V / REAL ( N )
	   Write(*, '(1I10, 2F10.4, 2I10 )') j, LJENSAMBLES(j)%BETA, LjEnsambles(j)%vn, LJENSAMBLEs(j)%lowerBetaEns, LJensambles(j)%higherBetaEns
        ENDDO
   
!    *******************************************************************
!    ** PARALLEL TAMPERING STARTED                                    **
!    *******************************************************************	
	numSwitch = 0
	SwitchAttempts = 0
		
	DO STEP = NEQU, ANNEAL
	   
	   DO j = 1, ENSMAX
              CALL DOMOVES ( N, SIGMA, DRMAX, LjEnsambles(j) )
           ENDDO

!          ** PERFORM PERIODIC OPERATIONS  **
	   IF ( MOD (STEP, PTSwitch ) .EQ. 0 ) THEN
	      DO k = 1, ENSMAX
!		 deltaB = LJENSAMBLES(k+1)%BETA - LJENSAMBLES(k)%BETA
!		 deltaV = LJENSAMBLES(k)%V - LJENSAMBLES(k+1)%V
!	         IF ( LJENSAMBLES(k+1)%V .LE. LJENSAMBLES(k)%V ) THEN
!		    tempB = LJENSAMBLES(k)%BETA
!		    LJENSAMBLES(k)%BETA = LJENSAMBLES(k+1)%BETA
!		    LJENSAMBLES(k+1)%BETA = tempB
!		    numSwitch = numSwitch + 1
!		 ELSEIF (EXP ( deltaB * deltaV ) .GT. ZBQLU01(DUMM)) THEN
!		    tempB = LJENSAMBLES(k)%BETA
!		    LJENSAMBLES(k)%BETA = LJENSAMBLES(k+1)%BETA
!		    LJENSAMBLES(k+1)%BETA = tempB
!		    numSwitch = numSwitch + 1
!                 ENDIF
		SwitchAttempts = SwitchAttempts + 1
		Call TrySwap (LjEnsambles, k, ensmax, swapSuccess )
		If ( swapSuccess .eq. .TRUE. ) then
		  numSwitch = numSwitch + 1
		endif
		swapSuccess = .False.
              ENDDO
	   ENDIF
           
           IF ( MOD ( STEP, IPRINT ) .EQ. 0 ) THEN
!          ** WRITE OUT RUNTIME INFORMATION **
              WRITE(*,'(1I8, 3F10.4)') STEP, LJensambles(20)%Vn, Pres, Temp
           ENDIF
           IF ( MOD ( STEP, ISAVE ) .EQ. 0 ) THEN
!          ** WRITE OUT THE CONFIGURATION AT INTERVALS **
              CALL WRITEenergy ( energyData, ENSMAX, LJensambles )
           ENDIF
        ENDDO
           	   	
!    *******************************************************************
!    ** ENDS THE LOOP OVER CYCLES                                     **
!    *******************************************************************   	  	   
        
	WRITE(*,'(//'' END OF PARALLEL TEMPERING         ''//)')
	WRITE(*,'("Number of attempted swap moves         ",I10)') SwitchAttempts
	WRITE(*,'("Number of successful swap moves         ",I10)') numSwitch
	
!    ** CHECKS FINAL VALUE OF THE POTENTIAL ENERGY IS CONSISTENT **

        DO j=1, ENSMAX
	   CALL SUMUP ( SIGMA, N, LjEnsambles(j) )
	   LJENSAMBLES(j)%VN = LJENSAMBLES(j)%V / REAL ( N )
        ENDDO

	WRITE(*, '('' V/N at Temp for ensamble -array index- lowerBeta    higherBeta'',I10)')
        DO j=1, ENSMAX
	   tempB = 1 / LJENSAMBLES(j)%BETA
	   WRITE(*,'(2F12.6,3I10)') LJENSAMBLES(j)%VN, tempB, j, LJENSAMBLES(j)%lowerBetaEns, LJENSAMBLES(j)%higherBetaEns
        ENDDO

!    ** WRITE OUT THE FINAL CONFIGURATION FROM THE RUN **

        DO j=1, ENSMAX
           CALL WRITCN ( SAVEFILE, N, LJENSAMBLES(j)%VN, LJENSAMBLES(j)%X, &
	                 & LJENSAMBLES(j)%Y, LJENSAMBLES(j)%Z )
        ENDDO

        WRITE(*,'(/'' END OF SIMULATION '')')

        STOP
        END

        
        SUBROUTINE READCN ( CNFILE, N, RX, RY, RZ )

!    *******************************************************************
!    ** SUBROUTINE TO READ IN THE CONFIGURATION FROM UNIT 10          **
!    *******************************************************************

	INTEGER	    N
	REAL*8,     dimension(N) :: RX, RY, RZ
        CHARACTER   CNFILE*(*)
        
        INTEGER, PARAMETER :: CNUNIT = 10
        
!   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'OLD')


        READ ( CNUNIT,* ) N
        DO I=1, N
           READ ( CNUNIT,* ) RX(I), RY(I), RZ(I)
        ENDDO

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( SAVEFILE, N, VN, tempX, tempY, tempZ )
        
!    *******************************************************************
!    ** SUBROUTINE TO WRITE OUT THE CONFIGURATION TO UNIT 10          **
!    *******************************************************************

	INTEGER	    N
	REAL*8,     dimension(N) :: tempX, tempY, tempZ
	REAL*8      VN
	CHARACTER   SAVEFILE*(*)
        
        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 11 )
        INTEGER     I

!   ********************************************************************

	OPEN ( UNIT = CNUNIT, FILE = SAVEFILE, STATUS = 'UNKNOWN', POSITION = 'APPEND')

        WRITE ( CNUNIT, * ) N
	WRITE ( CNUNIT, '(''V/N '',1F10.6)') VN
	DO I=1, N
           WRITE ( CNUNIT, * ) 'Ar ', tempX(I), tempY(I), tempZ(I)
        ENDDO

        CLOSE ( UNIT = CNUNIT )
        
        RETURN
        END
        
        
        
        SUBROUTINE WRITEenergy ( SAVEFILE, ENS, LjEnsambles )
        
!    *******************************************************************
!    ** SUBROUTINE TO WRITE OUT THE CONFIGURATION TO UNIT 10          **
!    *******************************************************************

	use LjSystem

	INTEGER	    ENS
	TYPE(LJENSAMBLE), dimension( ENS ) :: LJENSAMBLES
	CHARACTER   SAVEFILE*(*)
        
        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 12 )
        INTEGER     I

!   ********************************************************************

	OPEN ( UNIT = CNUNIT, FILE = SAVEFILE, STATUS = 'UNKNOWN', POSITION = 'APPEND')

	DO I=1, ENS
           WRITE ( CNUNIT, '(F10.4)', advance='no' ) LjEnsambles(i)%Vn
        ENDDO
        write( CNunit, *) ''

        CLOSE ( UNIT = CNUNIT )
        
        RETURN
        END





 

 
