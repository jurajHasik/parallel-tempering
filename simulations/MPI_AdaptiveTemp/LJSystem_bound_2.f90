Module LJSystem

	use mt_stream
	
	type LjEnsamble
		Real*8 maxDisplacement
		Real*8 boundary
		
		Real*8 sigma, LJepsilon
		Real*8 V, vn, w
		Real*8, allocatable, dimension(:) :: x,y,z

		Real*8, allocatable, dimension(:) :: staticBetaList
		integer, dimension(2) :: betaNeighbours
		integer myBetaId
		integer acceptance
		
		integer, dimension(3) :: toSend, toRecieve
		real*8, dimension(4) :: toSendReals, toRecieveReals
		
		real*8, allocatable, dimension(:) :: testingKelvinTemp
		real*8 runningVavg
	end type LjEnsamble

	contains

	Subroutine centreOfMass ( N, LjEns, x, y, z )

		Integer N
		Type(LjEnsamble) LjEns
		Real*8 x, y, z

		Real*8 sumX, sumY, sumZ
		Integer I
		
		sumX=0.0
		sumY=0.0
		sumZ=0.0
		Do i = 1, N
			sumX=sumX+LjEns%x(i)
			sumY=sumY+LjEns%y(i)
			sumZ=sumZ+LjEns%z(i)		
		Enddo

		x = sumX/dble(N)
		y = sumY/dble(N)
		z = sumZ/dble(N)

	RETURN   
        END subroutine centreOfMass

	Subroutine sweepOverReplica ( N, LjEns, mts)
	
!    	*******************************************************************  
!    	** DOES ONE CYCLE OF TRIAL MOVES OVER ALL ATOMS                  **
!    	*******************************************************************
		
		INTEGER     N, atom
		TYPE(LjEnsamble) LjEns
		Type(mt_state) mts
		
		REAL*8      RXIOLD, RYIOLD, RZIOLD, RXINEW, RYINEW, RZINEW, rnewsq
		REAL*8      VNEW, VOLD, DELTV, DELTVB
		REAL*8      WNEW, WOLD, DELTW		
		INTEGER		I
		logical	outOfBounds		
		
		Integer numberOfAcceptedMoves
		Real*8 currentAcceptanceRatio
		Real*8, parameter :: targetAcceptanceRatio = 0.5 ! 0.2 Mountain, Thirumalai
	
!    		*******************************************************************	
		numberOfAcceptedMoves = 0
		currentAcceptanceRatio = 0.0
		
		DO I = 1, N		
		RXIOLD = LjEns%x(i)
		RYIOLD = LjEns%y(i)
		RZIOLD = LjEns%z(i)

!          	** CALCULATE THE ENERGY OF I IN THE OLD CONFIGURATION **

		CALL energyOfAtom(RXIOLD,RYIOLD,RZIOLD,i,VOLD,WOLD, N, LjEns )
              
!          	** MOVE I **

              	RXINEW = RXIOLD + ( genrand_double1(mts) - 0.5 ) * LjEns%maxDisplacement
              	RYINEW = RYIOLD + ( genrand_double1(mts) - 0.5 ) * LjEns%maxDisplacement
              	RZINEW = RZIOLD + ( genrand_double1(mts) - 0.5 ) * LjEns%maxDisplacement

		rnewsq = rxinew*rxinew + ryinew*ryinew + rzinew*rzinew

		if ( rnewsq .lt. LjEns%boundary*LjEns%boundary ) then
!          	** CALCULATE THE ENERGY OF I IN THE NEW CONFIGURATION **

		CALL energyOfAtom(RXINEW,RYINEW,RZINEW,i,VNEW,WNEW, N, LjEns )

!          	** CHECK FOR ACCEPTANCE **

		DELTV  = VNEW - VOLD
		DELTW  = WNEW - WOLD
		DELTVB = LjEns%staticBetaList(LjEns%myBetaId) * DELTV

		IF ( DELTV .LE. 0.0 ) THEN
                    LjEns%V      = LjEns%V + DELTV
                    LjEns%W      = LjEns%W + DELTW
                    LjEns%X(i)  = RXINEW
                    LjEns%Y(i)  = RYINEW
                    LjEns%Z(i)  = RZINEW
                    numberOfAcceptedMoves = numberOfAcceptedMoves + 1
                 ELSEIF ( EXP ( - DELTVB ) .GT. genrand_double1(mts) ) THEN
                    LjEns%V      = LjEns%V + DELTV
                    LjEns%W      = LjEns%W + DELTW
                    LjEns%X(i)  = RXINEW
                    LjEns%Y(i)  = RYINEW
                    LjEns%Z(i)  = RZINEW
                    numberOfAcceptedMoves = numberOfAcceptedMoves + 1
               ENDIF
               
!          *************************************************************
!          ** ENDS LOOP OVER ATOMS                                    **
!          *************************************************************
		endif
           ENDDO
        LjEns%vn = LjEns%v / dble ( N )
           
!	** Optimization of maximal Displacement **

	currentAcceptanceRatio = dble(numberOfAcceptedMoves) / dble(N)
        if (currentAcceptanceRatio .gt. targetAcceptanceRatio) then
	  LjEns%maxDisplacement = LjEns%maxDisplacement*1.05
	else if (currentAcceptanceRatio .lt. targetAcceptanceRatio) then
	  LjEns%maxDisplacement = LjEns%maxDisplacement*0.95
        end if
        
        RETURN   
        END subroutine sweepOverReplica
        
        
        Subroutine sweepOverReplicaWOAccept ( N, LjEns, mts)
	
!    	*******************************************************************  
!    	** DOES ONE CYCLE OF TRIAL MOVES OVER ALL ATOMS                  **
!    	*******************************************************************
		
		INTEGER     N, atom
		TYPE(LjEnsamble) LjEns
		Type(mt_state) mts
		
		REAL*8      RXIOLD, RYIOLD, RZIOLD, RXINEW, RYINEW, RZINEW, rnewsq
		REAL*8      VNEW, VOLD, DELTV, DELTVB
		REAL*8      WNEW, WOLD, DELTW		
		INTEGER		I	
	
!    		*******************************************************************	
		
		DO I = 1, N		
		RXIOLD = LjEns%x(i)
		RYIOLD = LjEns%y(i)
		RZIOLD = LjEns%z(i)

!          	** CALCULATE THE ENERGY OF I IN THE OLD CONFIGURATION **

		CALL energyOfAtom(RXIOLD,RYIOLD,RZIOLD,i,VOLD,WOLD, N, LjEns )
              
!          	** MOVE I **

              	RXINEW = RXIOLD + ( genrand_double1(mts) - 0.5 ) * LjEns%maxDisplacement
              	RYINEW = RYIOLD + ( genrand_double1(mts) - 0.5 ) * LjEns%maxDisplacement
              	RZINEW = RZIOLD + ( genrand_double1(mts) - 0.5 ) * LjEns%maxDisplacement

		rnewsq = rxinew*rxinew + ryinew*ryinew + rzinew*rzinew

		if ( rnewsq .lt. LjEns%boundary*LjEns%boundary ) then
!          	** CALCULATE THE ENERGY OF I IN THE NEW CONFIGURATION **

		CALL energyOfAtom(RXINEW,RYINEW,RZINEW,i,VNEW,WNEW, N, LjEns )

!          	** CHECK FOR ACCEPTANCE **

		DELTV  = VNEW - VOLD
		DELTW  = WNEW - WOLD
		DELTVB = LjEns%staticBetaList(LjEns%myBetaId) * DELTV

		IF ( DELTV .LE. 0.0 ) THEN
                    LjEns%V      = LjEns%V + DELTV
                    LjEns%W      = LjEns%W + DELTW
                    LjEns%X(i)  = RXINEW
                    LjEns%Y(i)  = RYINEW
                    LjEns%Z(i)  = RZINEW
                 ELSEIF ( EXP ( - DELTVB ) .GT. genrand_double1(mts) ) THEN
                    LjEns%V      = LjEns%V + DELTV
                    LjEns%W      = LjEns%W + DELTW
                    LjEns%X(i)  = RXINEW
                    LjEns%Y(i)  = RYINEW
                    LjEns%Z(i)  = RZINEW
               ENDIF
               
!          *************************************************************
!          ** ENDS LOOP OVER ATOMS                                    **
!          *************************************************************
		endif
           ENDDO
        LjEns%vn = LjEns%v / dble ( N )
           
        RETURN   
        END subroutine sweepOverReplicaWOAccept 

 	
 	SUBROUTINE energyOfAtom ( RXI, RYI, RZI, I, V, W, N, LjEns )
 	
!    	*******************************************************************
!    	** RETURNS THE POTENTIAL ENERGY OF ATOM I WITH ALL OTHER ATOMS.  **
!    	**                                                               **
!    	** PRINCIPAL VARIABLES:                                          **
!    	**                                                               **
!    	** INTEGER I                 THE ATOM OF INTEREST                **
!    	** INTEGER N                 THE NUMBER OF ATOMS                 **
!    	** REAL    RX(N),RY(N),RZ(N) THE ATOM POSITIONS                  **
!    	** REAL    RXI,RYI,RZI       THE COORDINATES OF ATOM I           **
!    	** REAL    V                 THE POTENTIAL ENERGY OF ATOM I      **
!    	** REAL    W                 THE VIRIAL OF ATOM I                **
!    	**                                                               **
!    	** USAGE:                                                        **
!    	**                                                               **
!    	** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
!    	** DURING A TRIAL MOVE OF ATOM I. IT IS CALLED BEFORE AND        **
!    	** AFTER THE RANDOM DISPLACEMENT OF I.                           **
!    	*******************************************************************

	INTEGER     N
        REAL*8      V, W
	Type(LjEnsamble) LjEns
        REAL*8      RXI, RYI, RZI
	INTEGER     I   

        REAL*8      SIGSQ, SR2, SR6
        REAL*8      RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ
	INTEGER     J
        
!     ******************************************************************

        SIGSQ  = LjEns%SIGMA * LjEns%SIGMA

        V      = 0.0
        W      = 0.0

!    ** LOOP OVER ALL MOLECULES EXCEPT I  **

        DO J = 1, N
           IF ( I .NE. J ) THEN
              RXIJ  = RXI - LjEns%X(J)
              RYIJ  = RYI - LjEns%Y(J)
              RZIJ  = RZI - LjEns%Z(J)
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
	   
              SR2 = SIGSQ / RIJSQ
              SR6 = SR2 * SR2 * SR2
              VIJ = SR6 * ( SR6 - 1.0 )
              WIJ = SR6 * ( SR6 - 0.5 )
              V = V + VIJ
              W = W + WIJ

           ENDIF
        ENDDO

        V = 4.0 * LjEns%LJepsilon * V
        W = 48.0 * W / 3.0

        RETURN
        END subroutine energyOfAtom
        
        
        SUBROUTINE energyOfSystem ( N, LjEns )
        
!    *******************************************************************
!    ** CALCULATES THE TOTAL POTENTIAL ENERGY FOR A CONFIGURATION.    **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                 THE NUMBER OF ATOMS                 **
!    ** REAL    RX(N(,RY(N),RZ(N) THE POSITIONS OF THE ATOMS          **
!    ** REAL    V                 THE POTENTIAL ENERGY                **
!    ** REAL    W                 THE VIRIAL                          **
!    ** LOGICAL OVRLAP            TRUE FOR SUBSTANTIAL ATOM OVERLAP   **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THE SUBROUTINE RETURNS THE TOTAL POTENTIAL ENERGY AT THE      **
!    ** BEGINNING AND END OF THE RUN.                                 **
!    *******************************************************************

	INTEGER     N
        Type(LjEnsamble) LjEns

        REAL*8      SIGSQ, RXIJ, RYIJ, RZIJ
        REAL*8      RXI, RYI, RZI, VIJ, WIJ, SR2, SR6, RIJSQ
        INTEGER     I, j

!    *******************************************************************

        SIGSQ  = LjEns%SIGMA * LjEns%SIGMA

        LjEns%V      = 0.0
        LjEns%W      = 0.0

!    ** LOOP OVER ALL THE PAIRS IN THE LIQUID **
        DO I = 1, N - 1
           RXI = LjEns%x(I)
           RYI = LjEns%Y(I)
           RZI = LjEns%Z(I)
           DO j = I + 1, N
              RXIJ  = RXI - LjEns%X(j)
              RYIJ  = RYI - LjEns%Y(j)
              RZIJ  = RZI - LjEns%Z(j)
              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
              SR2 = SIGSQ / RIJSQ
              SR6 = SR2 * SR2 * SR2
              VIJ = SR6 * ( SR6 - 1.0 )
              WIJ = SR6 * ( SR6 - 0.5 )
              LjEns%V   = LjEns%V + VIJ
              LjEns%W   = LjEns%W + WIJ
           ENDDO
        ENDDO

        LjEns%V = 4.0 * LjEns%LJepsilon * LjEns%V
        LjEns%W = 48.0 * LjEns%W / 3.0
        
        RETURN
        ENd subroutine energyOfSystem
        
        
        SUBROUTINE InitLjSystemCoords ( N, LjEns, mts )
        
        
!    *******************************************************************
!    ** Intializes LJ system of N atoms to random CONFIGURATION.      **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                 THE NUMBER OF ATOMS                 **
!    **								      **
!    ** THE SUBROUTINE RETURNS LJ system in random configuration      **
!    *******************************************************************

	INTEGER     N
        Type(LjEnsamble) LjEns
        Type(mt_state) mts
        
        Integer I
        Real*8 diameter

!    *******************************************************************
	
	diameter = 3.5 * LjEns%sigma !* (dble(N)**(1.0d0/3.0d0))
	LjEns%boundary = 4.0 * LjEns%sigma
	Allocate( LjEns%X(N), LjEns%Y(N), LjEns%Z(N))
	Do I=1, N
	  LjEns%X(i)= (genrand_double1(mts)-0.5)*diameter
	  LjEns%Y(i)= (genrand_double1(mts)-0.5)*diameter
	  LjEns%Z(i)= (genrand_double1(mts)-0.5)*diameter
	enddo
        
        RETURN
        END SUBROUTINE InitLjSystemCoords
        

	SUBROUTINE InitLjSystemBeta ( LjEns, betaArray, numOfEnsambles, myId, intialMaxDisplacement, sigma, LJepsilon )
        
!    *******************************************************************
!    ** Intializes LJ system of N atoms to random CONFIGURATION.      **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                 THE NUMBER OF ATOMS                 **
!    **								      **
!    ** THE SUBROUTINE RETURNS LJ system in random configuration      **
!    *******************************************************************

	Type(LjEnsamble) LjEns
        Real*8 intialMaxDisplacement, sigma, LJepsilon	  
	integer myId, numOfEnsambles
	real*8, dimension(0:(numOfEnsambles-1)) :: betaArray
	        
!    *******************************************************************
	
	LjEns%sigma = sigma
	LjEns%LJepsilon = LJepsilon
	LjEns%maxDisplacement = intialMaxDisplacement
        allocate(LjEns%staticBetaList(0:(numOfEnsambles-1)))
        LjEns%staticBetaList = betaArray
        LjEns%myBetaId = myId
        LjEns%betaNeighbours(1) = myId-1
        LjEns%betaNeighbours(2) = myId+1
        LjEns%acceptance = 0
        
        LjEns%toSend(1) = LjEns%myBetaId
        LjEns%toSend(2) = LjEns%betaNeighbours(1)
        LjEns%toSend(3) = LjEns%betaNeighbours(2)
        
        RETURN
        END subroutine InitLjSystemBeta
        
end module
