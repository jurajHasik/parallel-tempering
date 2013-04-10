Module LJSystem
	
	type LjEnsamble
		Real*8 beta, maxDisplacement
		Real*8 V, vn, w
		Real*8, allocatable, dimension(:) :: x,y,z
		Integer lowerBetaId, myBetaId, higherBetaId
		
		Integer, dimension(3) :: neighbourList
		Real*8, dimension(4) :: higherBetaNeighbourStats  
	end type LjEnsamble

	contains

	Subroutine sweepOverReplica ( N, SIGMA, LjEns )
	
!    	*******************************************************************  
!    	** DOES ONE CYCLE OF TRIAL MOVES OVER ALL ATOMS                  **
!    	*******************************************************************
		
		INTEGER     N
		REAL*8      SIGMA
		TYPE(LjEnsamble) LjEns
		
		REAL*8      RXIOLD, RYIOLD, RZIOLD, RXINEW, RYINEW, RZINEW
		REAL*8      VNEW, VOLD, DELTV, DELTVB
		REAL*8      WNEW, WOLD, DELTW		
		INTEGER	I
		
		Integer numberOfAcceptedMoves
		Real*8 currentAcceptanceRatio
		Real*8, parameter :: targetAcceptanceRatio = 0.2 ! Mountain, Thirumalai
	
!    		*******************************************************************	
		numberOfAcceptedMoves = 0
		currentAcceptanceRatio = 0.0
		
		DO I = 1, N

              RXIOLD = LjEns%x(I)
              RYIOLD = LjEns%y(I)
              RZIOLD = LjEns%z(I)

!          	** CALCULATE THE ENERGY OF I IN THE OLD CONFIGURATION **

		CALL energyOfAtom(RXIOLD,RYIOLD,RZIOLD,I,SIGMA,VOLD,WOLD, N, LjEns )
              
!          	** MOVE I **

              	RXINEW = RXIOLD + ( 2.0*rmafun() - 1.0 ) * LjEns%maxDisplacement
              	RYINEW = RYIOLD + ( 2.0*rmafun() - 1.0 ) * LjEns%maxDisplacement
              	RZINEW = RZIOLD + ( 2.0*rmafun() - 1.0 ) * LjEns%maxDisplacement

!          	** CALCULATE THE ENERGY OF I IN THE NEW CONFIGURATION **

              CALL energyOfAtom(RXINEW,RYINEW,RZINEW,I,SIGMA,VNEW,WNEW, N, LjEns )

!          ** CHECK FOR ACCEPTANCE **

              DELTV  = VNEW - VOLD
              DELTW  = WNEW - WOLD
              DELTVB = LjEns%BETA * DELTV

              IF ( DELTV .LE. 0.0 ) THEN
                    LjEns%V      = LjEns%V + DELTV
                    LjEns%W      = LjEns%W + DELTW
                    LjEns%X(I)  = RXINEW
                    LjEns%Y(I)  = RYINEW
                    LjEns%Z(I)  = RZINEW
                    numberOfAcceptedMoves = numberOfAcceptedMoves + 1
                 ELSEIF ( EXP ( - DELTVB ) .GT. rmafun() ) THEN
                    LjEns%V      = LjEns%V + DELTV
                    LjEns%W      = LjEns%W + DELTW
                    LjEns%X(I)  = RXINEW
                    LjEns%Y(I)  = RYINEW
                    LjEns%Z(I)  = RZINEW
                    numberOfAcceptedMoves = numberOfAcceptedMoves + 1
               ENDIF
               
!          *************************************************************
!          ** ENDS LOOP OVER ATOMS                                    **
!          *************************************************************

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

 	
 	SUBROUTINE energyOfAtom ( RXI, RYI, RZI, I, SIGMA, V, W, N, LjEns )
 	
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
        REAL*8      SIGMA, V, W
	Type(LjEnsamble) LjEns
        REAL*8      RXI, RYI, RZI
	INTEGER     I
        
        REAL*8      SIGSQ, SR2, SR6
        REAL*8      RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ
	INTEGER     J
        
!     ******************************************************************

        SIGSQ  = SIGMA * SIGMA

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

        V = 4.0 * V
        W = 48.0 * W / 3.0

        RETURN
        END subroutine energyOfAtom
        
        
        SUBROUTINE energyOfSystem ( SIGMA, N, LjEns )
        
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
        REAL*8      SIGMA
        Type(LjEnsamble) LjEns

        REAL*8      SIGSQ, RXIJ, RYIJ, RZIJ
        REAL*8      RXI, RYI, RZI, VIJ, WIJ, SR2, SR6, RIJSQ
        INTEGER     I, j

!    *******************************************************************

        SIGSQ  = SIGMA * SIGMA

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

        LjEns%V = 4.0 * LjEns%V
        LjEns%W = 48.0 * LjEns%W / 3.0
        
        RETURN
        ENd subroutine energyOfSystem
        
        
        SUBROUTINE InitLjSystemCoords ( N, LjEns )
        
        
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
        
        Integer I 

!    *******************************************************************
	
	Allocate( LjEns%X(N), LjEns%Y(N), LjEns%Z(N))
	Do I=1, N
	  LjEns%X(i)= rmafun()-0.5
	  LjEns%Y(i)= rmafun()-0.5
	  LjEns%Z(i)= rmafun()-0.5
	enddo
        
        RETURN
        END SUBROUTINE InitLjSystemCoords
        

	SUBROUTINE InitLjSystemBeta ( LjEns, Beta, higherBetaId, myBetaId, lowerBetaId, intialMaxDisplacement )
        
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
        Real*8 Beta, intialMaxDisplacement
        Integer lowerBetaId, myBetaId, higherBetaId
        
        Integer I 

!    *******************************************************************
	
	LjEns%Beta = Beta
	LjEns%maxDisplacement = intialMaxDisplacement
	LjEns%lowerBetaId = lowerBetaId
	LjEns%myBetaId = myBetaId
	LjEns%higherBetaId = higherBetaId 
        
        RETURN
        END subroutine InitLjSystemBeta
        
        include 'rmafun.f'
        
end module
