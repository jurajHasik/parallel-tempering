Module LJSystem
	
	use RanGen	

	implicit none

	type LjEnsamble
		Real*8 beta, V, vn, w
		Real*8, allocatable, dimension(:) :: x,y,z
		Integer lowerBetaEns, higherBetaEns
	end type LjEnsamble

	contains

!    		*******************************************************************  
!    		** DOES ONE CYCLE OF TRIAL MOVES OVER ALL ATOMS                  **
!    		*******************************************************************	
		Subroutine DOMOVES ( N, SIGMA, DRMAX, LjEns )

		INTEGER     N
		REAL*8      SIGMA, DRMAX
		TYPE(LjEnsamble) LjEns
		
		REAL*8      RXIOLD, RYIOLD, RZIOLD, RXINEW, RYINEW, RZINEW
		REAL*8      VNEW, VOLD, DELTV, DELTVB
		REAL*8      WNEW, WOLD, DELTW		
		INTEGER	I
		REAL*8, PARAMETER :: DUMM = 1.0
	
!    		*******************************************************************	
	
		DO I = 1, N

              RXIOLD = LjEns%x(I)
              RYIOLD = LjEns%y(I)
              RZIOLD = LjEns%z(I)

!          	** CALCULATE THE ENERGY OF I IN THE OLD CONFIGURATION **

		CALL ENERGY(RXIOLD,RYIOLD,RZIOLD,I,SIGMA,VOLD,WOLD, N, LjEns )
              
!          	** MOVE I **

              	RXINEW = RXIOLD + ( 2.0 * ZBQLU01(DUMM) - 1.0 ) * DRMAX
              	RYINEW = RYIOLD + ( 2.0 * ZBQLU01(DUMM) - 1.0 ) * DRMAX
              	RZINEW = RZIOLD + ( 2.0 * ZBQLU01(DUMM) - 1.0 ) * DRMAX

!          	** CALCULATE THE ENERGY OF I IN THE NEW CONFIGURATION **

              CALL ENERGY(RXINEW,RYINEW,RZINEW,I,SIGMA,VNEW,WNEW, N, LjEns )

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
                 ELSEIF ( EXP ( - DELTVB ) .GT. ZBQLU01(DUMM) ) THEN
                    LjEns%V      = LjEns%V + DELTV
                    LjEns%W      = LjEns%W + DELTW
                    LjEns%X(I)  = RXINEW
                    LjEns%Y(I)  = RYINEW
                    LjEns%Z(I)  = RZINEW
               ENDIF
               
!          *************************************************************
!          ** ENDS LOOP OVER ATOMS                                    **
!          *************************************************************

           ENDDO
           LjEns%vn = LjEns%v / REAL ( N )
           
        RETURN   
        END subroutine domoves 

 	SUBROUTINE ENERGY ( RXI, RYI, RZI, I, SIGMA, V, W, N, LjEns )

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
        END subroutine Energy
        
        
        SUBROUTINE SUMUP ( SIGMA, N, LjEns )
        
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
        ENd subroutine sumup
        
        
        SUBROUTINE TrySwap ( LjEnsambles, I, ensmax, success )
	  
	  Integer I, toSwap, ensmax
	  TYPE(LJENSAMBLE), dimension( ensmax ) :: LJENSAMBLES
	  logical success
	  
	  Real*8 upOrDown
	  Integer tempIndexUp, tempIndexDown
	  Real*8 deltaBeta, deltaV, tempB
	  REAL*8, PARAMETER :: DUMM = 1.0
	  
	  upOrDown = ZBQLU01(DUMM)
	  if ( (upOrDown .LE. 0.5) .and. (LJENSAMBLES(i)%lowerBetaEns .gt. 0) ) then
		 toSwap = LJENSAMBLES(i)%lowerBetaEns
		 deltaBeta = LJENSAMBLES(toSwap)%BETA - LJENSAMBLES(i)%BETA
		 deltaV = LJENSAMBLES(i)%V - LJENSAMBLES(toSwap)%V
	         IF ( deltaV .ge. 0.0 ) THEN
		    tempB = LJEnSAMBLES(i)%BETA
		    LJENSAMBLES(i)%BETA = LJENSAMBLES(toSwap)%BETA
		    LJENSAMBLES(toSwap)%BETA = tempB
		    success = .TRUE.
!		 ** index arithmetics **
		    tempIndexDown = toSwap
		    tempIndexUp = LJENSAMBLES(i)%higherBetaEns
		    LJENSAMBLES(i)%lowerBetaEns = LJENSAMBLES(toSwap)%lowerBetaEns
		    LJENSAMBLES(tempIndexUp)%lowerBetaEns = toSwap 
		    LJENSAMBLES(i)%higherBetaEns = toSwap
		    LJENSAMBLES(LJENSAMBLES(toSwap)%lowerBetaEns)%higherBetaEns = i
		    LJENSAMBLES(toSwap)%lowerBetaEns = i
		    LJENSAMBLES(toSwap)%higherBetaEns = tempIndexUp
		 ELSEIF (EXP ( - deltaBeta * deltaV ) .GT. ZBQLU01(DUMM)) THEN
		    tempB = LJENSAMBLES(i)%BETA
		    LJENSAMBLES(i)%BETA = LJENSAMBLES(toSwap)%BETA
		    LJENSAMBLES(toSwap)%BETA = tempB
		    success = .TRUE.
!		 ** index arithmetics **
		     tempIndexDown = toSwap
		    tempIndexUp = LJENSAMBLES(i)%higherBetaEns
		    LJENSAMBLES(i)%lowerBetaEns = LJENSAMBLES(toSwap)%lowerBetaEns
		    LJENSAMBLES(tempIndexUp)%lowerBetaEns = toSwap 
		    LJENSAMBLES(i)%higherBetaEns = toSwap
		    LJENSAMBLES(LJENSAMBLES(toSwap)%lowerBetaEns)%higherBetaEns = i
		    LJENSAMBLES(toSwap)%lowerBetaEns = i
		    LJENSAMBLES(toSwap)%higherBetaEns = tempIndexUp
                 ENDIF
	  else if (LJENSAMBLES(i)%higherBetaEns .lt. 21 ) then
		 toSwap = LJENSAMBLES(i)%higherBetaEns
		 deltaBeta = LJENSAMBLES(toSwap)%BETA - LJENSAMBLES(i)%BETA
		 deltaV = LJENSAMBLES(i)%V - LJENSAMBLES(toSwap)%V
	         IF ( deltaV .le. 0.0 ) THEN
		    tempB = LJENSAMBLES(i)%BETA
		    LJENSAMBLES(i)%BETA = LJENSAMBLES(toSwap)%BETA
		    LJENSAMBLES(toSwap)%BETA = tempB
!		 ** index arithmetics **
		    tempIndexUp = toSwap
		    tempIndexDown = LJENSAMBLES(i)%lowerBetaEns
		    LJENSAMBLES(i)%higherBetaEns = LJENSAMBLES(toSwap)%higherBetaEns
		    LJENSAMBLES(tempIndexDown)%higherBetaEns = toSwap 
		    LJENSAMBLES(i)%lowerBetaEns = toSwap
		    LJENSAMBLES(LJENSAMBLES(toSwap)%higherBetaEns)%lowerBetaEns = i
		    LJENSAMBLES(toSwap)%higherBetaEns = i
		    LJENSAMBLES(toSwap)%lowerBetaEns = tempIndexDown
		    success = .TRUE.
		 ELSEIF (EXP ( - deltaBeta * deltaV ) .GT. ZBQLU01(DUMM)) THEN
		    tempB = LJENSAMBLES(i)%BETA
		    LJENSAMBLES(i)%BETA = LJENSAMBLES(toSwap)%BETA
		    LJENSAMBLES(toSwap)%BETA = tempB
		    success = .TRUE.
!		 ** index arithmetics **
		    tempIndexUp = toSwap
		    tempIndexDown = LJENSAMBLES(i)%lowerBetaEns
		    LJENSAMBLES(i)%higherBetaEns = LJENSAMBLES(toSwap)%higherBetaEns
		    LJENSAMBLES(tempIndexDown)%higherBetaEns = toSwap 
		    LJENSAMBLES(i)%lowerBetaEns = toSwap
		    LJENSAMBLES(LJENSAMBLES(toSwap)%higherBetaEns)%lowerBetaEns = i
		    LJENSAMBLES(toSwap)%higherBetaEns = i
		    LJENSAMBLES(toSwap)%lowerBetaEns = tempIndexDown
                 ENDIF
	  endif
	  
        end subroutine TrySwap

end module
