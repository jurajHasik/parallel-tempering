      SUBROUTINE RMASET(iuo,iud,iseed1,iseed2,cfile)
! Copyright Bernd Berg, Sep 21, 2000.
! INITIALIZING ROUTINE FOR RANMAR OR RAVMAR, MUST BE CALLED 
! BEFORE GENERATING PSEUDORANDOM NUMBERS WITH RANMAR (RAVMAR).
! RANGES: 0 <= IJ <= 31328  AND  0 <= KL <= 30081.
! FOR IUO<=5 THE INFORMATIVE MESSAGES ARE OFF.
      include 'implicit.sta'
      include 'constants.par'
      character*(*) cfile
      logical lexist
      COMMON/RASET1/U(97),C,CD,CM,I,J
!
      inquire(file=cfile,exist=lexist)
      IF(.not.lexist) THEN
        IJ=iseed1+1801
        KL=9373+iseed2 
        I=MOD(IJ/177, 177)+2 ! I=0 for IJ=0 and IJ=31329.
        J=MOD(IJ,     177)+2
        K=MOD(KL/169, 178)+1 ! K=0 for KL=0 and KL=30082.
        M=MOD(KL,     169)
        IF(IUO.GT.5) WRITE(IUO,*) 'RANMAR INITIALIZED.'
!
        DO II=1,97
          S=ZERO
          T=HALF
          DO JJ=1,24
            N=MOD(MOD(I*J,179)*K, 179)
            I=J
            J=K
            K=N
            M=MOD(53*M+1, 169)
            IF(MOD(M*N,64).GE.32) S=S+T
            T=HALF*T
          END DO
          U(II)=S
        END DO
!
        C =  (362436*ONE)/(16777216*ONE) ! 2**24=16777216
        CD= (7654321*ONE)/(16777216*ONE)
        CM=(16777213*ONE)/(16777216*ONE)
!
        I=97
        J=33
!
      ELSE
!
        IF(IUO.GT.5) WRITE(IUO,*) 'MARSAGLIA CONTINUATION.'
        OPEN(UNIT=IUD,FILE=cfile,STATUS='UNKNOWN',FORM='UNFORMATTED')
        REWIND IUD
        READ(IUD) U,C,CD,CM,I,J
        CLOSE(IUD)
      END IF
!
      RETURN
      END
