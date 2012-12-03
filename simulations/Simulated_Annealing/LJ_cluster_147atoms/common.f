
C *****************************************************************************
C **                                                                         **
C **  Common area definition for the NVT-MC-program.                         **
C **                                                                         **
C **  NMAX:     maximum allowed number of atoms in the system                **
C **  N:        number of atoms in the system (read in from conf file)       **
C **  RX:       x coordinates of the atoms                                   **
C **  RY:       y coordinates of the atoms                                   **
C **  RZ:       z coordinates of the atoms                                   **
C **                                                                         **
C *****************************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 2000 )
        REAL        RX(NMAX), RY(NMAX), RZ(NMAX)
        INTEGER N
        COMMON /COORDS/ RX, RY, RZ, N
