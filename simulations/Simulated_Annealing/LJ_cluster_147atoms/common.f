! *****************************************************************************
! **                                                                         **
! **  Common area definition for the NVT-MC-program.                         **
! **                                                                         **
! **  NMAX:     maximum allowed number of atoms in the system                **
! **  N:        number of atoms in the system (read in from conf file)       **
! **  RX:       x coordinates of the atoms                                   **
! **  RY:       y coordinates of the atoms                                   **
! **  RZ:       z coordinates of the atoms                                   **
! **                                                                         **
! *****************************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 125 )
        REAL*8        RX(NMAX), RY(NMAX), RZ(NMAX)
        INTEGER N
        COMMON /COORDS/ RX, RY, RZ, N
