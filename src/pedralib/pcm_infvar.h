C
C file: infvar.h (originally from SIRIUS)
C       contains INFormation on number and symmetry of MCSCF VARiables
C
      INTEGER MAXWOP, JWOP,
     &        NCONF, NWOPT, NVAR, JWOPSY, NWOP, NWOPH, NVARH, NCDETS
C     MAXWOP = maximum number of orbital rotations (dimension of JWOP)
      PARAMETER ( MAXWOP = 200000 )
      COMMON /INFVAR/ JWOP(2,MAXWOP),
     &                NCONF,NWOPT,NVAR,JWOPSY,NWOP(8),NWOPH,NVARH,NCDETS
