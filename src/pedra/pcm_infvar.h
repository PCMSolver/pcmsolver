!
! file: infvar.h (originally from SIRIUS)
!       contains INFormation on number and symmetry of MCSCF VARiables
!
      INTEGER MAXWOP, JWOP,                                              &
             NCONF, NWOPT, NVAR, JWOPSY, NWOP, NWOPH, NVARH, NCDETS
!     MAXWOP = maximum number of orbital rotations (dimension of JWOP)
      PARAMETER ( MAXWOP = 200000 )
      COMMON /PCM_INFVAR/ JWOP(2,MAXWOP),                                &
                     NCONF,NWOPT,NVAR,JWOPSY,NWOP(8),NWOPH,NVARH,NCDETS
