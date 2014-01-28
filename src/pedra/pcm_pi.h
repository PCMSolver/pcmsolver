! File: pi.h
!
      REAL*8     PI, PI2, SQRTPI, R2PI52
      PARAMETER (PI     = 3.14159265358979323846D00, PI2 = PI*PI,    &
                SQRTPI = 1.77245385090551602730D00,                  &
                R2PI52 = 5.91496717279561287782D00)
!     R2PI52 = sqrt(2 * sqrt(PI^5) ) -- used in calc. of 2-el. integrals
! -- end of pi.h --
