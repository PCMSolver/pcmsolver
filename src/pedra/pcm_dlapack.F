! ----------------------------------------------------------------------------
! Numerical diagonalization of 3x3 matrcies
! Copyright (C) 2006  Joachim Kopp
! ----------------------------------------------------------------------------
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
! ----------------------------------------------------------------------------


! ----------------------------------------------------------------------------
      SUBROUTINE DSYEVV3(A, Q, W)
! ----------------------------------------------------------------------------
! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
! matrix A using Cardano's method for the eigenvalues and an analytical
! method based on vector cross products for the eigenvectors.
! Only the diagonal and upper triangular parts of A need to contain
! meaningful values. However, all of A may be used as temporary storage
! and may hence be destroyed.
! ----------------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   Q: Storage buffer for eigenvectors
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------
! Dependencies:
!   DSYEVC3()
! ----------------------------------------------------------------------------
! Version history:
!   v1.1 (12 Mar 2012): Removed access to lower triangualr part of A
!     (according to the documentation, only the upper triangular part needs
!     to be filled)
!   v1.0: First released version
! ----------------------------------------------------------------------------
!     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

!     .. Parameters ..
      DOUBLE PRECISION EPS
      PARAMETER        ( EPS = 2.2204460492503131D-16 )

!     .. Local Variables ..
      DOUBLE PRECISION NORM, N1, N2, N1TMP, N2TMP
      DOUBLE PRECISION THRESH, ERROR, WMAX, F, T
      INTEGER          I, J

!     .. External Functions ..
      EXTERNAL         DSYEVC3

!     Calculate eigenvalues
      CALL DSYEVC3(A, W)

!     --- The rest of this subroutine can be omitted if only the eigenvalues are desired ---

      WMAX   = MAX(ABS(W(1)), ABS(W(2)), ABS(W(3)))
      THRESH = (8.0D0 * EPS * WMAX)**2

!     Prepare calculation of eigenvectors
      N1TMP   = A(1, 2)**2 + A(1, 3)**2
      N2TMP   = A(1, 2)**2 + A(2, 3)**2
      Q(1, 1) = A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
      Q(1, 2) = Q(1, 1)
      Q(2, 1) = A(1, 3) * A(1, 2) - A(2, 3) * A(1, 1)
      Q(2, 2) = Q(2, 1)
      Q(3, 2) = A(1, 2)**2

!     Calculate first eigenvector by the formula
!       v[0] = (A - lambda[0]).e1 x (A - lambda[0]).e2
      A(1, 1) = A(1, 1) - W(1)
      A(2, 2) = A(2, 2) - W(1)
      Q(1, 1) = Q(1, 2) + A(1, 3) * W(1)
      Q(2, 1) = Q(2, 2) + A(2, 3) * W(1)
      Q(3, 1) = A(1, 1) * A(2, 2) - Q(3, 2)
      NORM    = Q(1, 1)**2 + Q(2, 1)**2 + Q(3, 1)**2
      N1      = N1TMP + A(1, 1)**2
      N2      = N2TMP + A(2, 2)**2
      ERROR   = N1 * N2

!     If the first column is zero, then (1, 0, 0) is an eigenvector
      IF (N1 .LE. THRESH) THEN
        Q(1, 1) = 1.0D0
        Q(2, 1) = 0.0D0
        Q(3, 1) = 0.0D0
!     If the second column is zero, then (0, 1, 0) is an eigenvector
      ELSE IF (N2 .LE. THRESH) THEN
        Q(1, 1) = 0.0D0
        Q(2, 1) = 1.0D0
        Q(3, 1) = 0.0D0
!     If angle between A(*,1) and A(*,2) is too small, don't use
!     cross product, but calculate v ~ (1, -A0/A1, 0)
      ELSE IF (NORM .LT. (64.0D0 * EPS)**2 * ERROR) THEN
        T = ABS(A(1, 2))
        F = -A(1, 1) / A(1, 2)
        IF (ABS(A(2, 2)) .GT. T) THEN
          T = ABS(A(2, 2))
          F = -A(1, 2) / A(2, 2)
        END IF
        IF (ABS(A(2, 3)) .GT. T) THEN
          F = -A(1, 3) / A(2, 3)
        END IF
        NORM    = 1.0D0 / SQRT(1.0D0 + F**2)
        Q(1, 1) = NORM
        Q(2, 1) = F * NORM
        Q(3, 1) = 0.0D0
!     This is the standard branch
      ELSE
        NORM = SQRT(1.0D0 / NORM)
        DO 20, J = 1, 3
          Q(J, 1) = Q(J, 1) * NORM
   20   CONTINUE
      END IF
 
!     Prepare calculation of second eigenvector     
      T = W(1) - W(2)

!     Is this eigenvalue degenerate?
      IF (ABS(T) .GT. 8.0D0 * EPS * WMAX) THEN
!       For non-degenerate eigenvalue, calculate second eigenvector by
!       the formula
!         v[1] = (A - lambda[1]).e1 x (A - lambda[1]).e2
        A(1, 1) = A(1, 1) + T
        A(2, 2) = A(2, 2) + T
        Q(1, 2) = Q(1, 2) + A(1, 3) * W(2)
        Q(2, 2) = Q(2, 2) + A(2, 3) * W(2)
        Q(3, 2) = A(1, 1) * A(2, 2) - Q(3, 2)
        NORM    = Q(1, 2)**2 + Q(2, 2)**2 + Q(3, 2)**2
        N1      = N1TMP + A(1, 1)**2
        N2      = N2TMP + A(2, 2)**2
        ERROR   = N1 * N2

        IF (N1 .LE. THRESH) THEN
          Q(1, 2) = 1.0D0
          Q(2, 2) = 0.0D0
          Q(3, 2) = 0.0D0
        ELSE IF (N2 .LE. THRESH) THEN
          Q(1, 2) = 0.0D0
          Q(2, 2) = 1.0D0
          Q(3, 2) = 0.0D0
        ELSE IF (NORM .LT. (64.0D0 * EPS)**2 * ERROR) THEN
          T = ABS(A(1, 2))
          F = -A(1, 1) / A(1, 2)
          IF (ABS(A(2, 2)) .GT. T) THEN
            T = ABS(A(2, 2))
            F = -A(1, 2) / A(2, 2)
          END IF
          IF (ABS(A(2, 3)) .GT. T) THEN
            F = -A(1, 3) / A(2, 3)
          END IF
          NORM    = 1.0D0 / SQRT(1.0D0 + F**2)
          Q(1, 2) = NORM
          Q(2, 2) = F * NORM
          Q(3, 2) = 0.0D0
        ELSE
          NORM = SQRT(1.0D0 / NORM)
          DO 40, J = 1, 3 
            Q(J, 2) = Q(J, 2) * NORM
   40     CONTINUE
        END IF
      ELSE
!       For degenerate eigenvalue, calculate second eigenvector according to
!         v[1] = v[0] x (A - lambda[1]).e[i]
!   
!       This would really get to complicated if we could not assume all of A to
!       contain meaningful values.
        A(2, 1) = A(1, 2)
        A(3, 1) = A(1, 3)
        A(3, 2) = A(2, 3)
        A(1, 1) = A(1, 1) + W(1)
        A(2, 2) = A(2, 2) + W(1)
        DO 50 I = 1, 3
          A(I, I) = A(I, I) - W(2)
          N1      = A(1, I)**2 + A(2, I)**2 + A(3, I)**2
          IF (N1 .GT. THRESH) THEN
            Q(1, 2) = Q(2, 1) * A(3, I) - Q(3, 1) * A(2, I)
            Q(2, 2) = Q(3, 1) * A(1, I) - Q(1, 1) * A(3, I)
            Q(3, 2) = Q(1, 1) * A(2, I) - Q(2, 1) * A(1, I)
            NORM    = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
            IF (NORM .GT. (256.0D0 * EPS)**2 * N1) THEN
              NORM = SQRT(1.0D0 / NORM)
              DO 55 J = 1, 3
                Q(J, 2) = Q(J, 2) * NORM
   55         CONTINUE
              GO TO 60
            END IF
          END IF
   50   CONTINUE
   
!       This means that any vector orthogonal to v[0] is an EV.
   60   IF (I .EQ. 4) THEN
          DO 70 J = 1, 3
!           Find nonzero element of v[0] and swap it with the next one
            IF (Q(J, 1).NE.0.0D0) THEN
              NORM = 1.0D0 / SQRT(Q(J, 1)**2 + Q(1 + MOD(J,3), 1)**2)
              Q(J, 2)              = Q(1 + MOD(J,3), 1) * NORM
              Q(1 + MOD(J,3), 2)   = -Q(J, 1) * NORM
              Q(1 + MOD(J+1,3), 2) = 0.0D0
              GO TO 80
            END IF
   70     CONTINUE
        END IF
      END IF

!     Calculate third eigenvector according to
!       v[2] = v[0] x v[1]
   80 Q(1, 3) = Q(2, 1) * Q(3, 2) - Q(3, 1) * Q(2, 2)
      Q(2, 3) = Q(3, 1) * Q(1, 2) - Q(1, 1) * Q(3, 2)
      Q(3, 3) = Q(1, 1) * Q(2, 2) - Q(2, 1) * Q(1, 2)

      END SUBROUTINE
! End of subroutine DSYEVV3

! ----------------------------------------------------------------------------
      SUBROUTINE DSYEVC3(A, W)
! ----------------------------------------------------------------------------
! Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
! analytical algorithm.
! Only the diagonal and upper triangular parts of A are accessed. The access
! is read-only.
! ----------------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------
!     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION W(3)

!     .. Parameters ..
      DOUBLE PRECISION SQRT3
      PARAMETER        ( SQRT3 = 1.73205080756887729352744634151D0 )

!     .. Local Variables ..
      DOUBLE PRECISION M, C1, C0
      DOUBLE PRECISION DE, DD, EE, FF
      DOUBLE PRECISION P, SQRTP, Q, C, S, PHI
  
!     Determine coefficients of characteristic poynomial. We write
!           | A   D   F  |
!      A =  | D*  B   E  |
!           | F*  E*  C  |
      DE    = A(1,2) * A(2,3)
      DD    = A(1,2)**2
      EE    = A(2,3)**2
      FF    = A(1,3)**2
      M     = A(1,1) + A(2,2) + A(3,3)
      C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) )
     $         - (DD + EE + FF)
      C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3)
     $         - 2.0D0 * A(1,3)*DE

      P     = M**2 - 3.0D0 * C1
      Q     = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
      SQRTP = SQRT(ABS(P))

      PHI   = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1)
     $          + C0 * (Q + (27.0D0/4.0D0)*C0) )
      PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)

      C     = SQRTP * COS(PHI)
      S     = (1.0D0/SQRT3) * SQRTP * SIN(PHI)

      W(2) = (1.0D0/3.0D0) * (M - C)
      W(3) = W(2) + S
      W(1) = W(2) + C
      W(2) = W(2) - S

      END SUBROUTINE
! End of subroutine DSYEVC3
