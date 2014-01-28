! ----------------------------------------------------------------------------
! Numerical diagonalization of 3x3 matrcies
! Copyright (C) 2006  Joachim Kopp
! ----------------------------------------------------------------------------
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.

! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
! ----------------------------------------------------------------------------


! These subroutines were copied verbatim from the source code provided at:
! http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
! by Joachim Kopp and described in:

!      Joachim Kopp
!      Efficient numerical diagonalization of hermitian 3x3 matrices
!      Int. J. Mod. Phys. C 19 (2008) 523-548
!      arXiv.org: physics/0610206


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
    DOUBLE PRECISION :: A(3,3)
    DOUBLE PRECISION :: Q(3,3)
    DOUBLE PRECISION :: W(3)

!     .. Parameters ..
    DOUBLE PRECISION :: EPS
    PARAMETER        ( EPS = 2.2204460492503131D-16 )

!     .. Local Variables ..
    DOUBLE PRECISION :: NORM, N1, N2, N1TMP, N2TMP
    DOUBLE PRECISION :: THRESH, ERROR, WMAX, F, T
    INTEGER ::          I, J

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
    IF (N1 <= THRESH) THEN
        Q(1, 1) = 1.0D0
        Q(2, 1) = 0.0D0
        Q(3, 1) = 0.0D0
    !     If the second column is zero, then (0, 1, 0) is an eigenvector
    ELSE IF (N2 <= THRESH) THEN
        Q(1, 1) = 0.0D0
        Q(2, 1) = 1.0D0
        Q(3, 1) = 0.0D0
    !     If angle between A(*,1) and A(*,2) is too small, don't use
    !     cross product, but calculate v ~ (1, -A0/A1, 0)
    ELSE IF (NORM < (64.0D0 * EPS)**2 * ERROR) THEN
        T = ABS(A(1, 2))
        F = -A(1, 1) / A(1, 2)
        IF (ABS(A(2, 2)) > T) THEN
            T = ABS(A(2, 2))
            F = -A(1, 2) / A(2, 2)
        END IF
        IF (ABS(A(2, 3)) > T) THEN
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
        20 END DO
    END IF
     
!     Prepare calculation of second eigenvector
    T = W(1) - W(2)

!     Is this eigenvalue degenerate?
    IF (ABS(T) > 8.0D0 * EPS * WMAX) THEN
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

        IF (N1 <= THRESH) THEN
            Q(1, 2) = 1.0D0
            Q(2, 2) = 0.0D0
            Q(3, 2) = 0.0D0
        ELSE IF (N2 <= THRESH) THEN
            Q(1, 2) = 0.0D0
            Q(2, 2) = 1.0D0
            Q(3, 2) = 0.0D0
        ELSE IF (NORM < (64.0D0 * EPS)**2 * ERROR) THEN
            T = ABS(A(1, 2))
            F = -A(1, 1) / A(1, 2)
            IF (ABS(A(2, 2)) > T) THEN
                T = ABS(A(2, 2))
                F = -A(1, 2) / A(2, 2)
            END IF
            IF (ABS(A(2, 3)) > T) THEN
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
            40 END DO
        END IF
    ELSE
    !       For degenerate eigenvalue, calculate second eigenvector according to
    !         v[1] = v[0] x (A - lambda[1]).e[i]
    
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
            IF (N1 > THRESH) THEN
                Q(1, 2) = Q(2, 1) * A(3, I) - Q(3, 1) * A(2, I)
                Q(2, 2) = Q(3, 1) * A(1, I) - Q(1, 1) * A(3, I)
                Q(3, 2) = Q(1, 1) * A(2, I) - Q(2, 1) * A(1, I)
                NORM    = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
                IF (NORM > (256.0D0 * EPS)**2 * N1) THEN
                    NORM = SQRT(1.0D0 / NORM)
                    DO 55 J = 1, 3
                        Q(J, 2) = Q(J, 2) * NORM
                    55 END DO
                    GO TO 60
                END IF
            END IF
        50 END DO
           
    !       This means that any vector orthogonal to v[0] is an EV.
        60 IF (I == 4) THEN
            DO 70 J = 1, 3
            !           Find nonzero element of v[0] and swap it with the next one
                IF (Q(J, 1) /= 0.0D0) THEN
                    NORM = 1.0D0 / SQRT(Q(J, 1)**2 + Q(1 + MOD(J,3), 1)**2)
                    Q(J, 2)              = Q(1 + MOD(J,3), 1) * NORM
                    Q(1 + MOD(J,3), 2)   = -Q(J, 1) * NORM
                    Q(1 + MOD(J+1,3), 2) = 0.0D0
                    GO TO 80
                END IF
            70 END DO
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
    SUBROUTINE DSYEVH3(A, Q, W)
! ----------------------------------------------------------------------------
! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
! matrix A using Cardano's method for the eigenvalues and an analytical
! method based on vector cross products for the eigenvectors. However,
! if conditions are such that a large error in the results is to be
! expected, the routine falls back to using the slower, but more
! accurate QL algorithm. Only the diagonal and upper triangular parts of A need
! to contain meaningful values. Access to A is read-only.
! ----------------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   Q: Storage buffer for eigenvectors
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------
! Dependencies:
!   DSYEVC3(), DSYTRD3(), DSYEVQ3()
! ----------------------------------------------------------------------------
! Version history:
!   v1.2 (12 Mar 2012): Removed unused label to avoid gfortran warning,
!     removed unnecessary use of DREAL which led to gfortran error
!   v1.1: Simplified fallback condition --> speed-up
!   v1.0: First released version
! ----------------------------------------------------------------------------
!     .. Arguments ..
    DOUBLE PRECISION :: A(3,3)
    DOUBLE PRECISION :: Q(3,3)
    DOUBLE PRECISION :: W(3)

!     .. Parameters ..
    DOUBLE PRECISION :: EPS
    PARAMETER        ( EPS = 2.2204460492503131D-16 )

!     .. Local Variables ..
    DOUBLE PRECISION :: NORM
    DOUBLE PRECISION :: ERROR
    DOUBLE PRECISION :: T, U
    INTEGER ::          J

!     .. External Functions ..
    EXTERNAL         DSYEVC3, DSYEVQ3

!     Calculate eigenvalues
    CALL DSYEVC3(A, W)

!     --- The rest of this subroutine can be omitted if only the eigenvalues are desired ---

!     Prepare calculation of eigenvectors
!      N1      = A(1, 1)**2 + A(1, 2)**2 + A(1, 3)**2
!      N2      = A(1, 2)**2 + A(2, 2)**2 + A(2, 3)**2
    T       = MAX(ABS(W(1)), ABS(W(2)), ABS(W(3)))
    U       = MAX(T, T**2)
    ERROR   = 256.0D0 * EPS * U**2
!      ERROR   = 256.0D0 * EPS * (N1 + U) * (N2 + U)
    Q(1, 2) = A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
    Q(2, 2) = A(1, 3) * A(1, 2) - A(2, 3) * A(1, 1)
    Q(3, 2) = A(1, 2)**2

!     Calculate first eigenvector by the formula
!       v[0] = (A - lambda[0]).e1 x (A - lambda[0]).e2
    Q(1, 1) = Q(1, 2) + A(1, 3) * W(1)
    Q(2, 1) = Q(2, 2) + A(2, 3) * W(1)
    Q(3, 1) = (A(1,1) - W(1)) * (A(2,2) - W(1)) - Q(3,2)
    NORM    = Q(1, 1)**2 + Q(2, 1)**2 + Q(3, 1)**2

!     If vectors are nearly linearly dependent, or if there might have
!     been large cancellations in the calculation of A(I,I) - W(1), fall
!     back to QL algorithm
!     Note that this simultaneously ensures that multiple eigenvalues do
!     not cause problems: If W(1) = W(2), then A - W(1) * I has rank 1,
!     i.e. all columns of A - W(1) * I are linearly dependent.
    IF (NORM <= ERROR) THEN
        CALL DSYEVQ3(A, Q, W)
        RETURN
    !     This is the standard branch
    ELSE
        NORM = SQRT(1.0D0 / NORM)
        DO 20, J = 1, 3
            Q(J, 1) = Q(J, 1) * NORM
        20 END DO
    END IF
     
!     Calculate second eigenvector by the formula
!       v[1] = (A - lambda[1]).e1 x (A - lambda[1]).e2
    Q(1, 2) = Q(1, 2) + A(1, 3) * W(2)
    Q(2, 2) = Q(2, 2) + A(2, 3) * W(2)
    Q(3, 2) = (A(1,1) - W(2)) * (A(2,2) - W(2)) - Q(3, 2)
    NORM    = Q(1, 2)**2 + Q(2, 2)**2 + Q(3, 2)**2
    IF (NORM <= ERROR) THEN
        CALL DSYEVQ3(A, Q, W)
        RETURN
    ELSE
        NORM = SQRT(1.0D0 / NORM)
        DO 40, J = 1, 3
            Q(J, 2) = Q(J, 2) * NORM
        40 END DO
    END IF

!     Calculate third eigenvector according to
!       v[2] = v[0] x v[1]
    Q(1, 3) = Q(2, 1) * Q(3, 2) - Q(3, 1) * Q(2, 2)
    Q(2, 3) = Q(3, 1) * Q(1, 2) - Q(1, 1) * Q(3, 2)
    Q(3, 3) = Q(1, 1) * Q(2, 2) - Q(2, 1) * Q(1, 2)

    END SUBROUTINE
! End of subroutine DSYEVH3

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
    DOUBLE PRECISION :: A(3,3)
    DOUBLE PRECISION :: W(3)

!     .. Parameters ..
    DOUBLE PRECISION :: SQRT3
    PARAMETER        ( SQRT3 = 1.73205080756887729352744634151D0 )

!     .. Local Variables ..
    DOUBLE PRECISION :: M, C1, C0
    DOUBLE PRECISION :: DE, DD, EE, FF
    DOUBLE PRECISION :: P, SQRTP, Q, C, S, PHI
      
!     Determine coefficients of characteristic poynomial. We write
!           | A   D   F  |
!      A =  | D*  B   E  |
!           | F*  E*  C  |
    DE    = A(1,2) * A(2,3)
    DD    = A(1,2)**2
    EE    = A(2,3)**2
    FF    = A(1,3)**2
    M     = A(1,1) + A(2,2) + A(3,3)
    C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) &
    - (DD + EE + FF)
    C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) &
    - 2.0D0 * A(1,3)*DE

    P     = M**2 - 3.0D0 * C1
    Q     = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
    SQRTP = SQRT(ABS(P))

    PHI   = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1) &
    + C0 * (Q + (27.0D0/4.0D0)*C0) )
    PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)

    C     = SQRTP * COS(PHI)
    S     = (1.0D0/SQRT3) * SQRTP * SIN(PHI)

    W(2) = (1.0D0/3.0D0) * (M - C)
    W(3) = W(2) + S
    W(1) = W(2) + C
    W(2) = W(2) - S

    END SUBROUTINE
! End of subroutine DSYEVC3

! ----------------------------------------------------------------------------
    SUBROUTINE DSYTRD3(A, Q, D, E)
! ----------------------------------------------------------------------------
! Reduces a symmetric 3x3 matrix to real tridiagonal form by applying
! (unitary) Householder transformations:
!            [ D[1]  E[1]       ]
!    A = Q . [ E[1]  D[2]  E[2] ] . Q^T
!            [       E[2]  D[3] ]
! The function accesses only the diagonal and upper triangular parts of
! A. The access is read-only.
! ---------------------------------------------------------------------------
!     .. Arguments ..
    DOUBLE PRECISION :: A(3,3)
    DOUBLE PRECISION :: Q(3,3)
    DOUBLE PRECISION :: D(3)
    DOUBLE PRECISION :: E(2)

!     .. Parameters ..
    INTEGER ::          N
    PARAMETER        ( N = 3 )

!     .. Local Variables ..
    DOUBLE PRECISION :: U(N), P(N)
    DOUBLE PRECISION :: OMEGA, F
    DOUBLE PRECISION :: K, H, G
    INTEGER ::          I, J

!     Initialize Q to the identitity matrix
!     --- This loop can be omitted if only the eigenvalues are desired ---
    DO 10 I = 1, N
        Q(I,I) = 1.0D0
        DO 11, J = 1, I-1
            Q(I, J) = 0.0D0
            Q(J, I) = 0.0D0
        11 END DO
    10 END DO

!     Bring first row and column to the desired form
    H = A(1,2)**2 + A(1,3)**2
    IF (A(1,2) > 0.0D0) THEN
        G = -SQRT(H)
    ELSE
        G = SQRT(H)
    END IF
    E(1)  = G
    F     = G * A(1,2)
    U(2)  = A(1,2) - G
    U(3)  = A(1,3)

    OMEGA = H - F
    IF (OMEGA > 0.0D0) THEN
        OMEGA = 1.0D0 / OMEGA
        K     = 0.0D0
        DO 20 I = 2, N
            F    = A(2,I)*U(2) + A(I,3)*U(3)
            P(I) = OMEGA * F
            K    = K + U(I) * F
        20 END DO
        K = 0.5D0 * K * OMEGA**2

        DO 30 I = 2, N
            P(I) = P(I) - K * U(I)
        30 END DO

        D(1) = A(1,1)
        D(2) = A(2,2) - 2.0D0 * P(2) * U(2)
        D(3) = A(3,3) - 2.0D0 * P(3) * U(3)

    !       Store inverse Householder transformation in Q
    !       --- This loop can be omitted if only the eigenvalues are desired ---
        DO 40, J = 2, N
            F = OMEGA * U(J)
            DO 41 I = 2, N
                Q(I,J) = Q(I,J) - F * U(I)
            41 END DO
        40 END DO
                    
    !       Calculated updated A(2, 3) and store it in E(2)
        E(2) = A(2, 3) - P(2) * U(3) - U(2) * P(3)
    ELSE
        DO 50 I = 1, N
            D(I) = A(I, I)
        50 END DO
        E(2) = A(2, 3)
    END IF
          
    END SUBROUTINE
! End of subroutine DSYTRD3

! ----------------------------------------------------------------------------
    SUBROUTINE DSYEVQ3(A, Q, W)
! ----------------------------------------------------------------------------
! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
! matrix A using the QL algorithm with implicit shifts, preceded by a
! Householder reduction to real tridiagonal form.
! The function accesses only the diagonal and upper triangular parts of
! A. The access is read-only.
! ----------------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   Q: Storage buffer for eigenvectors
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------
! Dependencies:
!   DSYTRD3()
! ----------------------------------------------------------------------------
!     .. Arguments ..
    DOUBLE PRECISION :: A(3,3)
    DOUBLE PRECISION :: Q(3,3)
    DOUBLE PRECISION :: W(3)

!     .. Parameters ..
    INTEGER ::          N
    PARAMETER        ( N = 3 )

!     .. Local Variables ..
    DOUBLE PRECISION :: E(3)
    DOUBLE PRECISION :: G, R, P, F, B, S, C, T
    INTEGER ::          NITER
    INTEGER ::          L, M, I, J, K

!     .. External Functions ..
    EXTERNAL         DSYTRD3
          
!     Transform A to real tridiagonal form by the Householder method
    CALL DSYTRD3(A, Q, W, E)

!     Calculate eigensystem of the remaining real symmetric tridiagonal
!     matrix with the QL method

!     Loop over all off-diagonal elements
    DO 10 L = 1, N-1
        NITER = 0

    !       Iteration loop
        DO 11 I = 1, 50
        !         Check for convergence and exit iteration loop if off-diagonal
        !         element E(L) is zero
            DO 20 M = L, N-1
                G = ABS(W(M)) + ABS(W(M+1))
                IF (ABS(E(M)) + G == G) THEN
                    GO TO 30
                END IF
            20 END DO
            30 IF (M == L) THEN
                GO TO 10
            END IF

            NITER = NITER + 1
            IF (NITER >= 30) THEN
                PRINT *, 'DSYEVQ3: No convergence.'
                RETURN
            END IF

        !         Calculate G = D(M) - K
            G = (W(L+1) - W(L)) / (2.0D0 * E(L))
            R = SQRT(1.0D0 + G**2)
            IF (G >= 0.0D0) THEN
                G = W(M) - W(L) + E(L)/(G + R)
            ELSE
                G = W(M) - W(L) + E(L)/(G - R)
            END IF

            S = 1.0D0
            C = 1.0D0
            P = 0.0D0
            DO 40 J = M - 1, L, -1
                F = S * E(J)
                B = C * E(J)
                IF (ABS(F) > ABS(G)) THEN
                    C      = G / F
                    R      = SQRT(1.0D0 + C**2)
                    E(J+1) = F * R
                    S      = 1.0D0 / R
                    C      = C * S
                ELSE
                    S      = F / G
                    R      = SQRT(1.0D0 + S**2)
                    E(J+1) = G * R
                    C      = 1.0D0 / R
                    S      = S * C
                END IF

                G      = W(J+1) - P
                R      = (W(J) - G) * S + 2.0D0 * C * B
                P      = S * R
                W(J+1) = G + P
                G      = C * R - B

            !           Form eigenvectors
            !           --- This loop can be omitted if only the eigenvalues are desired ---
                DO 50 K = 1, N
                    T         = Q(K, J+1)
                    Q(K, J+1) = S * Q(K, J) + C * T
                    Q(K, J)   = C * Q(K, J) - S * T
                50 END DO
            40 END DO
            W(L) = W(L) - P
            E(L) = G
            E(M) = 0.0D0
        11 END DO
    10 END DO
      
    END SUBROUTINE
! End of subroutine DSYEVQ3
