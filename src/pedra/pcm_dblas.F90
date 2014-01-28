!  /* Deck blas */
    FUNCTION DASUM_(N,DX,INCX)

!     RETURNS SUM OF MAGNITUDES OF DOUBLE PRECISION DX.
!     DASUM_ = SUM FROM 0 TO N-1 OF DABS(DX(1+I*INCX))

    DOUBLE PRECISION :: DASUM_, DASUM
    DOUBLE PRECISION :: DX(1)
    DASUM_ = 0.D0
    DASUM  = 0.0D0
    IF(N <= 0)RETURN
    IF(INCX == 1)GOTO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

    NS = N*INCX
    DO 10 I=1,NS,INCX
        DASUM_ = DASUM + DABS(DX(I))
    10 END DO
    RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.


!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.

    20 M = MOD(N,6)
    IF( M == 0 ) GO TO 40
    DO 30 I = 1,M
        DASUM_ = DASUM + DABS(DX(I))
    30 END DO
    IF( N < 6 ) RETURN
    40 MP1 = M + 1
    DO 50 I = MP1,N,6
        DASUM_ = DASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2)) &
        + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
    50 END DO
    RETURN
    END FUNCTION 

    SUBROUTINE DAXPY_(N,DA,DX,INCX,DY,INCY)

!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION :: DX(1),DY(1),DA
    INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N
    DATA ZERO/0.0D0/

    IF (N <= 0)       RETURN
    IF (DA == ZERO) RETURN
    IF (INCX == 1 .AND. INCY == 1) GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

    IX = 1
    IY = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    IF(INCY < 0)IY = (-N+1)*INCY + 1
    DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
    10 END DO
    RETURN

    20 DO 30 I = 1, N
        DY(I) = DY(I) + DA*DX(I)
    30 END DO
    RETURN
    END SUBROUTINE DAXPY_

    FUNCTION DDOT_(N,DX,INCX,DY,INCY)

!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION :: DDOT_, DX(1),DY(1),DTEMP
    INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N
    DATA ZERO/0.0D0/

    DDOT_ = ZERO
    DTEMP = ZERO
    IF(N <= 0)RETURN
    IF(INCX == 1 .AND. INCY == 1)GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

    IX = 1
    IY = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    IF(INCY < 0)IY = (-N+1)*INCY + 1
    DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
    10 END DO
    DDOT_ = DTEMP
    RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

    20 M = MOD(N,5)
    IF( M == 0 ) GO TO 40
    DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
    30 END DO
    IF( N < 5 ) GO TO 60
    40 MP1 = M + 1
    DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) + &
        DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
    50 END DO
    60 DDOT_ = DTEMP
    RETURN
    END FUNCTION DDOT_

    FUNCTION DNRM2_(N,DX,INCX)

!     Forms the two-norm of a vector.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
! 30-Apr-1984 -- hjaaj -- based on DDOT_ from LINPACK
!     DNRM2_(N,DX,INCX) = DSQRT( DDOT_(N,DX,INCX,DX,INCX) )
!     This version does not use extended precision for intermediates
!     as the original LINPACK version does.
!     DNRM2_: JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION :: DNRM2_,DX(1),DTEMP
    INTEGER :: I,INCX,IX,M,MP1,N
    DATA ZERO/0.0D0/

    DNRM2_ = ZERO
    IF(N <= 0)RETURN
    DTEMP = ZERO
    IF(INCX == 1)GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

    IX = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DX(IX)
        IX = IX + INCX
    10 END DO
    DNRM2_ = DSQRT(DTEMP)
    RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

    20 M = MOD(N,5)
    IF( M == 0 ) GO TO 40
    DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DX(I)
    30 END DO
    IF( N < 5 ) GO TO 60
    40 MP1 = M + 1
    DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DX(I) + DX(I + 1)*DX(I + 1) + &
        DX(I + 2)*DX(I + 2) + DX(I + 3)*DX(I + 3) + DX(I + 4)*DX(I + 4)
    50 END DO
    60 DNRM2_ = DSQRT(DTEMP)
    RETURN
    END FUNCTION DNRM2_

    SUBROUTINE DSCAL_(N,DA,DX,INCX)

!     CONSTANT TIMES A VECTOR
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
! 30-Apr-1984 -- hjaaj -- based on DAXPY_ from LINPACK
!     DAXPY_: JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION :: DA,DX(*)
    INTEGER :: I,INCX,M,MP1,N

    IF(N <= 0)RETURN
    IF(INCX == 1)GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

    IX = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
    10 END DO
    RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

    20 CONTINUE
    M = MOD(N,4)
    IF( M == 0 ) GO TO 40
    DO 30 I = 1,M
        DX(I) = DA*DX(I)
    30 END DO
    IF( N < 4 ) RETURN
    40 CONTINUE
    MP1 = M + 1
    DO 50 I = MP1,N,4
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
    50 END DO
    RETURN
    END SUBROUTINE DSCAL_

    SUBROUTINE  DSWAP_ (N,DX,INCX,DY,INCY)

!     INTER_CHANGES TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION :: DX(1),DY(1),DTEMP
    INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N

    IF(N <= 0)RETURN
    IF(INCX == 1 .AND. INCY == 1)GO TO 20

!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1

    IX = 1
    IY = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    IF(INCY < 0)IY = (-N+1)*INCY + 1
    DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
    10 END DO
    RETURN

!       CODE FOR BOTH INCREMENTS EQUAL TO 1


!       CLEAN-UP LOOP

    20 M = MOD(N,3)
    IF( M == 0 ) GO TO 40
    DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
    30 END DO
    IF( N < 3 ) RETURN
    40 MP1 = M + 1
    DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
    50 END DO
    RETURN
    END SUBROUTINE 

    INTEGER FUNCTION IDAMAX_(N,DX,INCX)

!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
!     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION ::  DX(1),DMAX
    INTEGER :: I,INCX,IX,N

    IDAMAX_ = 0
    IF( N < 1 ) RETURN
    IDAMAX_ = 1
    IF(N == 1)RETURN
    IF(INCX == 1)GO TO 20

!        CODE FOR INCREMENT NOT EQUAL TO 1

    IX = 1
    DMAX = DABS(DX(1))
    IX = IX + INCX
    DO 10 I = 2,N
        IF(DABS(DX(IX)) <= DMAX) GO TO 5
        IDAMAX_ = I
        DMAX = DABS(DX(IX))
        5 IX = IX + INCX
    10 END DO
    RETURN

!        CODE FOR INCREMENT EQUAL TO 1

    20 DMAX = DABS(DX(1))
    DO 30 I = 2,N
        IF(DABS(DX(I)) <= DMAX) GO TO 30
        IDAMAX_ = I
        DMAX = DABS(DX(I))
    30 END DO
    RETURN
    END FUNCTION IDAMAX_

    SUBROUTINE DCOPY_(N,DX,INCX,DY,INCY)

!     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY.
!     FOR I = 0 TO N-1, COPY DX(LX+I*INCX) TO DY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.

    DOUBLE PRECISION :: DX(1),DY(1)
    IF(N <= 0)RETURN
    IF (INCX == INCY) THEN
        IF (INCX == 1) GOTO 20
        IF (INCX > 1) GOTO 60
    END IF
    5 CONTINUE

!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.

    IX = 1
    IY = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    IF(INCY < 0)IY = (-N+1)*INCY + 1
    DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
    10 END DO
    RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.

    20 M = MOD(N,7)
    IF( M == 0 ) GO TO 40
    DO 30 I = 1,M
        DY(I) = DX(I)
    30 END DO
    IF( N < 7 ) RETURN
    40 MP1 = M + 1
    DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
    50 END DO
    RETURN

!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.

    60 CONTINUE
    NS=N*INCX
    DO 70 I=1,NS,INCX
        DY(I) = DX(I)
    70 END DO
    RETURN
    END SUBROUTINE DCOPY_
