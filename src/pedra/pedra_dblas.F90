    module pedra_dblas

    implicit none

    public dasum
    public daxpy
    public dscal
    public dswap
    public dcopy
    public dzero
    public ddot
    public dnorm2
    public idamax
    public vector_product

    contains

    function dasum(n,dx,incx)

!     RETURNS SUM OF MAGNITUDES OF DOUBLE PRECISION DX.
!     DASUM_ = SUM FROM 0 TO N-1 OF DABS(DX(1+I*INCX))

    real(8)             :: dasum
    integer, intent(in) :: n
    real(8), intent(in) :: dx(n)
    integer, intent(in) :: incx

    real(8) :: temp_sum = 0.0d0
    integer :: m, mp1, i, ns

    dasum = 0.0d0
    
    if(n <= 0)return
    
    if (incx == 1) then
            m = mod(n, 6)
            if (m == 0) then
                    mp1 = m + 1
                    do i = mp1, n, 6
                        temp_sum = dasum + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2)) &
                      + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))
                    end do
            else
                    do i = 1, m
                        temp_sum = temp_sum + dabs(dx(i))
                    end do
            end if
    else
            ns = n * incx
            do i = 1, ns, incx
                temp_sum = temp_sum + dabs(dx(i))
            end do
    end if

    end function dasum

    subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     dy = da*dx + dy
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
    
    integer,    intent(in) :: n
    real(8),    intent(in) :: dx(n)
    real(8),    intent(in) :: da
    real(8), intent(inout) :: dy(n)
    integer,    intent(in) :: incx, incy

    integer :: i, ix, iy
    real(8), parameter :: zero = 0.0d0

    if (n <= 0)       return
    if (da == zero) return
    if (incx == 1 .and. incy == 1) then
            do i = 1, n
                dy(i) = dy(i) + da*dx(i)
            end do
    else
            ix = 1
            iy = 1
            if (incx < 0) then
                    ix = (-n+1)*incx + 1
            end if
            if (incy < 0) then
                    iy = (-n+1)*incy + 1
            end if
            do i = 1, n
                dy(iy) = dy(iy) + da*dx(ix)
                ix = ix + incx
                iy = iy + incy
            end do
    end if

    end subroutine daxpy

    function ddot(n,dx,incx,dy,incy)

!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

    real(8) :: ddot
    integer, intent(in) :: n, incx, incy
    real(8), intent(in) :: dx(n), dy(n)

    real(8) :: dtemp
    integer :: I, ix, iy, m, mp1
    real(8), parameter :: zero = 0.0d0

    DDOT = ZERO
    DTEMP = ZERO
    IF(N <= 0)RETURN
    IF(INCX == 1 .AND. INCY == 1)go to 20

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
    10 end do
    DDOT = DTEMP
    RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

    20 M = MOD(N,5)
    IF( M == 0 ) go to 40
    DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
    30 end do
    IF( N < 5 ) go to 60
    40 MP1 = M + 1
    DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) + &
        DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
    50 end do
    60 DDOT = DTEMP

    END FUNCTION DDOT
    
    subroutine dzero(dx, length)

!...................................................................
! Last revision 5-May-1984 by Hans Jorgen Aa. Jensen
!
!   Subroutine DZERO sets a real array of length *LENGTH*
!   to zero.
!...................................................................
    integer, intent(in)    :: length
    real(8), intent(inout) :: dx(length)

    integer :: i

    if (length <= 0) return

    do i = 1, length
        dx(i) = 0.0d0
    end do

    end subroutine dzero
    
    function dnorm2(n, dx, incx)

!     Forms the two-norm of a vector.
! 19-Sep-1988 -- hjaaj -- based on DNRM2_ from LINPACK
!     This version does not use extended precision for intermediates
!     as the LINPACK version does.C 1) DNORM2_ (emulate ESSL DNORM2: do not use \
! xtended precision for intermediates)

!     Equivalent to DNORM2_ in IBM's ESSL library.

!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     DNRM2_: JACK DONGARRA, LINPACK, 3/11/78.

    real(8)             :: dnorm2
    integer, intent(in) :: n
    real(8), intent(in) :: dx(n)
    integer, intent(in) :: incx
    real(8), parameter  :: zero = 0.0d0 

    real(8) :: dtemp
    integer :: m, mp1, i, ix

    if (n <= 0) then
        dnorm2 = zero
        return
    end if
    
    dtemp  = zero
    if (incx == 1) then
            m = mod(n, 5)
            if (m == 0) then
                    mp1 = m + 1
                    do i = mp1, n, 5
                        dtemp = dtemp + dx(i)*dx(i) + dx(i + 1)*dx(i + 1) + &
                        dx(i + 2)*dx(i + 2) + dx(i + 3)*dx(i + 3) + dx(i + 4)*dx(i + 4)
                    end do
            else
                    do i = 1, m
                        dtemp = dtemp + dx(i)*dx(i)
                    end do
            end if
    else
            ix = 1
            if (incx < 0) then
                    ix = (-n + 1)*incx + 1
            end if
            do i = 1, n
                dtemp = dtemp + dx(ix)*dx(ix)
                ix = ix + incx
            end do
    end if

    dnorm2 = sqrt(dtemp)

    end function dnorm2


    subroutine dscal(n,da,dx,incx)

!     CONSTANT TIMES A VECTOR
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
! 30-Apr-1984 -- hjaaj -- based on DAXPY_ from LINPACK
!     DAXPY_: JACK DONGARRA, LINPACK, 3/11/78.

    integer,    intent(in) :: n
    real(8),    intent(in) :: da
    real(8), intent(inout) :: dx(n)
    integer,    intent(in) :: incx
    integer :: i, ix, m, mp1

    IF(N <= 0)RETURN
    IF(INCX == 1)go to 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

    IX = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
    10 end do
    RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

    20 CONTINUE
    M = MOD(N,4)
    IF( M == 0 ) go to 40
    DO 30 I = 1,M
        DX(I) = DA*DX(I)
    30 end do
    IF( N < 4 ) RETURN
    40 CONTINUE
    MP1 = M + 1
    DO 50 I = MP1,N,4
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
    50 end do
    RETURN
    END SUBROUTINE DSCAL


    INTEGER FUNCTION IDAMAX(N,DX,INCX)

!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
!     JACK DONGARRA, LINPACK, 3/11/78.

    integer, intent(in) :: n, incx
    real(8), intent(in) :: dx(n)
    real(8) :: dmax
    INTEGER :: I, IX

    IDAMAX = 0
    IF( N < 1 ) RETURN
    IDAMAX = 1
    IF(N == 1)RETURN
    IF(INCX == 1)go to 20

!        CODE FOR INCREMENT NOT EQUAL TO 1

    IX = 1
    DMAX = DABS(DX(1))
    IX = IX + INCX
    DO 10 I = 2,N
        IF(DABS(DX(IX)) <= DMAX) go to 5
        IDAMAX = I
        DMAX = DABS(DX(IX))
        5 IX = IX + INCX
    10 end do
    RETURN

!        CODE FOR INCREMENT EQUAL TO 1

    20 DMAX = DABS(DX(1))
    DO 30 I = 2,N
        IF(DABS(DX(I)) <= DMAX) go to 30
        IDAMAX = I
        DMAX = DABS(DX(I))
    30 end do
    RETURN
    END FUNCTION IDAMAX

    SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!
!  -Reference BLAS level1 routine (version 3.4.0) --
!  -Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2011
!
!   .. Scalar Arguments ..
    INTEGER INCX,INCY,N
!   ..
!   .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
!   ..
!
!  ===================================================================
!
!   .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY,M,MP1
!   ..
!   .. Intrinsic Functions ..
    INTRINSIC MOD
!   ..
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!     code for both increments equal to 1
!
!
!     clean-up loop
!
       M = MOD(N,3)
       IF (M.NE.0) THEN
          DO I = 1,M
             DTEMP = DX(I)
             DX(I) = DY(I)
             DY(I) = DTEMP
          END DO
          IF (N.LT.3) RETURN
       END IF
       MP1 = M + 1
       DO I = MP1,N,3
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
          DTEMP = DX(I+1)
          DX(I+1) = DY(I+1)
          DY(I+1) = DTEMP
          DTEMP = DX(I+2)
          DX(I+2) = DY(I+2)
          DY(I+2) = DTEMP
       END DO
    ELSE
!
!     code for unequal increments or equal increments not equal
!       to 1
!
       IX = 1
       IY = 1
       IF (INCX.LT.0) IX = (-N+1)*INCX + 1
       IF (INCY.LT.0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          DTEMP = DX(IX)
          DX(IX) = DY(IY)
          DY(IY) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
       END DO
    END IF
    RETURN
    END subroutine

    SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)

!     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY.
!     FOR I = 0 TO N-1, COPY DX(LX+I*INCX) TO DY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.

    integer, intent(in) :: n, incx, incy
    real(8)             :: dx(n), dy(n)

    integer :: ix, iy, i, m, mp1, ns

    IF(N <= 0)RETURN
    IF (INCX == INCY) THEN
        IF (INCX == 1) go to 20
        IF (INCX > 1) go to 60
    end if
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.

    IX = 1
    IY = 1
    IF(INCX < 0)IX = (-N+1)*INCX + 1
    IF(INCY < 0)IY = (-N+1)*INCY + 1
    DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
    10 end do
    RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.

    20 M = MOD(N,7)
    IF( M == 0 ) go to 40
    DO 30 I = 1,M
        DY(I) = DX(I)
    30 end do
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
    50 end do
    RETURN

!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.

    60 CONTINUE
    NS=N*INCX
    DO 70 I=1,NS,INCX
        DY(I) = DX(I)
    70 end do
    RETURN
    END SUBROUTINE DCOPY

    subroutine vector_product(p1, p2, p3, dnorm3)
!            
! Calculates vector product and norm of resulting vector
!
    real(8), intent(in)  :: p1(3), p2(3) 
    real(8), intent(out) :: p3(3), dnorm3

    p3(1) = p1(2)*p2(3) - p1(3)*p2(2)
    p3(2) = p1(3)*p2(1) - p1(1)*p2(3)
    p3(3) = p1(1)*p2(2) - p1(2)*p2(1)
    dnorm3 = sqrt(p3(1)*p3(1) + p3(2)*p3(2) + p3(3)*p3(3))
    
    end subroutine vector_product

    end module pedra_dblas
