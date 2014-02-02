    module pedra_utils

    implicit none        

    public around
    public errwrk
    public get_point_group
    public wlkdin

    private

    contains

    subroutine dzero(dx, length)

!...................................................................
! Last revision 5-May-1984 by Hans Jorgen Aa. Jensen
!
!   Subroutine DZERO_ sets a real array of length *LENGTH*
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

    subroutine around(head, print_unit)

    character(len=*), intent(in) :: head
    integer,          intent(in) :: print_unit

    integer :: lhead, lng, ind, i
    
    lhead  = lnblnk(head)
    lng    = lhead + 2
    ind = max(1,(80 - lng)/2 + 1)
    write (print_unit,'(//150a)') (' ',i=1,ind), '+', ('-',i=1,lng), '+'
    write (print_unit,'(150a)')   (' ',i=1,ind), '! ', head(1:lhead), ' !'
    write (print_unit,'(150a)')   (' ',i=1,ind), '+', ('-',i=1,lng), '+'
    write (print_unit,'()')

    end subroutine around

    subroutine errwrk(string, lneed, lavail, print_unit)

! Version 6-Mar-1985 by hjaaj

    character(len=*), intent(in) :: string
    integer,          intent(in) :: lneed
    integer,          intent(in) :: lavail
    integer,          intent(in) :: print_unit
    
    if (lneed >= 0) then
        write(print_unit, 1010) string, lneed, lavail
    else
        write(print_unit, 1020) string, -lneed, lavail
    end if
    stop

    1010 FORMAT(/'  FATAL ERROR, insufficient core for ',A, &
    /T16,'( Need:',I10,', available:',I10,' )')
    1020 FORMAT(/'  FATAL ERROR, insufficient core for ',A, &
    //T16,'Need      :',I10,' + uncalculated amount', &
    /T16,'Available :',I10)

    end subroutine errwrk

    subroutine get_point_group(int_pgroup, char_pgroup, maxrep)

    integer,          intent(in)  :: int_pgroup
    character(len=3), intent(out) :: char_pgroup
    integer,          intent(out) :: maxrep

    integer :: nsymop  ! Number of generators

    select case (int_pgroup)
            case(7)
                   char_pgroup = 'D2h'
                   nsymop      = 3 
            case(6)
                   char_pgroup = 'C2v'
                   nsymop      = 2 
            case(5)
                   char_pgroup = 'D2'
                   nsymop      = 2 
            case(4)
                   char_pgroup = 'C2h'
                   nsymop      = 2
            case(3)
                   char_pgroup = 'Ci'
                   nsymop      = 1
            case(2)
                   char_pgroup = 'Ci'
                   nsymop      = 1
            case(1)
                   char_pgroup = 'Cs'
                   nsymop      = 1
            case(0) 
                   char_pgroup = 'C1'
                   nsymop      = 0
            case default
                    write(6, *) "Can't recognize your group..."
                    stop
    end select
    maxrep      = 2**nsymop - 1

    end subroutine get_point_group

    subroutine wlkdin(cor, tmass, n, angmom, tinert, omega, cepval, cepvec, docopy, planar, linear)
   !
   ! WLKDIN on the basis of ABACUS abawalk.F
   ! written and/or modified by Krzysztof Mozgawa 
   ! and Ville Weijo, 2010
   !
   ! RDR 051213 revision. Remove call to LAPACK diagonalization
   ! subroutine. Replace with a custom function specialized
   ! for 3x3 matrices (thus avoiding having to link BLAS/LAPACK)
   !
  
    use pedra_dblas, only: ddot
    use pedra_dlapack, only: dsyevj3

    integer,    intent(in) :: n
    real(8),    intent(in) :: cor(n, 3)
    real(8),    intent(in) :: tmass(n)
    real(8),    intent(in) :: angmom(3)
    real(8), intent(inout) :: tinert(3, 3), cepvec(3, 3)
    real(8),   intent(out) :: omega(3), cepval(3)
    logical,   intent(out) :: planar, linear
    logical,    intent(in) :: docopy
    
    integer :: i, j, k
    real(8) :: eigval(3), eigvec(3,3), tinver(3,3), eigvalinv(3,3), temp(3,3)
    ! Threshold to convert a numerical zero to a hard zero
    real(8), parameter :: tstlin = 1.0d-05
    real(8) :: average, r2
    
    tinert = 0.0d0 
    eigvec = 0.0d0
    eigval = 0.0d0
    eigvalinv = 0.0d0
    tinver = 0.0d0
    temp = 0.0d0
    r2 = 0.0d0
    ! Build tensor of inertia:
    ! T_ij = sum_l m_l((x_l)_k(x_l)_k delta_ij - (x_l)_i(x_l)_j)
    do i = 1, n ! Loop on centers
       ! Calculate diagonal elements
       r2 = tmass(i) * ddot(3, cor(i, 1), n, cor(i, 1), n)
       do j = 1, 3 ! Loop on coordinates
          ! Build diagonal elements
          tinert(j, j) = tinert(j, j) + r2
          do k = 1, 3 ! Loop on coordinates
             ! Build off-diagonal elements
             tinert(j, k) = tinert(j, k) - tmass(i)*cor(i,j)*cor(i,k)
          enddo
       enddo
    enddo

    ! Symmetrize the inertia tensor
    average = 0.0d0
    do i = 1, 3                                            
       do j = 1, 3
           average = 0.5d0 * (tinert(i, j) + tinert(j, i))
           tinert(i, j) = average
           tinert(j, i) = average
       enddo
    enddo

    temp = tinert
    
    ! Now diagonalize it, we want back both the eigenvalues and eigenvectors
    !call DGEEV('V','N',3,TEMP,3,EIGVAL,IEIGVAL,EIGVEC,3,1,1,WORK,15,INFO)
    call dsyevj3(temp, eigvec, eigval)
    ! Order by decreasing eigenvalue
    call order(eigvec, eigval, 3, 3)
   
    if ( abs(eigval(3)-eigval(2)-eigval(1)) .lt. tstlin) then
       planar = .true.
    else
       planar = .false.
    end if
    if (eigval(3) .lt. tstlin) then
       linear = .true.
       planar = .false.
       eigval(3) = 0.0d0
    else
       linear = .false.
       eigval(3) = 1d0/eigval(3)
    end if
    eigval(2) = 1.0d0/eigval(2)
    eigval(1) = 1.0d0/eigval(1)
    eigvalinv = 0.0d0
    eigvalinv(1,1) = eigval(1)
    eigvalinv(2,2) = eigval(2)
    eigvalinv(3,3) = eigval(3)
    
    tinver= matmul(matmul(eigvec, eigvalinv), transpose(eigvec))
    
    omega = matmul(tinver, angmom)
    
    if (docopy) then
       cepval = eigval
       cepvec = eigvec 
    endif
    
    end subroutine wlkdin

    subroutine order(evec, eval, n, nevec)
!
! Purpose: order the N values in EVAL and their associated vectors
!          in EVEC so EVAL(i+1) .ge. EVAL(i)
!
! Revisions:
!   29-Jul-1992 hjaaj (only dswap if nevec .gt. 0)
!    2-Nov-1984 hjaaj (new parameter NEVEC, EVEC(1:NEVEC,1:N))
!   27-Oct-1984 hjaaj (reduced number of swaps)
!
    use pedra_dblas, only: dswap

    integer,    intent(in) :: n
    real(8), intent(inout) :: evec(n), eval(n)

    integer :: beg, imin, nevec, i, j
    real(8) :: emin
    
    if (n.le.1) return
    beg = 1
    do i=1,n-1
      emin = eval(i)
      imin = i
      do j=i+1,n
        if (eval(j) .lt. emin) then
          emin = eval(j)
          imin = j
        endif
      end do
      if (imin.ne.i) then
        eval(imin)=eval(i)
        eval(i)=emin
        if (nevec .gt. 0) then
          call dswap(nevec,evec(beg),1,evec((imin-1)*nevec+1),1)
        endif
      endif
      beg = beg + nevec
    end do

    end subroutine order

    end module pedra_utils
