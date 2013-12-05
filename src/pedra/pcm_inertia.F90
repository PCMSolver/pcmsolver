!
!    WLKDIN on the basis of ABACUS abawalk.F
!    written and/or modified by Krzysztof Mozgawa 
!    and Ville Weijo, 2010
! 
!    051213 revision. Remove call to LAPACK diagonalization
!    subroutine. Replace with a custom function specialized
!    for 3x3 matrices (thus avoiding having to link BLAS/LAPACK)
!

subroutine pcm_wlkdin(cor, tmass, n, angmom, tinert, omega, cepval, cepvec, docopy, planar, linear)

  implicit none
  real(8), dimension(n,3), intent(in) :: cor
  real(8), dimension(n), intent(in) :: tmass
  integer, intent(in) :: n
  real(8), dimension(3), intent(in) :: angmom
  real(8), dimension(3,3), intent(out) :: tinert, cepvec
  real(8), dimension(3), intent(out) :: omega, cepval
  logical, intent(out) :: planar, linear
  logical, intent(in) :: docopy
  
  integer :: i, j, k, info
  real(8) :: eigval(3), eigvec(3,3), ieigval(3), work(15), tinver(3,3), eigvalinv(3,3), temp(3,3)
  real(8), parameter :: tstlin = 1d-10
  real(8) :: average

  tinert = 0d0 
  do i = 1, n
     do j = 1,3
        do k = 1,3
           tinert(j,k) = tinert(j,k)+ tmass(i)*cor(i,j)*cor(i,k)
        enddo
     enddo
  enddo
! Symmetrize the inertia tensor
  do i = 1, 3
     do j = 1, 3
         average = 0.5d0 * (tinert(i, j) + tinert(j, i))
         tinert(i, j) = average
         tinert(j, i) = average
     enddo
  enddo
  temp = tinert

! Now diagonalize it, we want back both the eigenvalues and eigenvectors
  call dsyevv3(temp, eigvec, eigval)

  if ( abs(eigval(3)-eigval(2)-eigval(1)) .lt. tstlin) then
     planar = .true.
  else
     planar = .false.
  end if
  if (eigval(3) .lt. tstlin) then
     linear = .true.
     planar = .false.
     eigval(3) = 0d0
  else
     linear = .false.
     eigval(3) = 1d0/eigval(3)
  end if
  eigval(2) = 1d0/eigval(2)
  eigval(1) = 1d0/eigval(1)
  eigvalinv = 0d0
  eigvalinv(1,1) = eigval(1)
  eigvalinv(2,2) = eigval(2)
  eigvalinv(3,3) = eigval(3)
  
  tinver= matmul(matmul(eigvec, eigvalinv), transpose(eigvec))

  omega = matmul(tinver, angmom)

  if (docopy) then
     cepval = eigval
     cepvec = eigvec 
  endif


end subroutine pcm_wlkdin
