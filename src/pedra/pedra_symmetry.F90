    module pedra_symmetry
! NOTE: the ibtfun module is gone, as we can safely use Fortran
!       standard intrisic functions. 
!       Mapping between intrinsics and ibtfun:
!           ibtand(i, j) <-> iand(i, j)             
!           ibtor(i, j)  <-> ior(i, j)
!           ibtshl(i, j) <-> ishft(i, j)
!           ibtshr(i, j) <-> ishft(i, -j) !WARNING!
!           ibtxor(i, j) <-> ieor(i, j)

    implicit none

    public get_point_group
    public get_pt

    type, public :: point_group
    ! A type containing all you need to know about the point group
 
      ! String with the group name:                                           
      ! D2h, C2v, D2, C2h, Ci, C2, Cs, C1
      character(len=3) :: group_name = 'C1'
      ! Integer representing the group
      ! 7, 6, 5, 4, 3, 2, 1, 0
      integer          :: group_int = 0
      ! Number of generators
      integer          :: nr_generators = 0
      ! maxrep = 2**nr_generators - 1
      integer          :: maxrep = 0 
      ! group%isymax(i, 1): behaviour of principal axes under basic operations
      !     (x-y-z)
      ! group%isymax(i, 2): behaviour of principal rotations under basic operations
      !     (Rx-Ry-Rz)
      integer          :: isymax(3, 2) = 0
      ! Symmetry operations in the Abelian groups
      !      zyx 
      !   0  000    E
      !   1  001   Oyz
      !   2  010   Oxz
      !   3  011   C2z
      !   4  100   Oxy
      !   5  101   C2y
      !   6  110   C2x
      !   7  111    i
      logical          :: lsymop(0:7) = .false.
      !
      integer          :: jsop(0:7) = 0
      ! 
      integer          :: nr_rotations = 0
      !
      integer          :: nr_reflections = 0
      !
      integer          :: nr_inversion = 0
    end type point_group

    private

    real(8) :: pt(0:7) = 0.0d0   ! Parity of a bitstring

    contains

    subroutine get_point_group(pgroup, int_pgroup)
    
    type(point_group), intent(inout) :: pgroup
    integer,              intent(in) :: int_pgroup
   
    ! OK, this is a dirty, hacky patch and my eyes are bleeding for
    ! it... Will fix as soon as I have a better idea
    character(len=1) :: kasym(3, 3) = ' '

    call init_pt

    select case (int_pgroup)
            case(7)
                   ! Generators: X Y Z
                   pgroup%group_name = 'D2h'
                   pgroup%group_int  = int_pgroup
                   pgroup%nr_generators = 3 
                   pgroup%maxrep = 2**(pgroup%nr_generators) - 1
                   kasym(1, 1) = 'X'
                   kasym(1, 2) = 'Y'
                   kasym(1, 3) = 'Z'
                   call init_group(pgroup, kasym)
            case(6)
                   ! Generators: X Y
                   pgroup%group_name = 'C2v'
                   pgroup%group_int  = int_pgroup
                   pgroup%nr_generators = 2
                   pgroup%maxrep = 2**(pgroup%nr_generators) - 1
                   kasym(1, 1) = 'X'
                   kasym(1, 2) = 'Y'
                   call init_group(pgroup, kasym)
            case(5)
                   ! Generators: XZ YZ
                   pgroup%group_name = 'D2'
                   pgroup%group_int  = int_pgroup
                   pgroup%nr_generators = 2
                   pgroup%maxrep = 2**(pgroup%nr_generators) - 1
                   kasym(1, 1) = 'X'
                   kasym(1, 2) = 'Y'
                   kasym(2, 1) = 'Z'
                   kasym(2, 2) = 'Z'
                   call init_group(pgroup, kasym)
            case(4)
                   ! Generators: Z XY
                   pgroup%group_name = 'C2h'
                   pgroup%group_int  = int_pgroup
                   pgroup%nr_generators = 2
                   pgroup%maxrep = 2**(pgroup%nr_generators) - 1
                   kasym(1, 1) = 'Z'
                   kasym(1, 2) = 'X'
                   kasym(2, 2) = 'Y'
                   call init_group(pgroup, kasym)
            case(3)
                   ! Generators: XYZ 
                   pgroup%group_name = 'Ci'
                   pgroup%group_int  = int_pgroup
                   pgroup%nr_generators = 1
                   pgroup%maxrep = 2**(pgroup%nr_generators) - 1
                   kasym(1, 1) = 'X'
                   kasym(2, 1) = 'Y'
                   kasym(3, 1) = 'Z'
                   call init_group(pgroup, kasym)
            case(2)
                   ! Generators: XY 
                   pgroup%group_name = 'C2'
                   pgroup%group_int  = int_pgroup
                   pgroup%nr_generators = 1
                   pgroup%maxrep = 2**(pgroup%nr_generators) - 1
                   kasym(1, 1) = 'X'
                   kasym(2, 1) = 'Y'
                   call init_group(pgroup, kasym)
            case(1)
                   ! Generators: Z
                   pgroup%group_name = 'Cs'
                   pgroup%group_int  = int_pgroup
                   pgroup%nr_generators = 1
                   pgroup%maxrep = 2**(pgroup%nr_generators) - 1
                   kasym(1, 1) = 'Z' 
                   call init_group(pgroup, kasym)
            case(0)
                   ! Generators: none
                   pgroup%group_name = 'C1'
                   pgroup%group_int  = int_pgroup
                   pgroup%nr_generators = 0
                   pgroup%maxrep = 2**(pgroup%nr_generators) - 1
                   call init_group(pgroup, kasym)
            case default
                    write(6, *) "Can't recognize your group..."
                    stop
    end select

    end subroutine get_point_group

    real(8) function get_pt(bit_rep)

    integer, intent(in) :: bit_rep

    get_pt = pt(bit_rep)

    end function get_pt

    subroutine init_pt()
!
!     PT is the parity of a bitstring:
!       1 for an even number of ones: 000,011,110,101
!      -1 for an odd  number of ones: 001,010,100,111
!
      pt(0) =  1.0d0
      pt(1) = -1.0d0
      pt(2) = -1.0d0
      pt(3) =  1.0d0
      pt(4) = -1.0d0
      pt(5) =  1.0d0
      pt(6) =  1.0d0
      pt(7) = -1.0d0

    end subroutine init_pt

    subroutine init_group(group, kasym)

    type(point_group), intent(inout) :: group
    character(len=1),     intent(in) :: kasym(3, 3)

    integer            :: i, j, k, iaxis, i0, i1, i2, ind
    integer            :: igen(3)
    integer :: irots(3) = [3, 5, 6] ! The third, fifth and sixth ops
    integer :: irefl(3) = [4, 2, 1] ! The fourth, second and first ops

    do i = 1, 3
       group%isymax(i, 1) = 0
       igen(i)      = 0
    end do

!
!   Determine:
!   igen(i)   - basic operations
!   group%isymax(i, 1) - behavior of principal axes under basic operations
!   ================================================================
!
    iaxis = 0
    if (group%nr_generators > 0) then
       write(121201, '(a)') "Symmetry operations"
       write(121201, '(a, i2)') "  Number of generators = ", group%nr_generators 
    end if
    do j = 1, group%nr_generators 
      do i = 1, 3
        if (kasym(i, j) /= ' ') then
          k = ichar(kasym(i, j)) - ichar('W') ! This letter needs to be uppercase 
          igen(j)      = igen(j)      + 2**(k-1)
          group%isymax(k, 1) = group%isymax(k, 1) + 2**(j-1)
        end if
      end do
      iaxis = ior(iaxis, igen(j))
    end do 
!
!   Determine:
!   group%isymax(i, 2) - behaviour of principal rotations under basic operations
!   ======================================================================
!
    group%isymax(1, 2) = ieor(group%isymax(2, 1), group%isymax(3, 1))
    group%isymax(2, 2) = ieor(group%isymax(3, 1), group%isymax(1, 1))
    group%isymax(3, 2) = ieor(group%isymax(1, 1), group%isymax(2, 1))

    group%lsymop(0) = .true.
    group%jsop(0) = 0
    do i = 1, group%maxrep 
        i0  = iand(1,i)             * igen(1)
        i1  = iand(1, ishft(i, -1)) * igen(2)
        i2  = iand(1, ishft(i, -2)) * igen(3)
        ind = ieor(ieor(i0, i1), i2)
        group%lsymop(ind) = .true.
    end do
!
!   Identity
!
    ind = 0
    group%jsop(ind) = 0
!
!   Rotations
!
    group%nr_rotations = 0
    do i = 1, 3
       if(group%lsymop(irots(i))) then
          ind          = ind + 1
          group%jsop(ind)    = irots(i)
          group%nr_rotations = group%nr_rotations + 1
       endif
    end do
!
!   Inversion
!
    group%nr_inversion = 0
    if (group%lsymop(7)) then
      ind          = ind + 1
      group%jsop(ind)    = 7
      group%nr_inversion = 1
    endif
!
!   Reflections
!
    group%nr_reflections = 0
    do i = 1,3
       if (group%lsymop(irefl(i))) then             
         ind            = ind + 1
         group%jsop(ind)      = irefl(i)
         group%nr_reflections = group%nr_reflections + 1
       endif
    end do

    end subroutine init_group
     
    end module pedra_symmetry
