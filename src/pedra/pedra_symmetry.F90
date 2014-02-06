    module pedra_symmetry
! NOTE: the ibtfun module is gone, as we can safely use Fortran
!       standard intrisic functions. 
!       Mapping between intrinsics and ibtfun:
!           ibtand(i, j) <-> iand(i, j)             
!           ibtor(i, j)  <-> ior(i, j)
!           ibtshl(i, j) <-> ishft(i, j)
!           ibtshr(i, j) <-> ishft(i, -j) !WARNING!
!           ibtxor(i, j) <-> ieor(i, j)

    use, intrinsic :: iso_c_binding

    implicit none

    public get_point_group
    public get_pt

    type, public :: point_group
    ! A type containing all you need to know about the point group
 
      ! String with the group name:
      ! C1, Cs, C2, Ci, C2h, D2, C2v, D2h
      character(len=3) :: group_name
      ! Integer representing the group
      ! 0, 1, 2, 3, 4, 5, 6, 7
      integer          :: group_int
      ! Number of generators
      integer          :: nr_generators = 0
      ! Number of not-trivial symmetry operations (2**nr_generators - 1)
      integer          :: maxrep = 0
      ! group%isymax(i, 1): behaviour of principal axes under basic operations
      !     (x-y-z)
      ! group%isymax(i, 2): behaviour of principal rotations under basic operations
      !     (Rx-Ry-Rz)
      integer          :: isymax(3, 2) = 0
      ! Symmetry operations in the Abelian groups.
      ! Bitstring: 1 coordinate changes sign under operation;
      !            0 coordinate does not change sign.
      ! Of course, that's also the binary representation of
      ! numbers from 0 to 7!
      !      zyx 
      !   0  000    E
      !   1  001   Oyz
      !   2  010   Oxz
      !   3  011   C2z
      !   4  100   Oxy
      !   5  101   C2y
      !   6  110   C2x
      !   7  111    i
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
   
    call init_pt
                   
    select case (int_pgroup)
            case(0)
                   ! Generators: none
                   call init_C1(pgroup)
            case(1)
                   ! Generators: Z
                   call init_Cs(pgroup)
            case(2)
                   ! Generators: XY 
                   call init_C2(pgroup)
            case(3)
                   ! Generators: XYZ 
                   call init_Ci(pgroup)
            case(4)
                   ! Generators: Z XY
                   call init_C2h(pgroup)
            case(5)
                   ! Generators: XZ YZ
                   call init_D2(pgroup)
            case(6)
                   ! Generators: X Y
                   call init_C2v(pgroup)
            case(7)
                   ! Generators: X Y Z
                   call init_D2h(pgroup)
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
    
    subroutine init_C1(group)

    type(point_group), intent(inout) :: group
    
    integer :: isymax(3, 2) = reshape([0, 0, 0, 0, 0, 0], [3, 2])
    integer :: jsop(0:7) = [0, 0, 0, 0, 0, 0, 0, 0]
    
    group = point_group("C1 ", 0, 0, 0, isymax, jsop, 0, 0, 0)
    
    end subroutine init_C1
    
    subroutine init_Cs(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2) = reshape([0, 0, 1, 1, 1, 0], [3, 2])
    integer :: jsop(0:7) = [0, 4, 0, 0, 0, 0, 0, 0]

    group = point_group("Cs ", 1, 1, 1, isymax, jsop, 0, 1, 0)
     
    end subroutine init_Cs
    
    subroutine init_C2(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2) = reshape([1, 1, 0, 1, 1, 0], [3, 2])
    integer :: jsop(0:7) = [0, 3, 0, 0, 0, 0, 0, 0]

    group = point_group("C2 ", 2, 1, 1, isymax, jsop, 1, 0, 0)
     
    end subroutine init_C2
    
    subroutine init_Ci(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2) = reshape([1, 1, 1, 0, 0, 0], [3, 2])
    integer :: jsop(0:7) = [0, 7, 0, 0, 0, 0, 0, 0]

    group = point_group("Ci ", 3, 1, 1, isymax, jsop, 0, 0, 1)
     
    end subroutine init_Ci
    
    subroutine init_C2h(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2) = reshape([2, 2, 1, 3, 3, 0], [3, 2])
    integer :: jsop(0:7) = [0, 3, 7, 4, 0, 0, 0, 0]

    group = point_group("C2h", 4, 2, 3, isymax, jsop, 1, 1, 1)
     
    end subroutine init_C2h
    
    subroutine init_D2(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2) = reshape([1, 2, 3, 1, 2, 3], [3, 2])
    integer :: jsop(0:7) = [0, 3, 5, 6, 0, 0, 0, 0]

    group = point_group("D2 ", 5, 2, 3, isymax, jsop, 3, 0, 0)
     
    end subroutine init_D2
    
    subroutine init_C2v(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2) = reshape([1, 2, 0, 2, 1, 3], [3, 2])
    integer :: jsop(0:7) = [0, 3, 2, 1, 0, 0, 0, 0]

    group = point_group("C2v", 6, 2, 3, isymax, jsop, 1, 2, 0)
     
    end subroutine init_C2v
    
    subroutine init_D2h(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2) = reshape([1, 2, 4, 6, 5, 3], [3, 2])
    integer :: jsop(0:7) = [0, 3, 5, 6, 7, 4, 2, 1]

    group = point_group("D2h", 6, 3, 7, isymax, jsop, 3, 3, 1)
     
    end subroutine init_D2h

    end module pedra_symmetry
