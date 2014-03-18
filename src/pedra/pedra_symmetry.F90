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
! MINI-MANUAL
! A. Indexing of symmetry operations and their mapping to a bitstring:
!      zyx         Parity
!   0  000    E      1.0
!   1  001   Oyz    -1.0
!   2  010   Oxz    -1.0
!   3  011   C2z     1.0
!   4  100   Oxy    -1.0
!   5  101   C2y     1.0
!   6  110   C2x     1.0
!   7  111    i     -1.0
! B. Indexing of irreps for the Abelian groups
! C1:  A  <-> 0
! Cs:  A' <-> 0; A'' <-> 1
! C2:  A  <-> 0; B   <-> 1
! Ci:  Ag <-> 0; Au  <-> 1
! C2h: Ag <-> 0; Au  <-> 1; Bu  <-> 2; Bg  <-> 3
! D2:  A  <-> 0; B3  <-> 1; B2  <-> 2; B1  <-> 3
! C2v: A1 <-> 0; B1  <-> 1; B2  <-> 2; A2  <-> 3
! D2h: Ag <-> 0; B3u <-> 1; B2u <-> 2; B1g <-> 3; B1u <-> 4; B2g <-> 5; B3g <-> 6; Au <-> 7
! C. Indexing of isymax matrix
!  The isymax array contains the irrep to which the linear functions
!  (first column) and the rotations (second column) belong. The indexing
!  is given above.
! D. Indexing of the jsop array
!  The jsop array contains the position at which the operation
!  i (index given in point A) appears
!       

    type, public :: point_group
    ! A type containing all you need to know about the point group
 
      ! String with the group name:
      ! C1, Cs, C2, Ci, C2h, D2, C2v, D2h
      character(len=3) :: group_name
      ! Integer representing the group
      ! 0, 1, 2, 3, 4, 5, 6, 7
      integer          :: group_int
      ! Number of generators
      integer          :: nr_generators
      ! Number of not-trivial symmetry operations (2**nr_generators - 1)
      integer          :: maxrep
      ! group%isymax(i, 1): behaviour of principal axes under basic operations
      !     (x-y-z)
      ! group%isymax(i, 2): behaviour of principal rotations under basic operations
      !     (Rx-Ry-Rz)
      integer          :: isymax(3, 2)
      ! Symmetry operations in the Abelian groups.
      ! Bitstring: 1 coordinate changes sign under operation;
      !            0 coordinate does not change sign.
      ! Of course, that's also the binary representation of
      ! numbers from 0 to 7!
      integer          :: jsop(0:7)
      ! 
      integer          :: nr_rotations
      !
      integer          :: nr_reflections
      !
      integer          :: nr_inversion
    end type point_group

    private

    contains

    subroutine get_point_group(pgroup, int_pgroup)
    
    type(point_group), intent(inout) :: pgroup
    integer,              intent(in) :: int_pgroup

    ! Cs
    call build_point_group(1, 1, 0, 0)
    ! C2
    call build_point_group(1, 3, 0, 0)
    ! Ci
    call build_point_group(1, 7, 0, 0)
    ! C2h
    call build_point_group(2, 4, 7, 0)
    ! D2
    call build_point_group(2, 3, 6, 0)
    ! C2v
    call build_point_group(2, 1, 2, 0)
    ! D2h
    call build_point_group(3, 4, 2, 1)
   
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

    real(8) :: pt(0:7)

!
!   PT is the parity of a bitstring:
!     1 for an even number of ones: 000,011,110,101
!    -1 for an odd  number of ones: 001,010,100,111
!
    pt(0) =  1.0d0
    pt(1) = -1.0d0
    pt(2) = -1.0d0
    pt(3) =  1.0d0
    pt(4) = -1.0d0
    pt(5) =  1.0d0
    pt(6) =  1.0d0
    pt(7) = -1.0d0

    get_pt = pt(bit_rep)

    end function get_pt

! Fields: group name, group integer, number of generators,
!         number of nontrivial operations, isymax, jsop,
!         number of rotations, number of reflections,
!         number of inversions.
    
    subroutine init_C1(group)

    type(point_group), intent(inout) :: group
    
    integer :: isymax(3, 2), jsop(0:7)
    
    isymax = reshape([0, 0, 0, 0, 0, 0], [3, 2])
    jsop   = [0, 0, 0, 0, 0, 0, 0, 0]
    
    group = point_group("C1 ", 0, 0, 0, isymax, jsop, 0, 0, 0)
    
    end subroutine init_C1
    
    subroutine init_Cs(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2), jsop(0:7)

    isymax = reshape([0, 0, 1, 1, 1, 0], [3, 2])
    jsop   = [0, 4, 0, 0, 0, 0, 0, 0]

    group = point_group("Cs ", 1, 1, 1, isymax, jsop, 0, 1, 0)
     
    end subroutine init_Cs
    
    subroutine init_C2(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2), jsop(0:7)

    isymax = reshape([1, 1, 0, 1, 1, 0], [3, 2])
    jsop   = [0, 3, 0, 0, 0, 0, 0, 0]

    group = point_group("C2 ", 2, 1, 1, isymax, jsop, 1, 0, 0)
     
    end subroutine init_C2
    
    subroutine init_Ci(group)

    type(point_group), intent(inout) :: group

    integer :: isymax(3, 2), jsop(0:7)
    
    isymax = reshape([1, 1, 1, 0, 0, 0], [3, 2])
    jsop   = [0, 7, 0, 0, 0, 0, 0, 0]

    group = point_group("Ci ", 3, 1, 1, isymax, jsop, 0, 0, 1)
     
    end subroutine init_Ci
    
    subroutine init_C2h(group)

    type(point_group), intent(inout) :: group
    
    integer :: isymax(3, 2), jsop(0:7)

    isymax = reshape([2, 2, 1, 3, 3, 0], [3, 2])
    jsop   = [0, 3, 7, 4, 0, 0, 0, 0]

    group = point_group("C2h", 4, 2, 3, isymax, jsop, 1, 1, 1)
     
    end subroutine init_C2h
    
    subroutine init_D2(group)

    type(point_group), intent(inout) :: group
    
    integer :: isymax(3, 2), jsop(0:7)

    isymax = reshape([1, 2, 3, 1, 2, 3], [3, 2])
    jsop   = [0, 3, 5, 6, 0, 0, 0, 0]

    group = point_group("D2 ", 5, 2, 3, isymax, jsop, 3, 0, 0)
     
    end subroutine init_D2
    
    subroutine init_C2v(group)

    type(point_group), intent(inout) :: group
    
    integer :: isymax(3, 2), jsop(0:7)

    isymax = reshape([1, 2, 0, 2, 1, 3], [3, 2])
    jsop   = [0, 3, 2, 1, 0, 0, 0, 0]

    group = point_group("C2v", 6, 2, 3, isymax, jsop, 1, 2, 0)
     
    end subroutine init_C2v
    
    subroutine init_D2h(group)

    type(point_group), intent(inout) :: group
    
    integer :: isymax(3, 2), jsop(0:7)

    isymax = reshape([1, 2, 4, 6, 5, 3], [3, 2])
    jsop   = [0, 3, 5, 6, 7, 4, 2, 1]

    group = point_group("D2h", 6, 3, 7, isymax, jsop, 3, 3, 1)
     
    end subroutine init_D2h

    subroutine build_point_group(nr_gen, gen1, gen2, gen3)

    integer, intent(in)  :: nr_gen                   
    integer, intent(in)  :: gen1, gen2, gen3
! Local variables
    integer              :: maxrep
    integer              :: isymax(3, 2)
    integer              :: igen(3)
! Integer representation of the rotations bitmaps                   
    integer, parameter   :: irots(3) = [3, 5, 6]
! Integer representation of the reflections bitmaps                   
    integer, parameter   :: irefl(3) = [4, 2, 1] 
! Parity of the symmetry operations bitmaps                   
    integer, parameter   :: jpar(0:7) = [1, -1, -1, 1, -1, 1, 1, -1]
    integer              :: i, j, k, l, i0, i1, i2, ind, ipos
    integer              :: nrots, nrefl, ninvc, igroup, print_unit
    integer              :: char_tab(0:7, 0:7)
    logical              :: lsymop(0:7)
    integer              :: jsop(0:7), ipar(0:7)
    character(3)         :: group, rep(0:7)                   
    character(3), parameter :: groups(0:7) = ['C1 ', 'C2 ', 'Cs ', &
                                              'Ci ', 'D2 ', 'C2v', &
                                              'C2h', 'D2h']
    character(3), parameter :: symop(0:7) = [' E ', 'Oyz', 'Oxz', &
                                             'C2z', 'Oxy', 'C2y', &
                                             'C2x', ' i ']
!
! Builds point group given the generators
!   
    print_unit = 121201
    isymax = 0
    igen = 0
    maxrep = 2**nr_gen - 1
! igen contains the bitmap for the generators                   
    igen = [gen1, gen2, gen3]
! Build isymax(:, 1)
!  determine to which irrep the translations belong to                                     
! x translations
! (igen(1) OR C2x) AND (igen(1) OR C2y) AND (igen(1) OR C2z)
    isymax(1, 1) = iand(iand(ior(igen(1), 6), ior(igen(1), 5)), ior(igen(1), 3))
! y translations
! (igen(2) OR C2x) AND (igen(2) OR C2y) AND (igen(2) OR C2z)
    isymax(2, 1) = iand(iand(ior(igen(2), 6), ior(igen(2), 5)), ior(igen(2), 3))
! z translations
! (igen(3) OR C2x) AND (igen(3) OR C2y) AND (igen(3) OR C2z)
    isymax(3, 1) = iand(iand(ior(igen(3), 6), ior(igen(3), 5)), ior(igen(3), 3))

! Build isymax(:, 2)
!  determine to which irrep the rotations belong to        
!  R_x = (y XOR z) and cyclic permutations 
    isymax(1, 2) = ieor(isymax(2, 1), isymax(3, 1))
    isymax(2, 2) = ieor(isymax(3, 1), isymax(1, 1))
    isymax(3, 2) = ieor(isymax(1, 1), isymax(2, 1))

    do i = 1, nr_gen
      write(print_unit, *) "igen(",i,")=",igen(i)
    end do      
    do i = 1, 3 
      write(print_unit, *) "isymax(",i,",1)=", isymax(i, 1)
    end do      
    do i = 1, 3
      write(print_unit, *) "isymax(",i,",2)=", isymax(i, 2)
    end do
! Build the character table
    lsymop = .false.
! Activate all symmetry operations of the group
    lsymop(0) = .true. 
    jsop(0) = 0
    ipar(0) = 1
    do i = 1, maxrep
      i0 = iand(1, i) * igen(1)
      i1 = iand(1, ishft(i, -1)) * igen(2)
      i2 = iand(1, ishft(i, -2)) * igen(3)
      ind = ieor(ieor(i0, i1),i2)
      lsymop(ind) = .true.
      ipar(i) = jpar(ind)
    end do
! List group operations in preferred order
! Identity, E      
    ind = 0
    jsop(ind) = 0
! Rotations
    nrots = 0
    do i = 1, 3
      if (lsymop(irots(i))) then
        ind = ind + 1
        jsop(ind) = irots(i)
        nrots = nrots + 1
      end if
    end do
! Inversion
    ninvc = 0
    if (lsymop(7)) then
      ind = ind + 1
      jsop(ind) = 7
      ninvc = 1
    end if
! Reflections
    nrefl = 0
    do i = 1, 3
      if (lsymop(irefl(i))) then
        ind = ind + 1
        jsop(ind) = irefl(i)
        nrefl = nrefl + 1
      end if
    end do
! Classify group
! ==============
! tsaue - Here I have devised a highly empirical formula, but it works !!!
    igroup = min(7, nint((4 * nrots + 8 * ninvc + 6 * nrefl) / 3.0))
    group  = groups(igroup)
    char_tab = 0
! Now generate the character table
    do i = 0, maxrep
      ! The character of the identity is always +1 
      char_tab(0, i) = 1
      do j = 1, nr_gen
        char_tab(igen(j), i) = nint(get_pt(iand(ishft(i,-(j-1)), 1)))
        do k = 1, (j-1)
          ind = ieor(igen(j), igen(k))
          char_tab(ind, i) = char_tab(igen(j), i) * char_tab(igen(k), i)
          do l = 1, (k-1)
            char_tab(ieor(ind, igen(l)), i) = char_tab(ind, i) * char_tab(igen(l), i)
          end do
        end do
      end do
    end do
! Classify irrep
    do i = 0, maxrep
      rep(i) = 'A  '
      ipos = 2
! 1. Rotational symmetry
      if (nrots == 3) then
        ind = (1 - char_tab(jsop(1), i)) + (1 - char_tab(jsop(2), i))/2
        if (ind /= 0) then
          rep(i)(1:1) = 'B'
          rep(i)(2:2) = char(ichar('0') + ind)
          ipos = 3
        end if
      else if (nrots == 1) then
        if (char_tab(jsop(1), i) == -1) then
          rep(i)(1:1) = 'B'
        end if
        if (nrefl == 2) then
          if (iand(ishft(jsop(1), -1), 1) == 1) then
            ind = 2
          else
            ind = 3
          end if
          if (char_tab(jsop(ind), i) == 1) then
            rep(i)(2:2) = '1'
          else
            rep(i)(2:2) = '2'
          end if
        end if
      else if (nrefl == 1) then
! 2. Mirror symmetry
          if (char_tab(jsop(1), i) == 1) then
            rep(i)(2:2) = ''''
          else if (char_tab(jsop(1), i) == -1) then
            rep(i)(2:2) = '"'
          end if
      end if
! 3. Inversion symmetry
      if (ninvc == 1) then
        ind = nrots + 1
        if (char_tab(jsop(ind), i) == 1) then
          rep(i)(ipos:ipos) = 'g'
        else
          rep(i)(ipos:ipos) = 'u'
        end if
      end if
    end do  
! Output
! 1. Group name and generators      
    write(print_unit, '(a, a3)') 'Point group: ', group
    if (nr_gen > 0) then
      write(print_unit, '(/3x, a/)') '* The point group was generated by:'
      do i = 1, nr_gen
        if (symop(igen(i))(1:1) == 'C') then
          write(print_unit, '(6x, 3a)') 'Rotation about the ', symop(igen(i))(3:3),'-axis'
        else if (symop(igen(i))(1:1) == 'O') then
          write(print_unit, '(6x, 3a)') 'Reflection in the ', symop(igen(i))(2:3),'-plane'
        else
          write(print_unit, '(6x, a)') 'Inversion center'
        end if
      end do
! 2. Group multiplication table      
      write(print_unit,'(/3x, a/)') '* Group multiplication table'
      write(print_unit,'(8x, a1, 8(1x, a3, 1x))') '|', (symop(jsop(i)), i = 0, maxrep)
      write(print_unit,'(3x,a6,8a5)') '-----+', ('-----', i = 0, maxrep)
      do i = 0, maxrep
        write(print_unit,'(4x, a3, 1x, a1, 8(1x, a3, 1x))') symop(jsop(i)), '|', &
                  (symop(ieor(jsop(i), jsop(j))), j = 0, maxrep)
      end do
! 3. Character table
     write(print_unit,'(/3x, a/)') '* Character table'
     write(print_unit,'(8x, a1, 8(1x, a3, 1x))') '|', (symop(jsop(j)), j = 0, maxrep)
     write(print_unit,'(3x, a6, 8a5)') '-----+', ('-----', i = 0, maxrep)
     do i = 0, maxrep
       write(print_unit,'(4x, a3, 1x, a1, 8(1x, i3, 1x))') rep(i), '|', (char_tab(jsop(j), i), j=0, maxrep)
     end do
! 4. Direct product table
      write(print_unit,'(/3x, a/)') '* Direct product table'
      write(print_unit,'(8x, a1, 8(1x, a3, 1x))') '|', (rep(i), i = 0, maxrep)
      write(print_unit,'(3x, a6, 8a5)') '-----+', ('-----', i = 0, maxrep)
      do i = 0, maxrep
        write(print_unit,'(3x, 1x, a3, 1x, a1, 8(1x, a3, 1x))') rep(i), '|', (rep(ieor(i, j)), j = 0, maxrep)
      end do
    end if

    end subroutine build_point_group

    end module pedra_symmetry
