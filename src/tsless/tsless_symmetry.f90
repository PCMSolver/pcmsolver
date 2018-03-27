!
! PCMSolver, an API for the Polarizable Continuum Model
! Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
!
! This file is part of PCMSolver.
!
! PCMSolver is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PCMSolver is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
!
! For information on the complete list of contributors to the
! PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
!

module tsless_symmetry
! NOTE: the ibtfun module is gone, as we can safely use Fortran
!       standard intrinsic functions.
!       Mapping between intrinsics and ibtfun:
!           ibtand(i, j) <-> iand(i, j)
!           ibtor(i, j)  <-> ior(i, j)
!           ibtshl(i, j) <-> ishft(i, j)
!           ibtshr(i, j) <-> ishft(i, -j) !WARNING!
!           ibtxor(i, j) <-> ieor(i, j)

use, intrinsic :: iso_c_binding
use tsless_precision

implicit none

public get_pt
public build_point_group
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
! B. Indexing of isymax matrix
!  The isymax array contains the irrep to which the linear functions
!  (first column) and the rotations (second column) belong. The indexing
!  is given above.
! C. Indexing of the jsop array
!  The jsop array contains the position at which the operation
!  i (index given in point A) appears
!

type, public :: point_group
! A type containing all you need to know about the point group

  ! String with the group name:
  ! C1, C2, Cs, Ci, D2, C2v, C2h, D2h
  character(len=3) :: group_name
  ! Integer representing the group
  ! 0, 1, 2, 3, 4, 5, 6, 7
  integer(kind=regint_k)          :: group_int
  ! Number of generators
  integer(kind=regint_k)          :: nr_generators
  ! Number of not-trivial symmetry operations (2**nr_generators - 1)
  integer(kind=regint_k)          :: maxrep
  ! group%isymax(i, 1): behaviour of principal axes under basic operations
  !     (x-y-z)
  ! group%isymax(i, 2): behaviour of principal rotations under basic operations
  !     (Rx-Ry-Rz)
  integer(kind=regint_k)          :: isymax(3, 2)
  ! Symmetry operations in the Abelian groups.
  ! Bitstring: 1 coordinate changes sign under operation;
  !            0 coordinate does not change sign.
  ! Of course, that's also the binary representation of
  ! numbers from 0 to 7!
  integer(kind=regint_k)          :: jsop(0:7)
  !
  integer(kind=regint_k)          :: nr_rotations
  !
  integer(kind=regint_k)          :: nr_reflections
  !
  integer(kind=regint_k)          :: nr_inversion
end type point_group

private
contains

    !> \brief returns parity of a bitstring
    !> \param[in] bit_rep bitstring represenation of symmetry operator
    !> PT is the parity of a bitstring:
    !>   1 for an even number of ones: 000,011,110,101
    !>  -1 for an odd  number of ones: 001,010,100,111
    real(kind=dp) function get_pt(bit_rep)

    integer(kind=regint_k), intent(in) :: bit_rep

    real(kind=dp) :: pt(0:7)

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

    !> \brief Builds point group given the generators
    !> \param[in] nr_gen number of generators
    !> \param[in] gen1 first generator
    !> \param[in] gen2 second generator
    !> \param[in] gen3 third generator
    !> \param[in] print_unit logical unit for printing
    !> Originally written by Trond Saue for DALTON/DIRAC
    !> Copied and adapted by Roberto Di Remigio
    !> The generators are encoded as bitstrings
    function build_point_group(nr_gen, gen1, gen2, gen3, printer) result(pg)

    use, intrinsic :: iso_fortran_env, only: output_unit

    !> Passed variables
    integer(kind=regint_k), intent(in) :: nr_gen
    integer(kind=regint_k), intent(in) :: gen1, gen2, gen3
    integer, optional,      intent(in) :: printer
    !> Output variables
    type(point_group) :: pg
    !> Local variables
    integer                             :: print_out
    integer(kind=regint_k)              :: maxrep
    integer(kind=regint_k)              :: isymax(3, 2)
    integer(kind=regint_k)              :: igen(3)
    !> Integer representation of the rotations bitmaps
    integer(kind=regint_k), parameter   :: irots(3) = [3, 5, 6]
    integer(kind=regint_k), parameter   :: rots(3) = [6, 5, 3]
    !> Integer representation of the reflections bitmaps
    integer(kind=regint_k), parameter   :: irefl(3) = [4, 2, 1]
    !> Parity of the symmetry operations bitmaps
    integer(kind=regint_k), parameter   :: jpar(0:7) = [1, -1, -1, 1, -1, 1, 1, -1]
    integer(kind=regint_k)              :: i, j, k, l, i0, i1, i2, ind, ipos, bitmap
    integer(kind=regint_k)              :: nrots, nrefl, ninvc, igroup
    integer(kind=regint_k)              :: char_tab(0:7, 0:7)
    logical                             :: lsymop(0:7)
    integer(kind=regint_k)              :: jsop(0:7), ipar(0:7)
    character(3)                        :: group, rep(0:7)
    character(3), parameter :: groups(0:7) = ['C1 ', 'C2 ', 'Cs ', &
                                              'Ci ', 'D2 ', 'C2v', &
                                              'C2h', 'D2h']
    character(3), parameter :: symop(0:7) = [' E ', 'Oyz', 'Oxz', &
                                             'C2z', 'Oxy', 'C2y', &
                                             'C2x', ' i ']

    if (present(printer)) then
       print_out = printer
    else
       print_out = output_unit
    end if

    isymax = 0
    igen = 0
    maxrep = 2**nr_gen - 1
    ! igen contains the bitmap for the generators
    igen = [gen1, gen2, gen3]
    ! Build isymax(:, 1)
    !  determine to which irrep the translations belong to
    ! Loop over Cartesian axes
    do i = 1, 3
      bitmap = 0
      ! Loop over generators
      do j = 1, nr_gen
        ! Apply generators on Cartesian axes rots(i) and check the character
        if (nint(get_pt(ior(igen(j), rots(i)))) == -1) then
          ! Set the bitmap
          bitmap = ibset(bitmap, j)
        end if
      end do
      ! Right-shift the bitmap and assign to isymax
      isymax(i, 1) = ishft(bitmap, -1)
    end do

    ! Build isymax(:, 2)
    !  determine to which irrep the rotations belong to
    !  R_x = (y XOR z) and cyclic permutations
    isymax(1, 2) = ieor(isymax(2, 1), isymax(3, 1))
    isymax(2, 2) = ieor(isymax(3, 1), isymax(1, 1))
    isymax(3, 2) = ieor(isymax(1, 1), isymax(2, 1))

    ! Build the character table
    lsymop = .false.
    ! Activate all symmetry operations of the group
    lsymop(0) = .true.
    jsop(0) = 0
    ipar(0) = 1
    do i = 1, maxrep
      i0 = iand(1_regint_k, i) * igen(1)
      i1 = iand(1_regint_k, ishft(i, -1)) * igen(2)
      i2 = iand(1_regint_k, ishft(i, -2)) * igen(3)
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
        char_tab(igen(j), i) = nint(get_pt(iand(ishft(i,-(j-1)), 1_regint_k)))
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
          if (iand(ishft(jsop(1), -1), 1_regint_k) == 1) then
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
    write(print_out, '(a, a3)') 'Point group: ', group
    if (nr_gen > 0) then
      write(print_out, '(/3x, a/)') '* The point group was generated by:'
      do i = 1, nr_gen
        if (symop(igen(i))(1:1) == 'C') then
          write(print_out, '(6x, 3a)') 'Rotation about the ', symop(igen(i))(3:3),'-axis'
        else if (symop(igen(i))(1:1) == 'O') then
          write(print_out, '(6x, 3a)') 'Reflection in the ', symop(igen(i))(2:3),'-plane'
        else
          write(print_out, '(6x, a)') 'Inversion center'
        end if
      end do
      ! 2. Group multiplication table
      write(print_out,'(/3x, a/)') '* Group multiplication table'
      write(print_out,'(8x, a1, 8(1x, a3, 1x))') '|', (symop(jsop(i)), i = 0, maxrep)
      write(print_out,'(3x,a6,8a5)') '-----+', ('-----', i = 0, maxrep)
      do i = 0, maxrep
        write(print_out,'(4x, a3, 1x, a1, 8(1x, a3, 1x))') symop(jsop(i)), '|', &
                  (symop(ieor(jsop(i), jsop(j))), j = 0, maxrep)
      end do
      ! 3. Character table
      write(print_out,'(/3x, a/)') '* Character table'
      write(print_out,'(8x, a1, 8(1x, a3, 1x))') '|', (symop(jsop(j)), j = 0, maxrep)
      write(print_out,'(3x, a6, 8a5)') '-----+', ('-----', i = 0, maxrep)
      do i = 0, maxrep
        write(print_out,'(4x, a3, 1x, a1, 8(1x, i3, 1x))') rep(i), '|', (char_tab(jsop(j), i), j=0, maxrep)
      end do
      ! 4. Direct product table
      write(print_out,'(/3x, a/)') '* Direct product table'
      write(print_out,'(8x, a1, 8(1x, a3, 1x))') '|', (rep(i), i = 0, maxrep)
      write(print_out,'(3x, a6, 8a5)') '-----+', ('-----', i = 0, maxrep)
      do i = 0, maxrep
        write(print_out,'(3x, 1x, a3, 1x, a1, 8(1x, a3, 1x))') rep(i), '|', (rep(ieor(i, j)), j = 0, maxrep)
      end do
    end if
    ! Fields: group name, group integer(kind=regint_k), number of generators,
    !         number of nontrivial operations, isymax, jsop,
    !         number of rotations, number of reflections,
    !         number of inversions.
    pg = point_group(group, igroup, nr_gen, maxrep, isymax, jsop, nrots, nrefl, ninvc)

    end function build_point_group

    end module tsless_symmetry
