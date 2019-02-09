!
! PCMSolver, an API for the Polarizable Continuum Model
! Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
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

module pcmsolver

  use, intrinsic :: iso_c_binding

  implicit none

  type, bind(C) :: PCMInput
    character(kind=c_char, len=1) :: cavity_type(8)
    integer(c_int)                :: patch_level = 0
    real(c_double)                :: coarsity = 0.0
    real(c_double)                :: area = 0.0
    character(kind=c_char, len=1) :: radii_set(8)
    real(c_double)                :: min_distance = 0.0
    integer(c_int)                :: der_order = 0
    logical(c_bool)               :: scaling = .false.
    character(kind=c_char, len=1) :: restart_name(20)
    real(c_double)                :: min_radius = 0.0
    character(kind=c_char, len=1) :: solver_type(7)
    real(c_double)                :: correction = 0.0
    character(kind=c_char, len=1) :: solvent(16)
    real(c_double)                :: probe_radius = 0.0
    character(kind=c_char, len=1) :: equation_type(11)
    character(kind=c_char, len=1) :: inside_type(7)
    real(c_double)                :: outside_epsilon = 0.0
    character(kind=c_char, len=1) :: outside_type(22)
  end type PCMInput

  enum, bind(C)
    enumerator :: PCMSOLVER_READER_OWN = 0, PCMSOLVER_READER_HOST = 1
  end enum

interface pcmsolver_new
  function pcmsolver_new(input_reading, nr_nuclei, charges, coordinates, symmetry_info, &
                         host_input, writer) result(context) bind(C)
    import
    integer(c_int), intent(in), value :: input_reading
    integer(c_int), intent(in), value :: nr_nuclei
    real(c_double), intent(in)        :: charges(*)
    real(c_double), intent(in)        :: coordinates(*)
    integer(c_int), intent(in)        :: symmetry_info(*)
    type(PCMInput), intent(in)        :: host_input
    type(c_funptr), intent(in), value :: writer
    type(c_ptr) :: context
  end function
end interface

interface pcmsolver_new_v1112
  function pcmsolver_new_v1112(input_reading, nr_nuclei, charges, coordinates, symmetry_info, &
                               parsed_fname, host_input, writer) result(context) bind(C)
    import
    integer(c_int), intent(in), value :: input_reading
    integer(c_int), intent(in), value :: nr_nuclei
    real(c_double), intent(in)        :: charges(*)
    real(c_double), intent(in)        :: coordinates(*)
    integer(c_int), intent(in)        :: symmetry_info(*)
    character(kind=c_char, len=1), intent(in) :: parsed_fname(*)
    type(PCMInput), intent(in)        :: host_input
    type(c_funptr), intent(in), value :: writer
    type(c_ptr) :: context
  end function
end interface

interface pcmsolver_delete
  subroutine pcmsolver_delete(context) bind(C)
    import
    type(c_ptr), value :: context
  end subroutine
end interface

interface pcmsolver_is_compatible_library
  function pcmsolver_is_compatible_library() result(compatible) bind(C)
    import
    logical(c_bool) :: compatible
  end function
end interface

interface pcmsolver_print
  subroutine pcmsolver_print(context) bind(C)
    import
    type(c_ptr), value :: context
  end subroutine
end interface

interface pcmsolver_citation
   subroutine pcmsolver_citation(writer) bind(C)
     import
     type(c_funptr), intent(in), value :: writer
   end subroutine
end interface

interface pcmsolver_get_cavity_size
  function pcmsolver_get_cavity_size(context) result(nr_points) bind(C)
    import
    type(c_ptr), value :: context
    integer(c_int)  :: nr_points
  end function
end interface

interface pcmsolver_get_irreducible_cavity_size
  function pcmsolver_get_irreducible_cavity_size(context) result(nr_points_irr) bind(C)
    import
    type(c_ptr), value :: context
    integer(c_int)  :: nr_points_irr
  end function
end interface

interface pcmsolver_get_centers
  subroutine pcmsolver_get_centers(context, centers) bind(C)
    import
    type(c_ptr), value :: context
    real(c_double), intent(inout) :: centers(*)
  end subroutine
end interface

interface pcmsolver_get_center
  subroutine pcmsolver_get_center(context, its, center) bind(C)
    import
    type(c_ptr), value :: context
    integer(c_int), value, intent(in) :: its
    real(c_double), intent(inout) :: center(*)
  end subroutine
end interface

interface pcmsolver_get_areas
  subroutine pcmsolver_get_areas(context, areas) bind(C)
    import
    type(c_ptr), value :: context
    real(c_double), intent(inout) :: areas(*)
  end subroutine
end interface

interface pcmsolver_compute_asc
  subroutine pcmsolver_compute_asc(context, mep_name, asc_name, irrep) bind(C)
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: mep_name(*), asc_name(*)
    integer(c_int), value, intent(in) :: irrep
  end subroutine
end interface

interface pcmsolver_compute_response_asc
  subroutine pcmsolver_compute_response_asc(context, mep_name, asc_name, irrep) bind(C)
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: mep_name(*), asc_name(*)
    integer(c_int), value, intent(in) :: irrep
  end subroutine
end interface

interface pcmsolver_compute_polarization_energy
  function pcmsolver_compute_polarization_energy(context, mep_name, asc_name) result(energy) bind(C)
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: mep_name(*), asc_name(*)
    real(c_double) :: energy
  end function
end interface

interface pcmsolver_get_asc_dipole
  function pcmsolver_get_asc_dipole(context, asc_name, dipole) result(mu) bind(C)
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: asc_name(*)
    real(c_double), intent(inout) :: dipole(*)
    real(c_double) :: mu
  end function
end interface

interface pcmsolver_get_surface_function
  subroutine pcmsolver_get_surface_function(context, f_size, values, name) bind(C)
    import
    type(c_ptr), value :: context
    integer(c_int), value, intent(in) :: f_size
    real(c_double), intent(inout) :: values(*)
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine
end interface

interface pcmsolver_set_surface_function
  subroutine pcmsolver_set_surface_function(context, f_size, values, name) bind(C)
    import
    type(c_ptr), value :: context
    integer(c_int), value, intent(in) :: f_size
    real(c_double), intent(in) :: values(*)
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine
end interface

interface pcmsolver_print_surface_function
  subroutine pcmsolver_print_surface_function(context, name) bind(C)
    import
    type(c_ptr), value :: context
    character(c_char), intent(in) :: name(*)
  end subroutine
end interface

interface pcmsolver_save_surface_functions
  subroutine pcmsolver_save_surface_functions(context) bind(C)
    import
    type(c_ptr), value :: context
  end subroutine
end interface

interface pcmsolver_save_surface_function
  subroutine pcmsolver_save_surface_function(context, name) bind(C)
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine
end interface

interface pcmsolver_load_surface_function
  subroutine pcmsolver_load_surface_function(context, name) bind(C)
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine
end interface

interface pcmsolver_write_timings
  subroutine pcmsolver_write_timings(context) bind(C)
    import
    type(c_ptr), value :: context
  end subroutine
end interface

contains

  ! \brief Convert a Fortran string into a C string.
  ! \param[in] string_f03 a Fortran character string.
  ! \return array_c Null-terminated C string in a character array.
  pure function pcmsolver_fstring_to_carray(string_f03) result(array_c)
    character(len=*), intent(in) :: string_f03
    character(kind=c_char, len=1) :: array_c(len(string_f03)+1)

    integer :: i

    do i = 1, len(string_f03)
        array_c(i) = string_f03(i:i)
    end do
    array_c(i) = c_null_char
  end function

end module pcmsolver
