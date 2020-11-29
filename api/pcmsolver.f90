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

  private :: pcmsolver_fstring_to_carray

interface
  function pcmsolver_new(input_reading, nr_nuclei, charges, coordinates, symmetry_info, &
       host_input, writer) result(context) &
       bind(C)
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

  function pcmsolver_new_v1112(input_reading, nr_nuclei, charges, coordinates, symmetry_info, &
       parsed_fname, host_input, writer) result(context) &
       bind(C)
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

  function pcmsolver_new_read_host(nr_nuclei, charges, coordinates, symmetry_info, writer) result(context) &
       bind(C)
    import
    integer(c_int), intent(in), value :: nr_nuclei
    real(c_double), intent(in)        :: charges(*)
    real(c_double), intent(in)        :: coordinates(*)
    integer(c_int), intent(in)        :: symmetry_info(*)
    type(c_funptr), intent(in), value :: writer
    type(c_ptr) :: context
  end function

  subroutine pcmsolver_set_bool_option_c(context, param, val) &
       bind(C, name="pcmsolver_set_bool_option")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: param(*)
    logical(c_bool), intent(in), value :: val
  end subroutine

  subroutine pcmsolver_set_int_option_c(context, param, val) &
       bind(C, name="pcmsolver_set_int_option")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: param(*)
    integer(c_int), intent(in), value :: val
  end subroutine

  subroutine pcmsolver_set_double_option_c(context, param, val) &
       bind(C, name="pcmsolver_set_double_option")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: param(*)
    real(c_double), intent(in), value :: val
  end subroutine

  subroutine pcmsolver_set_string_option_c(context, param, val) &
       bind(C, name="pcmsolver_set_string_option")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: param(*)
    character(kind=c_char, len=1), intent(in) :: val(*)
  end subroutine

  subroutine pcmsolver_refresh(context) &
       bind(C)
    import
    type(c_ptr), value :: context
  end subroutine

  subroutine pcmsolver_delete(context) &
       bind(C)
    import
    type(c_ptr), value :: context
  end subroutine

  function pcmsolver_is_compatible_library() result(compatible) &
       bind(C)
    import
    logical(c_bool) :: compatible
  end function

  subroutine pcmsolver_print(context) &
       bind(C)
    import
    type(c_ptr), value :: context
  end subroutine

  subroutine pcmsolver_citation(writer) &
       bind(C)
     import
     type(c_funptr), intent(in), value :: writer
   end subroutine

   function pcmsolver_get_cavity_size(context) result(nr_points) &
        bind(C)
    import
    type(c_ptr), value :: context
    integer(c_int)  :: nr_points
  end function

  function pcmsolver_get_irreducible_cavity_size(context) result(nr_points_irr) &
       bind(C)
    import
    type(c_ptr), value :: context
    integer(c_int)  :: nr_points_irr
  end function

  subroutine pcmsolver_get_centers(context, centers) &
       bind(C)
    import
    type(c_ptr), value :: context
    real(c_double), intent(inout) :: centers(*)
  end subroutine

  subroutine pcmsolver_get_center(context, its, center) &
       bind(C)
    import
    type(c_ptr), value :: context
    integer(c_int), value, intent(in) :: its
    real(c_double), intent(inout) :: center(*)
  end subroutine

  subroutine pcmsolver_get_areas(context, areas) &
       bind(C)
    import
    type(c_ptr), value :: context
    real(c_double), intent(inout) :: areas(*)
  end subroutine

  subroutine pcmsolver_compute_asc_c(context, mep_name, asc_name, irrep) &
       bind(C, name="pcmsolver_compute_asc")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: mep_name(*), asc_name(*)
    integer(c_int), value, intent(in) :: irrep
  end subroutine

  subroutine pcmsolver_compute_response_asc_c(context, mep_name, asc_name, irrep) &
       bind(C, name="pcmsolver_compute_response_asc")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: mep_name(*), asc_name(*)
    integer(c_int), value, intent(in) :: irrep
  end subroutine

  function pcmsolver_compute_polarization_energy_c(context, mep_name, asc_name) result(energy) &
       bind(C, name="pcmsolver_compute_polarization_energy")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: mep_name(*), asc_name(*)
    real(c_double) :: energy
  end function

  function pcmsolver_get_asc_dipole_c(context, asc_name, dipole) result(mu) &
       bind(C, name="pcmsolver_get_asc_dipole")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: asc_name(*)
    real(c_double), intent(inout) :: dipole(*)
    real(c_double) :: mu
  end function

  subroutine pcmsolver_get_surface_function_c(context, f_size, values, name) &
       bind(C, name="pcmsolver_get_surface_function")
    import
    type(c_ptr), value :: context
    integer(c_int), value, intent(in) :: f_size
    real(c_double), intent(inout) :: values(*)
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine

  subroutine pcmsolver_set_surface_function_c(context, f_size, values, name) &
       bind(C, name="pcmsolver_set_surface_function")
    import
    type(c_ptr), value :: context
    integer(c_int), value, intent(in) :: f_size
    real(c_double), intent(in) :: values(*)
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine

  subroutine pcmsolver_print_surface_function_c(context, name) &
       bind(C, name="pcmsolver_print_surface_function")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine

  subroutine pcmsolver_save_surface_functions(context) &
       bind(C)
    import
    type(c_ptr), value :: context
  end subroutine

  subroutine pcmsolver_save_surface_function_c(context, name) &
       bind(C, name="pcmsolver_save_surface_function")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine

  subroutine pcmsolver_load_surface_function_c(context, name) &
       bind(C, name="pcmsolver_load_surface_function")
    import
    type(c_ptr), value :: context
    character(kind=c_char, len=1), intent(in) :: name(*)
  end subroutine

  subroutine pcmsolver_write_timings(context) &
       bind(C)
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

  subroutine pcmsolver_compute_asc(context, mep_name, asc_name, irrep)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: mep_name, asc_name
    integer(c_int), intent(in) :: irrep

    call pcmsolver_compute_asc_c(context, pcmsolver_fstring_to_carray(mep_name), &
                                          pcmsolver_fstring_to_carray(asc_name), irrep)
  end subroutine

  subroutine pcmsolver_compute_response_asc(context, mep_name, asc_name, irrep)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: mep_name, asc_name
    integer(c_int), intent(in) :: irrep

    call pcmsolver_compute_response_asc_c(context, pcmsolver_fstring_to_carray(mep_name), &
         pcmsolver_fstring_to_carray(asc_name), irrep)
  end subroutine

  function pcmsolver_compute_polarization_energy(context, mep_name, asc_name) result(energy)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: mep_name, asc_name
    real(c_double) :: energy

    energy =  pcmsolver_compute_polarization_energy_c(context, &
                                                      pcmsolver_fstring_to_carray(mep_name), &
                                                      pcmsolver_fstring_to_carray(asc_name))
  end function

  function pcmsolver_get_asc_dipole(context, asc_name, dipole) result(mu)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: asc_name
    real(c_double), intent(inout) :: dipole(*)
    real(c_double) :: mu

    mu = pcmsolver_get_asc_dipole_c(context, pcmsolver_fstring_to_carray(asc_name), dipole)
  end function

  subroutine pcmsolver_set_surface_function(context, f_size, values, name)
    type(c_ptr), value :: context
    integer(c_int), intent(in) :: f_size
    real(c_double), intent(inout) :: values(*)
    character(kind=c_char, len=*), intent(in) :: name

    call pcmsolver_set_surface_function_c(context, f_size, values, pcmsolver_fstring_to_carray(name))
  end subroutine

  subroutine pcmsolver_get_surface_function(context, f_size, values, name)
    type(c_ptr), value :: context
    integer(c_int), intent(in) :: f_size
    real(c_double), intent(inout) :: values(*)
    character(kind=c_char, len=*), intent(in) :: name

    call pcmsolver_get_surface_function_c(context, f_size, values, pcmsolver_fstring_to_carray(name))
  end subroutine

  subroutine pcmsolver_print_surface_function(context, name)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: name

    call pcmsolver_print_surface_function_c(context, pcmsolver_fstring_to_carray(name))
  end subroutine

  subroutine pcmsolver_save_surface_function(context, name)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: name

    call pcmsolver_save_surface_function_c(context, pcmsolver_fstring_to_carray(name))
  end subroutine

  subroutine pcmsolver_load_surface_function(context, name)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: name

    call pcmsolver_load_surface_function_c(context, pcmsolver_fstring_to_carray(name))
  end subroutine

  subroutine pcmsolver_set_bool_option(context, param, val)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: param
    logical(c_bool), intent(in) :: val

    call pcmsolver_set_bool_option_c(context, pcmsolver_fstring_to_carray(param), val)
  end subroutine

  subroutine pcmsolver_set_int_option(context, param, val)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: param
    integer(c_int), intent(in) :: val

    call pcmsolver_set_int_option_c(context, pcmsolver_fstring_to_carray(param), val)
  end subroutine

  subroutine pcmsolver_set_double_option(context, param, val)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: param
    real(c_double), intent(in) :: val

    call pcmsolver_set_double_option_c(context, pcmsolver_fstring_to_carray(param), val)
  end subroutine

  subroutine pcmsolver_set_string_option(context, param, val)
    type(c_ptr), value :: context
    character(kind=c_char, len=*), intent(in) :: param
    character(kind=c_char, len=*), intent(in) :: val

    call pcmsolver_set_string_option_c(context, pcmsolver_fstring_to_carray(param), &
                                                pcmsolver_fstring_to_carray(val))
  end subroutine
end module
