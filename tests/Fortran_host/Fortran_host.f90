!
! PCMSolver, an API for the Polarizable Continuum Model
! Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
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

program pcm_fortran_host

      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
      use pcmsolver
      use utilities
      use testing

      implicit none

      type(c_ptr) :: pcm_context
      integer(c_int) :: nr_nuclei
      real(c_double), allocatable :: charges(:)
      real(c_double), allocatable :: coordinates(:)
      integer(c_int) :: symmetry_info(4)
      type(PCMInput) :: host_input
      logical :: log_open, log_exist
      ! Shows two different, but equivalent ways of defining labels for surface functions
      character(kind=c_char, len=1) :: mep_lbl(7) = (/'N', 'u', 'c', 'M', 'E', 'P', c_null_char/)
      character(kind=c_char, len=1) :: asc_lbl(7) = (/'N', 'u', 'c', 'A', 'S', 'C', c_null_char/)
      character(kind=c_char, len=1) :: asc_B3g_lbl(7) = (/'O', 'I', 'T', 'A', 'S', 'C', c_null_char/)
      character(kind=c_char, len=1) :: asc_neq_B3g_lbl(7) = (/'A', 'S', 'C', 'B', '3', 'g', c_null_char/)
      real(c_double), allocatable :: grid(:), mep(:), asc_Ag(:), asc_B3g(:), asc_neq_B3g(:), areas(:)
      integer(c_int) :: irrep
      integer(c_int) :: grid_size, irr_grid_size
      real(c_double) :: energy
      ! Reference values for scalar quantities
      integer(c_int), parameter :: ref_size = 576, ref_irr_size = 72
      real(c_double), parameter :: ref_energy = -0.437960027982

      if (.not. pcmsolver_is_compatible_library()) then
        write(error_unit, *) 'PCMSolver library not compatible!'
        stop
      end if

      ! Open a file for the output...
      inquire(file = 'Fortran_host.out', opened = log_open, &
        exist = log_exist)
      if (log_exist) then
        open(unit = output_unit,                   &
          file = 'Fortran_host.out',       &
          status = 'unknown',       &
          form = 'formatted',       &
          access = 'sequential')
        close(unit = output_unit, status = 'delete')
      end if
      open(unit = output_unit,                      &
        file = 'Fortran_host.out',          &
        status = 'new',              &
        form = 'formatted',          &
        access = 'sequential')
      rewind(output_unit)
      write(output_unit, *) 'Starting a PCMSolver calculation'

      nr_nuclei = 6_c_int
      allocate(charges(nr_nuclei))
      allocate(coordinates(3*nr_nuclei))

      ! Use C2H4 in D2h symmetry
      charges = (/ 6.0_c_double, 1.0_c_double, 1.0_c_double, &
                   6.0_c_double, 1.0_c_double, 1.0_c_double /)
      coordinates = (/ 0.0_c_double,  0.0_c_double,       1.257892_c_double, &
                      0.0_c_double,  1.745462_c_double,  2.342716_c_double, &
                      0.0_c_double, -1.745462_c_double,  2.342716_c_double, &
                      0.0_c_double,  0.0_c_double,      -1.257892_c_double, &
                      0.0_c_double,  1.745462_c_double, -2.342716_c_double, &
                      0.0_c_double, -1.745462_c_double, -2.342716_c_double /)

      ! This means the molecular point group has three generators:
      ! the Oxy, Oxz and Oyz planes
      symmetry_info = (/ 3, 4, 2, 1 /)

      host_input = pcmsolver_input()

      pcm_context = pcmsolver_new(PCMSOLVER_READER_HOST,           &
                                  nr_nuclei, charges, coordinates, &
                                  symmetry_info, host_input,       &
                                  c_funloc(host_writer))

      call pcmsolver_print(pcm_context)

      grid_size = pcmsolver_get_cavity_size(pcm_context)
      irr_grid_size = pcmsolver_get_irreducible_cavity_size(pcm_context)
      allocate(grid(3*grid_size))
      grid = 0.0_c_double
      call pcmsolver_get_centers(pcm_context, grid)
      allocate(areas(grid_size))
      call pcmsolver_get_areas(pcm_context, areas)

      allocate(mep(grid_size))
      mep = 0.0_c_double
      mep = nuclear_mep(nr_nuclei, charges, reshape(coordinates, (/ 3_c_int, nr_nuclei /)), &
                        grid_size, reshape(grid, (/ 3_c_int, grid_size /)))
      call pcmsolver_set_surface_function(pcm_context, grid_size, mep, mep_lbl)
      ! This is the Ag irreducible representation (totally symmetric)
      irrep = 0_c_int
      call pcmsolver_compute_asc(pcm_context, mep_lbl, asc_lbl, irrep)
      allocate(asc_Ag(grid_size))
      asc_Ag = 0.0_c_double
      call pcmsolver_get_surface_function(pcm_context, grid_size, asc_Ag, asc_lbl)

      energy = pcmsolver_compute_polarization_energy(pcm_context, mep_lbl, asc_lbl)

      write(output_unit, '(A, F20.12)') 'Polarization energy = ', energy

      allocate(asc_neq_B3g(grid_size))
      asc_neq_B3g = 0.0_c_double
      ! This is the B3g irreducible representation
      irrep = 3_c_int
      call pcmsolver_compute_response_asc(pcm_context, mep_lbl, asc_neq_B3g_lbl, irrep)
      call pcmsolver_get_surface_function(pcm_context, grid_size, asc_neq_B3g, asc_neq_B3g_lbl)

      ! Equilibrium ASC in B3g symmetry.
      ! This is an internal check: the relevant segment of the vector
      ! should be the same as the one calculated using pcmsolver_compute_response_asc
      allocate(asc_B3g(grid_size))
      asc_B3g = 0.0_c_double
      ! This is the B3g irreducible representation
      irrep = 3_c_int
      call pcmsolver_compute_asc(pcm_context, mep_lbl, asc_B3g_lbl, irrep)
      call pcmsolver_get_surface_function(pcm_context, grid_size, asc_B3g, asc_B3g_lbl)

      ! Check that everything calculated is OK
      ! Cavity size
      if (grid_size .ne. ref_size) then
        write(error_unit, *) 'Error in the cavity size, please file an issue on: https://github.com/PCMSolver/pcmsolver'
        stop
      else
        write(output_unit, *) 'Test on cavity size: PASSED'
      end if
      ! Irreducible cavity size
      if (irr_grid_size .ne. ref_irr_size) then
        write(error_unit, *) 'Error in the irreducible cavity size, please file an issue on: https://github.com/PCMSolver/pcmsolver'
        stop
      else
        write(output_unit, *) 'Test on irreducible cavity size: PASSED'
      end if
      ! Polarization energy
      if (.not. check_unsigned_error(energy, ref_energy, 1.0e-7_c_double)) then
        write(error_unit, *) 'Error in the polarization energy, please file an issue on: https://github.com/PCMSolver/pcmsolver'
        stop
      else
        write(output_unit, *) 'Test on polarization energy: PASSED'
      end if
      ! Surface functions
      call test_surface_functions(grid_size, mep, asc_Ag, asc_B3g, asc_neq_B3g, areas)

      call pcmsolver_save_surface_function(pcm_context, mep_lbl)
      call pcmsolver_load_surface_function(pcm_context, mep_lbl)

      call pcmsolver_write_timings(pcm_context)

      call pcmsolver_delete(pcm_context)

      deallocate(charges)
      deallocate(coordinates)
      deallocate(grid)
      deallocate(mep)
      deallocate(asc_Ag)
      deallocate(asc_B3g)
      deallocate(asc_neq_B3g)
      deallocate(areas)

      close(output_unit)

end program pcm_fortran_host
