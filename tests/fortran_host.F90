module utilities

  implicit none

  public pcmsolver_input
  public nuclear_mep
  public check_unsigned_error

  private

contains

  ! Performs syntactic checks on PCMSolver input and fills the data structures
  ! holding input data
  function pcmsolver_input() result(host_input)

    use, intrinsic :: iso_c_binding
    use pcmsolver, only: PCMInput

    type(PCMInput) :: host_input

    ! These parameters would be set by the host input reading subroutine(s)
    ! Notice that the strings have a maximum pre-set length.
    ! This is to ensure C interoperability
    ! Length and area parameters are all assumed to be in Angstrom,
    ! the module will convert to Bohr internally
    character(kind=c_char, len=6) :: pcmmod_cavity_type = 'gepol'//c_null_char
    integer :: pcmmod_patch_level = 2
    real(c_double) :: pcmmod_coarsity = 0.5
    real(c_double) :: pcmmod_cavity_area = 0.2
    real(c_double) :: pcmmod_min_distance = 0.1
    integer :: pcmmod_der_order = 4
    logical(c_bool) :: pcmmod_scaling = .true.
    character(kind=c_char, len=6) :: pcmmod_radii_set = 'bondi'//c_null_char
    character(kind=c_char, len=11) :: pcmmod_restart_name = 'cavity.npz'//c_null_char
    real(c_double) :: pcmmod_min_radius = 100.0
    character(kind=c_char, len=7) :: pcmmod_solver_type = 'iefpcm'//c_null_char
    character(kind=c_char, len=6) :: pcmmod_solvent = 'water'//c_null_char
    character(kind=c_char, len=11) :: pcmmod_equation_type = 'secondkind'//c_null_char
    real(c_double) :: pcmmod_correction = 0.0
    real(c_double) :: pcmmod_probe_radius = 1.0
    character(kind=c_char, len=7) :: pcmmod_inside_type = 'vacuum'//c_null_char
    character(kind=c_char, len=18) :: pcmmod_outside_type = 'uniformdielectric'//c_null_char
    real(c_double) :: pcmmod_outside_epsilon = 1.0

    character(kind=c_char, len=1) :: cavity_type(8)
    character(kind=c_char, len=1) :: radii_set(8)
    character(kind=c_char, len=1) :: restart_name(20)
    character(kind=c_char, len=1) :: solver_type(7)
    character(kind=c_char, len=1) :: solvent(16)
    character(kind=c_char, len=1) :: equation_type(11)
    character(kind=c_char, len=1) :: inside_type(7)
    character(kind=c_char, len=1) :: outside_type(22)

    call pcmsolver_f2c_string(pcmmod_cavity_type, cavity_type, 6_c_int)
    host_input%cavity_type  = cavity_type
    host_input%patch_level  = int(pcmmod_patch_level, kind=c_int)
    host_input%coarsity     = pcmmod_coarsity
    host_input%area         = pcmmod_cavity_area
    host_input%min_distance = pcmmod_min_distance
    host_input%der_order    = int(pcmmod_der_order, kind=c_int)
    host_input%scaling      = pcmmod_scaling
    call pcmsolver_f2c_string(pcmmod_radii_set, radii_set, 6_c_int)
    host_input%radii_set    = radii_set
    call pcmsolver_f2c_string(pcmmod_restart_name, restart_name, 11_c_int)
    host_input%restart_name = restart_name
    host_input%min_radius   = pcmmod_min_radius

    call pcmsolver_f2c_string(pcmmod_solver_type, solver_type, 7_c_int)
    host_input%solver_type   = solver_type
    call pcmsolver_f2c_string(pcmmod_solvent, solvent, 6_c_int)
    host_input%solvent       = solvent
    call pcmsolver_f2c_string(pcmmod_equation_type, equation_type, 11_c_int)
    host_input%equation_type = equation_type
    host_input%correction    = pcmmod_correction
    host_input%probe_radius  = pcmmod_probe_radius

    call pcmsolver_f2c_string(pcmmod_inside_type, inside_type, 7_c_int)
    host_input%inside_type     = inside_type
    host_input%outside_epsilon = pcmmod_outside_epsilon
    call pcmsolver_f2c_string(pcmmod_outside_type, outside_type, 18_c_int)
    host_input%outside_type    = outside_type

  end function pcmsolver_input

  !> \brief calculates nuclear molecular electrostatic potential (MEP) at cavity points
  !> \author Roberto Di Remigio
  !> \date 2015
  !> \param nr_nuclei number of atomic centers
  !> \param charges atomic charges
  !> \param coordinates coordinates of the atomic centers
  !> \param grid_size the number of cavity points
  !> \param grid the cavity points
  !>
  !> Calculates the nuclear MEP given a list of points.
  pure function nuclear_mep(nr_nuclei, charges, coordinates, grid_size, grid) result(n_mep)

    use, intrinsic :: iso_c_binding
    integer(c_int), intent(in) :: nr_nuclei
    real(c_double), intent(in) :: charges(nr_nuclei)
    real(c_double), intent(in) :: coordinates(3, nr_nuclei)
    integer(c_size_t), intent(in) :: grid_size
    real(c_double), intent(in) :: grid(3, grid_size)
    real(c_double) :: n_mep(grid_size)

    real(c_double) :: dist
    integer(c_int) :: i
    integer(c_size_t) :: ipoint

    n_mep = 0.0_c_double
    LoopOnAtoms: do i = 1, nr_nuclei
      LoopOnPoints: do ipoint = 1, grid_size
        dist = (coordinates(1, i) - grid(1, ipoint))**2 + &
          (coordinates(2, i) - grid(2, ipoint))**2 + &
          (coordinates(3, i) - grid(3, ipoint))**2
        dist = sqrt(dist)
        n_mep(ipoint) = n_mep(ipoint) + (charges(i) / dist)
      end do LoopOnPoints
    end do LoopOnAtoms

  end function nuclear_mep

  !> \brief Compares calculated and reference values within a threshold
  !> \author Roberto Di Remigio
  !> \date 2015
  !> \param calculated calculated value
  !> \param reference  reference value
  !> \param threshold  comparison threshold
  !> \return true if the calculated and reference values are compatible, false otherwise
  elemental function check_unsigned_error(calculated, reference, threshold) result(ok)

    use, intrinsic :: iso_c_binding, only: c_double, c_bool
    real(c_double), intent(in) :: calculated, reference, threshold
    logical(c_bool) :: ok
    real(c_double) :: err

    err = calculated - reference
    ok = (abs(err) <= abs(calculated) * threshold)

  end function check_unsigned_error

end module utilities

program pcm_fortran_host

      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
      use pcmsolver
      use utilities
      
      implicit none
      
      type(c_ptr) :: pcm_context
      integer(c_int) :: nr_nuclei
      real(c_double), allocatable :: charges(:)
      real(c_double), allocatable :: coordinates(:)
      integer(c_int) :: symmetry_info(4)
      type(PCMInput) :: host_input
      logical :: log_open, log_exist
      character(kind=c_char, len=7) :: mep_lbl, asc_lbl, asc_B3g_lbl, asc_neq_B3g_lbl
      real(c_double), allocatable :: grid(:), mep(:), asc_Ag(:), asc_B3g(:), asc_neq_B3g(:)
      integer(c_int) :: irrep
      integer(c_size_t) :: grid_size, irr_grid_size
      real(c_double) :: energy
      ! Reference values for scalar quantities
      integer(c_size_t), parameter :: ref_size = 576, ref_irr_size = 72
      real(c_double), parameter :: ref_energy = -0.437960027982

      if (.not. pcmsolver_is_compatible_library()) then
        write(error_unit, *) 'PCMSolver library not compatible!'
        stop -1
      end if

      ! Open a file for the output...
      inquire(file = 'fortran_host.log', opened = log_open, &
        exist = log_exist)
      if (log_exist) then
        open(unit = output_unit,                   &
          file = 'fortran_host.log',       &
          status = 'unknown',       &
          form = 'formatted',       &
          access = 'sequential')
        close(unit = output_unit, status = 'delete')
      end if
      open(unit = output_unit,                      &
        file = 'fortran_host.log',          &
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
                                  symmetry_info, host_input)

      call pcmsolver_print(pcm_context)

      grid_size = pcmsolver_get_cavity_size(pcm_context)
      irr_grid_size = pcmsolver_get_irreducible_cavity_size(pcm_context)
      allocate(grid(3*grid_size))
      grid = 0.0_c_double
      call pcmsolver_get_centers(pcm_context, grid)

      allocate(mep(grid_size))
      mep = 0.0_c_double
      mep = nuclear_mep(nr_nuclei, charges, reshape(coordinates, (/ 3, nr_nuclei /)), &
                        grid_size, reshape(grid, (/ 3_c_size_t, grid_size /)))
      mep_lbl = 'NucMEP'//c_null_char
      call pcmsolver_set_surface_function(pcm_context, grid_size, mep, mep_lbl)
      asc_lbl = 'NucASC'//c_null_char
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
      asc_neq_B3g_lbl = 'OITASC'//c_null_char
      ! This is the B3g irreducible representation
      irrep = 3_c_int
      call pcmsolver_compute_response_asc(pcm_context, mep_lbl, asc_neq_B3g_lbl, irrep)
      call pcmsolver_get_surface_function(pcm_context, grid_size, asc_neq_B3g, asc_neq_B3g_lbl)

      ! Equilibrium ASC in B3g symmetry.
      ! This is an internal check: the relevant segment of the vector
      ! should be the same as the one calculated using pcmsolver_compute_response_asc
      allocate(asc_B3g(grid_size))
      asc_B3g = 0.0_c_double
      asc_B3g_lbl = 'ASCB3g'//c_null_char
      ! This is the B3g irreducible representation
      irrep = 3_c_size_t
      call pcmsolver_compute_asc(pcm_context, mep_lbl, asc_B3g_lbl, irrep)
      call pcmsolver_get_surface_function(pcm_context, grid_size, asc_B3g, asc_B3g_lbl)

      ! Check that everything calculated is OK
      ! Cavity size
      if (grid_size .ne. ref_size) then
        write(error_unit, *) 'Error in the cavity size, please file an issue on: https://github.com/PCMSolver/pcmsolver'
        stop -1
      else 
        write(output_unit, *) 'Test on cavity size: PASSED'
      end if
      ! Irreducible cavity size
      if (irr_grid_size .ne. ref_irr_size) then
        write(error_unit, *) 'Error in the irreducible cavity size, please file an issue on: https://github.com/PCMSolver/pcmsolver'
        stop -1
      else 
        write(output_unit, *) 'Test on irreducible cavity size: PASSED'
      end if
      ! Polarization energy
      if (.not. check_unsigned_error(energy, ref_energy, 1.0e-7_c_double)) then
        write(error_unit, *) 'Error in the polarization energy, please file an issue on: https://github.com/PCMSolver/pcmsolver'
        stop -1
      else 
        write(output_unit, *) 'Test on polarization energy: PASSED'
      end if
      ! Surface functions
      call test_surface_functions(grid_size, mep, asc_Ag, asc_B3g, asc_neq_B3g)

      call pcmsolver_write_timings(pcm_context)

      call pcmsolver_delete(pcm_context)

      deallocate(charges)
      deallocate(coordinates)
      deallocate(grid)
      deallocate(mep)
      deallocate(asc_Ag)
      deallocate(asc_B3g)
      deallocate(asc_neq_B3g)

      close(output_unit)

end program pcm_fortran_host

subroutine host_writer(message, message_length) bind(c, name='host_writer')

  use, intrinsic :: iso_c_binding, only: c_char, c_size_t
  use, intrinsic :: iso_fortran_env, only: output_unit

  character(kind=c_char) :: message(*)
  integer(c_size_t), intent(in), value :: message_length

  character(len=message_length) :: f_message

  call pcmsolver_c2f_string(message, f_message, message_length)
  write(output_unit, '(1000A)') f_message

end subroutine host_writer

subroutine test_surface_functions(grid_size, mep, asc_Ag, asc_B3g, asc_neq_B3g)

  use, intrinsic :: iso_c_binding, only: c_bool, c_double, c_size_t
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use utilities, only: check_unsigned_error

  interface mep_reference
    pure function mep_reference() result(mep)
      use, intrinsic :: iso_c_binding, only: c_double, c_size_t
      real(c_double) :: mep(576)
    end function mep_reference
  end interface mep_reference
  
  interface asc_Ag_reference
    pure function asc_Ag_reference() result(asc_Ag)
      use, intrinsic :: iso_c_binding, only: c_double, c_size_t
      real(c_double) :: asc_Ag(576)
    end function asc_Ag_reference
  end interface asc_Ag_reference
  
  interface asc_B3g_reference
    pure function asc_B3g_reference() result(asc_B3g)
      use, intrinsic :: iso_c_binding, only: c_double, c_size_t
      real(c_double) :: asc_B3g(576)
    end function asc_B3g_reference
  end interface asc_B3g_reference
  
  interface asc_neq_B3g_reference
    pure function asc_neq_B3g_reference() result(asc_neq_B3g)
      use, intrinsic :: iso_c_binding, only: c_double, c_size_t
      real(c_double) :: asc_neq_B3g(576)
    end function asc_neq_B3g_reference
  end interface asc_neq_B3g_reference

  integer(c_size_t), intent(in) :: grid_size
  real(c_double), intent(in) :: mep(grid_size), asc_Ag(grid_size)
  real(c_double), intent(in) :: asc_B3g(grid_size), asc_neq_B3g(grid_size)
  logical(c_bool) :: check(grid_size)
  integer(c_size_t) :: ipoint

  check = check_unsigned_error(mep, mep_reference(), 1.0e-07_c_double)
  do ipoint = 1, grid_size
    if (.not. check(ipoint)) then
      write(error_unit, *) 'Error in MEP, please file an issue on: https://github.com/PCMSolver/pcmsolver'
      stop -1
    end if
  end do
  write(output_unit, *) 'Test on MEP: PASSED'

  check = check_unsigned_error(asc_Ag, asc_Ag_reference(), 1.0e-07_c_double)
  do ipoint = 1, grid_size
    if (.not. check(ipoint)) then
      write(error_unit, *) 'Error in ASC Ag, please file an issue on: https://github.com/PCMSolver/pcmsolver'
      stop -1
    end if
  end do
  write(output_unit, *) 'Test on ASC in Ag symmetry: PASSED'

  check = check_unsigned_error(asc_B3g, asc_B3g_reference(), 1.0e-07_c_double)
  do ipoint = 1, grid_size
    if (.not. check(ipoint)) then
       write(error_unit, *) 'Error in ASC B3g, please file an issue on: https://github.com/PCMSolver/pcmsolver'
       stop -1
    end if
  end do
  write(output_unit, *) 'Test on ASC in B3g symmetry: PASSED'

  check = check_unsigned_error(asc_neq_B3g, asc_neq_B3g_reference(), 1.0e-07_c_double)
  do ipoint = 1, grid_size
    if (.not. check(ipoint)) then
      write(error_unit, *) 'Error in nonequilibrium ASC B3g, please file an issue on: https://github.com/PCMSolver/pcmsolver'
      stop -1
    end if
  end do
  write(output_unit, *) 'Test on nonequilibrium ASC in B3g symmetry: PASSED'

end subroutine test_surface_functions

pure function mep_reference() result(mep)

  use, intrinsic :: iso_c_binding, only: c_double, c_size_t
  real(c_double) :: mep(576)

  mep = (/ 3.372240447717, &
    3.395095348797, &
    3.375472982012, &
    3.407830928930, &
    3.389673049936, &
    3.444734322778, &
    3.434450634681, &
    3.567420410364, &
    3.519403897623, &
    3.521500360692, &
    3.737204637428, &
    3.665836385359, &
    3.638475909577, &
    3.634697522671, &
    3.378728206062, &
    3.400511009364, &
    3.387155190391, &
    3.430448505960, &
    3.414568660060, &
    3.502974792852, &
    3.485472605945, &
    3.646792406411, &
    3.600687428544, &
    3.596042090542, &
    4.043143323646, &
    4.042627673887, &
    4.000207071572, &
    3.943830109299, &
    3.901090644538, &
    3.971528385050, &
    3.922106005834, &
    3.832443093990, &
    3.781559803821, &
    3.745079273102, &
    3.726956468823, &
    4.071479016092, &
    4.039193036616, &
    3.982347524785, &
    3.939905876613, &
    3.970318175141, &
    3.926366695705, &
    3.883604616653, &
    3.829224917150, &
    3.786789784941, &
    3.137381224862, &
    3.001405458723, &
    3.186612154619, &
    3.000297605619, &
    3.091578048072, &
    3.363118040231, &
    3.133206449617, &
    3.182855203960, &
    3.381998093999, &
    3.087995448210, &
    2.999536488978, &
    3.192470183633, &
    3.085172887610, &
    3.190455940555, &
    3.445493468024, &
    3.663096576477, &
    3.766481150018, &
    3.271760664144, &
    3.368194930794, &
    3.540250609025, &
    3.824911768840, &
    3.397511985055, &
    3.568518356734, &
    3.268466398614, &
    3.329351334592, &
    3.421372568835, &
    3.333950027352, &
    3.407773545920, &
    3.372240447717, &
    3.395095348797, &
    3.375472982012, &
    3.407830928930, &
    3.389673049936, &
    3.444734322778, &
    3.434450634681, &
    3.567420410364, &
    3.519403897623, &
    3.521500360692, &
    3.737204637428, &
    3.665836385359, &
    3.638475909577, &
    3.634697522671, &
    3.378728206062, &
    3.400511009364, &
    3.387155190391, &
    3.430448505960, &
    3.414568660060, &
    3.502974792852, &
    3.485472605945, &
    3.646792406411, &
    3.600687428544, &
    3.596042090542, &
    4.043143323646, &
    4.042627673887, &
    4.000207071572, &
    3.943830109299, &
    3.901090644538, &
    3.971528385050, &
    3.922106005834, &
    3.832443093990, &
    3.781559803821, &
    3.745079273102, &
    3.726956468823, &
    4.071479016092, &
    4.039193036616, &
    3.982347524785, &
    3.939905876613, &
    3.970318175141, &
    3.926366695705, &
    3.883604616653, &
    3.829224917150, &
    3.786789784941, &
    3.137381224862, &
    3.001405458723, &
    3.186612154619, &
    3.000297605619, &
    3.091578048072, &
    3.363118040231, &
    3.133206449617, &
    3.182855203960, &
    3.381998093999, &
    3.087995448210, &
    2.999536488978, &
    3.192470183633, &
    3.085172887610, &
    3.190455940555, &
    3.445493468024, &
    3.663096576477, &
    3.766481150018, &
    3.271760664144, &
    3.368194930794, &
    3.540250609025, &
    3.824911768840, &
    3.397511985055, &
    3.568518356734, &
    3.268466398614, &
    3.329351334592, &
    3.421372568835, &
    3.333950027352, &
    3.407773545920, &
    3.372240447717, &
    3.395095348797, &
    3.375472982012, &
    3.407830928930, &
    3.389673049936, &
    3.444734322778, &
    3.434450634681, &
    3.567420410364, &
    3.519403897623, &
    3.521500360692, &
    3.737204637428, &
    3.665836385359, &
    3.638475909577, &
    3.634697522671, &
    3.378728206062, &
    3.400511009364, &
    3.387155190391, &
    3.430448505960, &
    3.414568660060, &
    3.502974792852, &
    3.485472605945, &
    3.646792406411, &
    3.600687428544, &
    3.596042090542, &
    4.043143323646, &
    4.042627673887, &
    4.000207071572, &
    3.943830109299, &
    3.901090644538, &
    3.971528385050, &
    3.922106005834, &
    3.832443093990, &
    3.781559803821, &
    3.745079273102, &
    3.726956468823, &
    4.071479016092, &
    4.039193036616, &
    3.982347524785, &
    3.939905876613, &
    3.970318175141, &
    3.926366695705, &
    3.883604616653, &
    3.829224917150, &
    3.786789784941, &
    3.137381224862, &
    3.001405458723, &
    3.186612154619, &
    3.000297605619, &
    3.091578048072, &
    3.363118040231, &
    3.133206449617, &
    3.182855203960, &
    3.381998093999, &
    3.087995448210, &
    2.999536488978, &
    3.192470183633, &
    3.085172887610, &
    3.190455940555, &
    3.445493468024, &
    3.663096576477, &
    3.766481150018, &
    3.271760664144, &
    3.368194930794, &
    3.540250609025, &
    3.824911768840, &
    3.397511985055, &
    3.568518356734, &
    3.268466398614, &
    3.329351334592, &
    3.421372568835, &
    3.333950027352, &
    3.407773545920, &
    3.372240447717, &
    3.395095348797, &
    3.375472982012, &
    3.407830928930, &
    3.389673049936, &
    3.444734322778, &
    3.434450634681, &
    3.567420410364, &
    3.519403897623, &
    3.521500360692, &
    3.737204637428, &
    3.665836385359, &
    3.638475909577, &
    3.634697522671, &
    3.378728206062, &
    3.400511009364, &
    3.387155190391, &
    3.430448505960, &
    3.414568660060, &
    3.502974792852, &
    3.485472605945, &
    3.646792406411, &
    3.600687428544, &
    3.596042090542, &
    4.043143323646, &
    4.042627673887, &
    4.000207071572, &
    3.943830109299, &
    3.901090644538, &
    3.971528385050, &
    3.922106005834, &
    3.832443093990, &
    3.781559803821, &
    3.745079273102, &
    3.726956468823, &
    4.071479016092, &
    4.039193036616, &
    3.982347524785, &
    3.939905876613, &
    3.970318175141, &
    3.926366695705, &
    3.883604616653, &
    3.829224917150, &
    3.786789784941, &
    3.137381224862, &
    3.001405458723, &
    3.186612154619, &
    3.000297605619, &
    3.091578048072, &
    3.363118040231, &
    3.133206449617, &
    3.182855203960, &
    3.381998093999, &
    3.087995448210, &
    2.999536488978, &
    3.192470183633, &
    3.085172887610, &
    3.190455940555, &
    3.445493468024, &
    3.663096576477, &
    3.766481150018, &
    3.271760664144, &
    3.368194930794, &
    3.540250609025, &
    3.824911768840, &
    3.397511985055, &
    3.568518356734, &
    3.268466398614, &
    3.329351334592, &
    3.421372568835, &
    3.333950027352, &
    3.407773545920, &
    3.372240447717, &
    3.395095348797, &
    3.375472982012, &
    3.407830928930, &
    3.389673049936, &
    3.444734322778, &
    3.434450634681, &
    3.567420410364, &
    3.519403897623, &
    3.521500360692, &
    3.737204637428, &
    3.665836385359, &
    3.638475909577, &
    3.634697522671, &
    3.378728206062, &
    3.400511009364, &
    3.387155190391, &
    3.430448505960, &
    3.414568660060, &
    3.502974792852, &
    3.485472605945, &
    3.646792406411, &
    3.600687428544, &
    3.596042090542, &
    4.043143323646, &
    4.042627673887, &
    4.000207071572, &
    3.943830109299, &
    3.901090644538, &
    3.971528385050, &
    3.922106005834, &
    3.832443093990, &
    3.781559803821, &
    3.745079273102, &
    3.726956468823, &
    4.071479016092, &
    4.039193036616, &
    3.982347524785, &
    3.939905876613, &
    3.970318175141, &
    3.926366695705, &
    3.883604616653, &
    3.829224917150, &
    3.786789784941, &
    3.137381224862, &
    3.001405458723, &
    3.186612154619, &
    3.000297605619, &
    3.091578048072, &
    3.363118040231, &
    3.133206449617, &
    3.182855203960, &
    3.381998093999, &
    3.087995448210, &
    2.999536488978, &
    3.192470183633, &
    3.085172887610, &
    3.190455940555, &
    3.445493468024, &
    3.663096576477, &
    3.766481150018, &
    3.271760664144, &
    3.368194930794, &
    3.540250609025, &
    3.824911768840, &
    3.397511985055, &
    3.568518356734, &
    3.268466398614, &
    3.329351334592, &
    3.421372568835, &
    3.333950027352, &
    3.407773545920, &
    3.372240447717, &
    3.395095348797, &
    3.375472982012, &
    3.407830928930, &
    3.389673049936, &
    3.444734322778, &
    3.434450634681, &
    3.567420410364, &
    3.519403897623, &
    3.521500360692, &
    3.737204637428, &
    3.665836385359, &
    3.638475909577, &
    3.634697522671, &
    3.378728206062, &
    3.400511009364, &
    3.387155190391, &
    3.430448505960, &
    3.414568660060, &
    3.502974792852, &
    3.485472605945, &
    3.646792406411, &
    3.600687428544, &
    3.596042090542, &
    4.043143323646, &
    4.042627673887, &
    4.000207071572, &
    3.943830109299, &
    3.901090644538, &
    3.971528385050, &
    3.922106005834, &
    3.832443093990, &
    3.781559803821, &
    3.745079273102, &
    3.726956468823, &
    4.071479016092, &
    4.039193036616, &
    3.982347524785, &
    3.939905876613, &
    3.970318175141, &
    3.926366695705, &
    3.883604616653, &
    3.829224917150, &
    3.786789784941, &
    3.137381224862, &
    3.001405458723, &
    3.186612154619, &
    3.000297605619, &
    3.091578048072, &
    3.363118040231, &
    3.133206449617, &
    3.182855203960, &
    3.381998093999, &
    3.087995448210, &
    2.999536488978, &
    3.192470183633, &
    3.085172887610, &
    3.190455940555, &
    3.445493468024, &
    3.663096576477, &
    3.766481150018, &
    3.271760664144, &
    3.368194930794, &
    3.540250609025, &
    3.824911768840, &
    3.397511985055, &
    3.568518356734, &
    3.268466398614, &
    3.329351334592, &
    3.421372568835, &
    3.333950027352, &
    3.407773545920, &
    3.372240447717, &
    3.395095348797, &
    3.375472982012, &
    3.407830928930, &
    3.389673049936, &
    3.444734322778, &
    3.434450634681, &
    3.567420410364, &
    3.519403897623, &
    3.521500360692, &
    3.737204637428, &
    3.665836385359, &
    3.638475909577, &
    3.634697522671, &
    3.378728206062, &
    3.400511009364, &
    3.387155190391, &
    3.430448505960, &
    3.414568660060, &
    3.502974792852, &
    3.485472605945, &
    3.646792406411, &
    3.600687428544, &
    3.596042090542, &
    4.043143323646, &
    4.042627673887, &
    4.000207071572, &
    3.943830109299, &
    3.901090644538, &
    3.971528385050, &
    3.922106005834, &
    3.832443093990, &
    3.781559803821, &
    3.745079273102, &
    3.726956468823, &
    4.071479016092, &
    4.039193036616, &
    3.982347524785, &
    3.939905876613, &
    3.970318175141, &
    3.926366695705, &
    3.883604616653, &
    3.829224917150, &
    3.786789784941, &
    3.137381224862, &
    3.001405458723, &
    3.186612154619, &
    3.000297605619, &
    3.091578048072, &
    3.363118040231, &
    3.133206449617, &
    3.182855203960, &
    3.381998093999, &
    3.087995448210, &
    2.999536488978, &
    3.192470183633, &
    3.085172887610, &
    3.190455940555, &
    3.445493468024, &
    3.663096576477, &
    3.766481150018, &
    3.271760664144, &
    3.368194930794, &
    3.540250609025, &
    3.824911768840, &
    3.397511985055, &
    3.568518356734, &
    3.268466398614, &
    3.329351334592, &
    3.421372568835, &
    3.333950027352, &
    3.407773545920, &
    3.372240447717, &
    3.395095348797, &
    3.375472982012, &
    3.407830928930, &
    3.389673049936, &
    3.444734322778, &
    3.434450634681, &
    3.567420410364, &
    3.519403897623, &
    3.521500360692, &
    3.737204637428, &
    3.665836385359, &
    3.638475909577, &
    3.634697522671, &
    3.378728206062, &
    3.400511009364, &
    3.387155190391, &
    3.430448505960, &
    3.414568660060, &
    3.502974792852, &
    3.485472605945, &
    3.646792406411, &
    3.600687428544, &
    3.596042090542, &
    4.043143323646, &
    4.042627673887, &
    4.000207071572, &
    3.943830109299, &
    3.901090644538, &
    3.971528385050, &
    3.922106005834, &
    3.832443093990, &
    3.781559803821, &
    3.745079273102, &
    3.726956468823, &
    4.071479016092, &
    4.039193036616, &
    3.982347524785, &
    3.939905876613, &
    3.970318175141, &
    3.926366695705, &
    3.883604616653, &
    3.829224917150, &
    3.786789784941, &
    3.137381224862, &
    3.001405458723, &
    3.186612154619, &
    3.000297605619, &
    3.091578048072, &
    3.363118040231, &
    3.133206449617, &
    3.182855203960, &
    3.381998093999, &
    3.087995448210, &
    2.999536488978, &
    3.192470183633, &
    3.085172887610, &
    3.190455940555, &
    3.445493468024, &
    3.663096576477, &
    3.766481150018, &
    3.271760664144, &
    3.368194930794, &
    3.540250609025, &
    3.824911768840, &
    3.397511985055, &
    3.568518356734, &
    3.268466398614, &
    3.329351334592, &
    3.421372568835, &
    3.333950027352, &
    3.407773545920 /)

end function mep_reference

pure function asc_Ag_reference() result(asc_Ag)

  use, intrinsic :: iso_c_binding, only: c_double, c_size_t
  real(c_double) :: asc_Ag(576)

  asc_Ag = (/ -0.002085687695, &
    -0.002038528244, &
    -0.003385497599, &
    -0.002692240783, &
    -0.004538826296, &
    -0.006256935957, &
    -0.004427462623, &
    -0.003425534363, &
    -0.005219445681, &
    -0.003263010379, &
    -0.003395946702, &
    -0.005008248255, &
    -0.003558712718, &
    -0.002122854689, &
    -0.002963787793, &
    -0.000574041191, &
    -0.004442973688, &
    -0.001758110073, &
    -0.005103031982, &
    -0.006578238930, &
    -0.004371465134, &
    -0.005630175684, &
    -0.004743123268, &
    -0.003036162925, &
    -0.004710940904, &
    -0.003766809935, &
    -0.002673280233, &
    -0.003419105699, &
    -0.003972574975, &
    -0.000139297682, &
    -0.001733573296, &
    -0.006375867893, &
    -0.005858755974, &
    -0.003944967467, &
    -0.002302196052, &
    -0.001738802028, &
    -0.001256552044, &
    -0.001186089073, &
    -0.001509544849, &
    -0.002045416728, &
    -0.006489467452, &
    -0.007469535044, &
    -0.005709355473, &
    -0.003486000220, &
    -0.001866862181, &
    -0.002888446030, &
    -0.003571046727, &
    -0.002736100446, &
    -0.004427310453, &
    -0.004319609953, &
    -0.001591471378, &
    -0.003166193922, &
    -0.004049882627, &
    -0.002856209044, &
    -0.003491796907, &
    -0.004382752224, &
    -0.002529931148, &
    -0.004101912641, &
    -0.003418039787, &
    -0.004713178801, &
    -0.002553583687, &
    -0.001762365922, &
    -0.003467423195, &
    -0.003775482249, &
    -0.002733943137, &
    -0.003213103669, &
    -0.005602431233, &
    -0.002231146168, &
    -0.003192121862, &
    -0.000865638847, &
    -0.002106560786, &
    -0.000700644487, &
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000, & 
    0.000000000000 /)

end function asc_Ag_reference

pure function asc_B3g_reference() result(asc_B3g)

  use, intrinsic :: iso_c_binding, only: c_double, c_size_t
  real(c_double) :: asc_B3g(576)

  asc_B3g = (/ 0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    -0.050532399154, &
    -0.009551973174, &
    -0.071962042607, &
    -0.010643147094, &
    -0.089783373467, &
    -0.026341776966, &
    -0.091266072957, &
    -0.010662710196, &
    -0.027356074353, &
    -0.075647252100, &
    -0.012239510683, &
    -0.023321910317, &
    -0.023579909169, &
    -0.055103832662, &
    -0.029685704578, &
    -0.002608543194, &
    -0.041219826083, &
    -0.007025525578, &
    -0.047994090188, &
    -0.025677753050, &
    -0.044269146440, &
    -0.020378168898, &
    -0.024350440826, &
    -0.034219175277, &
    -0.079259123585, &
    -0.132830670580, &
    -0.211794178476, &
    -0.130275453807, &
    -0.092727593074, &
    -0.000900384356, &
    -0.010616419191, &
    -0.040408973071, &
    -0.041902118584, &
    -0.031068183493, &
    -0.056669476873, &
    -0.042970206780, &
    -0.032940883479, &
    -0.032495171991, &
    -0.042068531698, &
    -0.010473528927, &
    -0.042347806882, &
    -0.046158437117, &
    -0.046013993667, &
    -0.040859782010, &
    -0.008632259612, &
    -0.011817693694, &
    -0.015399367973, &
    -0.011144478844, &
    -0.017582319314, &
    -0.014744833006, &
    -0.007394676012, &
    -0.013954217983, &
    -0.015037574239, &
    -0.012429516061, &
    -0.013965993779, &
    -0.017524617968, &
    -0.010987983891, &
    -0.016694727259, &
    -0.010877261182, &
    -0.032651752502, &
    -0.012526764138, &
    -0.009237600215, &
    -0.017692253347, &
    -0.013515632316, &
    -0.017770630907, &
    -0.018351242050, &
    -0.025639062891, &
    -0.010865330530, &
    -0.013561584119, &
    -0.002834710324, &
    -0.009548932251, &
    -0.002495370219, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000 /)

end function asc_B3g_reference

pure function asc_neq_B3g_reference() result(asc_neq_B3g)

  use, intrinsic :: iso_c_binding, only: c_double, c_size_t
  real(c_double) :: asc_neq_B3g(576)

  asc_neq_B3g = (/ 0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    -0.050532399154, &
    -0.009551973174, &
    -0.071962042607, &
    -0.010643147094, &
    -0.089783373467, &
    -0.026341776966, &
    -0.091266072957, &
    -0.010662710196, &
    -0.027356074353, &
    -0.075647252100, &
    -0.012239510683, &
    -0.023321910317, &
    -0.023579909169, &
    -0.055103832662, &
    -0.029685704578, &
    -0.002608543194, &
    -0.041219826083, &
    -0.007025525578, &
    -0.047994090188, &
    -0.025677753050, &
    -0.044269146440, &
    -0.020378168898, &
    -0.024350440826, &
    -0.034219175277, &
    -0.079259123585, &
    -0.132830670580, &
    -0.211794178476, &
    -0.130275453807, &
    -0.092727593074, &
    -0.000900384356, &
    -0.010616419191, &
    -0.040408973071, &
    -0.041902118584, &
    -0.031068183493, &
    -0.056669476873, &
    -0.042970206780, &
    -0.032940883479, &
    -0.032495171991, &
    -0.042068531698, &
    -0.010473528927, &
    -0.042347806882, &
    -0.046158437117, &
    -0.046013993667, &
    -0.040859782010, &
    -0.008632259612, &
    -0.011817693694, &
    -0.015399367973, &
    -0.011144478844, &
    -0.017582319314, &
    -0.014744833006, &
    -0.007394676012, &
    -0.013954217983, &
    -0.015037574239, &
    -0.012429516061, &
    -0.013965993779, &
    -0.017524617968, &
    -0.010987983891, &
    -0.016694727259, &
    -0.010877261182, &
    -0.032651752502, &
    -0.012526764138, &
    -0.009237600215, &
    -0.017692253347, &
    -0.013515632316, &
    -0.017770630907, &
    -0.018351242050, &
    -0.025639062891, &
    -0.010865330530, &
    -0.013561584119, &
    -0.002834710324, &
    -0.009548932251, &
    -0.002495370219, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000, &
    0.000000000000 /)  

end function asc_neq_B3g_reference

