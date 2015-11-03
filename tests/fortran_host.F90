program pcm_fortran_host

      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env, only: output_unit
      use pcmsolver
      
      implicit none
      
      interface pcmsolver_input
        function pcmsolver_input() result(host_input)
          use pcmsolver, only: PCMInput
          type(PCMInput) :: host_input
        end function pcmsolver_input
      end interface pcmsolver_input

      interface nuclear_mep
        pure function nuclear_mep(nr_nuclei, charges, coordinates, grid_size, grid) result(n_mep)
          use, intrinsic :: iso_c_binding
          integer(c_int), intent(in) :: nr_nuclei
          real(c_double), intent(in) :: charges(nr_nuclei)
          real(c_double), intent(in) :: coordinates(3, nr_nuclei)
          integer(c_size_t), intent(in) :: grid_size
          real(c_double), intent(in) :: grid(3, grid_size)
          real(c_double) :: n_mep(grid_size)
        end function nuclear_mep
      end interface nuclear_mep
 
      type(c_ptr) :: pcm_context
      integer(c_int) :: nr_nuclei
      real(c_double), allocatable :: charges(:)
      real(c_double), allocatable :: coordinates(:)
      integer(c_int) :: symmetry_info(4)
      type(PCMInput) :: host_input
      logical :: log_open, log_exist
      character(kind=c_char, len=7) :: mep_lbl, asc_lbl, asc_B3g_lbl, oit_asc_lbl
      real(c_double), allocatable :: grid(:), mep(:), asc_Ag(:), asc_B3g(:), oit_asc(:)
      integer(c_int) :: irrep
      integer(c_size_t) :: grid_size, irr_grid_size, ipoint
      real(c_double) :: energy
      character(8)  :: for_title  = '(60X, A)'
      character(36) :: for_header = '(A, T27, A, T62, A, T97, A, T132, A)'
      character(20) :: for_data   = '(I6, 4(20X, F15.12))'
      real(c_double) :: tot_asc

      if (.not. pcmsolver_is_compatible_library()) then
        print *, 'PCMSolver library not compatible!'
        stop
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

      allocate(oit_asc(grid_size))
      oit_asc = 0.0_c_double
      oit_asc_lbl = 'OITASC'//c_null_char
      ! This is the B3g irreducible representation
      irrep = 3_c_size_t
      call pcmsolver_compute_response_asc(pcm_context, mep_lbl, oit_asc_lbl, irrep)
      call pcmsolver_get_surface_function(pcm_context, grid_size, oit_asc, oit_asc_lbl)

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

      write(output_unit, for_title) 'Converged MEP and ASC'
      write(output_unit, for_header) 'Finite element #', 'MEP', 'ASC Ag', 'ASC B3g', 'Nonequilibrium ASC B3g'
      tot_asc = 0.0_c_double
      do ipoint = 1, grid_size
        tot_asc = tot_asc + asc_Ag(ipoint)
        write(output_unit, for_data) ipoint, mep(ipoint), asc_Ag(ipoint), asc_B3g(ipoint), oit_asc(ipoint)
      end do
      write(output_unit, '(A, F20.12)') 'Sum of nuclear apparent charges ', tot_asc

      call pcmsolver_write_timings(pcm_context)

      call pcmsolver_delete(pcm_context)

      deallocate(charges)
      deallocate(coordinates)
      deallocate(grid)
      deallocate(mep)
      deallocate(asc_Ag)
      deallocate(asc_B3g)
      deallocate(oit_asc)

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
