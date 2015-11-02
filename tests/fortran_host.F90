program pcm_fortran_host

      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env, only: output_unit
      use pcmsolver
      
      implicit none
 
      type(c_ptr) :: pcm_context
      integer(c_int) :: nr_nuclei
      real(c_double), allocatable :: charges(:)
      real(c_double), allocatable :: coordinates(:)
      integer(c_int) :: symmetry_info(4)
      type(PCMInput), external :: pcmsolver_input
      type(PCMInput) :: host_input
      logical :: log_open, log_exist

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
      coordinates= (/ 0.0_c_double,  0.0_c_double,       1.257892_c_double, &
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


      call pcmsolver_write_timings(pcm_context)

      call pcmsolver_delete(pcm_context)

      deallocate(charges)
      deallocate(coordinates)

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
  character(kind=c_char, len=6) :: pcmmod_cavity_type = 'gepol'//c_null_char
  integer :: pcmmod_patch_level = 2
  real(c_double) :: pcmmod_coarsity = 0.5
  real(c_double) :: pcmmod_cavity_area = 0.3
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
