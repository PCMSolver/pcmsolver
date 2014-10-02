module pcmsolver_host_input
   
use, intrinsic :: iso_c_binding

implicit none

public host_input

! cavity specification
character(kind=c_char, len=8), public, save :: pcmsolver_cavity_type = 'gepol  '//c_null_char
integer(c_int), public, save  :: pcmsolver_patch_level = 2
real(c_double), public, save  :: pcmsolver_coarsity = 0.5
real(c_double), public, save  :: pcmsolver_cavity_area = 0.3
real(c_double), public, save  :: pcmsolver_min_distance = 0.1
integer(c_int), public, save  :: pcmsolver_der_order = 4
logical(c_bool), public, save :: pcmsolver_scaling = .true. 
character(kind=c_char, len=8), public, save :: pcmsolver_radii_set = 'bondi  '//c_null_char
character(kind=c_char, len=20), public, save :: pcmsolver_restart_name = 'cavity.npz         '//c_null_char
real(c_double), public, save  :: pcmsolver_min_radius = 100.0
! solver specification
character(kind=c_char, len=7), public, save :: pcmsolver_solver_type = 'iefpcm'//c_null_char
character(kind=c_char, len=16), public, save :: pcmsolver_solvent = 'water          '//c_null_char
character(kind=c_char, len=11), public, save :: pcmsolver_equation_type = 'secondkind'//c_null_char
real(c_double), public, save  :: pcmsolver_correction = 0.0
real(c_double), public, save  :: pcmsolver_probe_radius = 1.0
! green specification
character(kind=c_char, len=7), public, save :: pcmsolver_inside_type = 'vacuum'//c_null_char
character(kind=c_char, len=22), public, save :: pcmsolver_outside_type = 'uniformdielectric    '//c_null_char
real(c_double), public, save :: pcmsolver_outside_epsilon = 1.0 

type, bind(c) :: cavityInput
        character(kind=c_char) :: cavity_type(8)
        integer(c_int) :: patch_level
        real(c_double) :: coarsity
        real(c_double) :: area
        real(c_double) :: min_distance
        integer(c_int) :: der_order
        logical(c_bool) :: scaling
        character(kind=c_char) :: restart_name(20)
        character(kind=c_char) :: radii_set(8)
        real(c_double) :: min_radius
end type

type, bind(c) :: solverInput
        character(kind=c_char) :: solver_type(7)
        character(kind=c_char) :: solvent(16)
        character(kind=c_char) :: equation_type(11) 
        real(c_double)         :: correction
        real(c_double)         :: probe_radius
end type

type, bind(c) :: greenInput
        character(kind=c_char) :: inside_type(7)
        real(c_double)         :: outside_epsilon
        character(kind=c_char) :: outside_type(22)
end type

contains

   subroutine host_input(cavity, solver, green) bind(c, name='host_input')
   ! Performs syntactic checks on PCMSolver input and fills the data structures
   ! holding input data
   ! WARNING ! Do NOT change the order in which strings are passed since the
   ! module relies on that
                                                         
   type(cavityInput) :: cavity
   type(solverInput) :: solver
   type(greenInput)  :: green

   cavity%patch_level  = pcmmod_patch_level
   cavity%coarsity     = pcmmod_coarsity
   cavity%area         = pcmmod_cavity_area
   cavity%min_distance = pcmmod_min_distance
   cavity%der_order    = pcmmod_der_order
   cavity%scaling      = pcmmod_scaling
   cavity%min_radius   = pcmmod_min_radius
   ! Pass the strings relevant to the cavity input section
   call push_input_string(pcmmod_cavity_type//c_null_char)
   call push_input_string(pcmmod_radii_set//c_null_char)
   call push_input_string(pcmmod_restart_name//c_null_char)
   
   solver%correction        = pcmmod_correction
   solver%probe_radius      = pcmmod_probe_radius
   ! Pass the strings relevant to the solver input section
   call push_input_string(pcmmod_solver_type//c_null_char)
   call push_input_string(pcmmod_solvent//c_null_char)
   call push_input_string(pcmmod_equation_type//c_null_char)

   green%outside_epsilon = pcmmod_outside_epsilon
   ! Pass the strings relevant to the green input section
   call push_input_string(pcmmod_inside_type//c_null_char)
   call push_input_string(pcmmod_outside_type//c_null_char)
                                                         
   end subroutine host_input

end module pcmsolver_host_input
