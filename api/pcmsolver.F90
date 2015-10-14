module pcmsolver

    use, intrinsic :: iso_c_binding

    implicit none

    private

    public pcmsolver_new
    public pcmsolver_delete
    public pcmsolver_is_compatible_library
    public pcmsolver_print
    public pcmsolver_get_cavity_size
    public pcmsolver_get_irreducible_cavity_size
    public pcmsolver_get_centers
    public pcmsolver_get_center
    public pcmsolver_compute_asc
    public pcmsolver_compute_response_asc
    public pcmsolver_compute_polarization_energy
    public pcmsolver_get_surface_function
    public pcmsolver_set_surface_function
    public pcmsolver_save_surface_functions
    public pcmsolver_save_surface_function
    public pcmsolver_load_surface_function
    public pcmsolver_write_timings

    type, public, bind(C) :: PCMInput
        character(kind=c_char, len=1) :: cavity_type(8)
        integer(c_int)                :: patch_level
        real(c_double)                :: coarsity
        real(c_double)                :: area
        real(c_double)                :: min_distance
        integer(c_int)                :: der_order
        logical(c_bool)               :: scaling
        character(kind=c_char, len=1) :: radii_set(8)
        character(kind=c_char, len=1) :: restart_name(20)
        real(c_double)                :: min_radius
        character(kind=c_char, len=1) :: solver_type(7)
        character(kind=c_char, len=1) :: solvent(16)
        character(kind=c_char, len=1) :: equation_type(11)
        real(c_double)                :: correction
        real(c_double)                :: probe_radius
        character(kind=c_char, len=1) :: inside_type(7)
        real(c_double)                :: outside_epsilon
        character(kind=c_char, len=1) :: outside_type(22)
    end type PCMInput

    public PCMSOLVER_READER_OWN
    public PCMSOLVER_READER_HOST

    enum, bind(C)
        enumerator :: PCMSOLVER_READER_OWN = 0, PCMSOLVER_READER_HOST = 1
    end enum

    interface pcmsolver_new
        function pcmsolver_new(input_reading, nr_nuclei, charges, coordinates, symmetry_info, host_input) result(context) bind(C)
            import
            integer(c_int), intent(in), value :: input_reading
            integer(c_int), intent(in), value :: nr_nuclei
            real(c_double), intent(in)        :: charges(*)
            real(c_double), intent(in)        :: coordinates(*)
            integer(c_int), intent(in)        :: symmetry_info(*)
            type(PCMInput), intent(in), value :: host_input
            type(c_ptr) :: context
        end function pcmsolver_new
    end interface pcmsolver_new

    interface pcmsolver_delete
        subroutine pcmsolver_delete(context) bind(C)
            import
            type(c_ptr), value :: context
        end subroutine pcmsolver_delete
    end interface pcmsolver_delete

    interface pcmsolver_is_compatible_library
        function pcmsolver_is_compatible_library() result(compatible) bind(C)
            import
            logical(c_bool) :: compatible
        end function pcmsolver_is_compatible_library
    end interface pcmsolver_is_compatible_library

    interface pcmsolver_print
        subroutine pcmsolver_print(context) bind(C)
            import
            type(c_ptr), value :: context
        end subroutine pcmsolver_print
    end interface pcmsolver_print

    interface pcmsolver_get_cavity_size
        function pcmsolver_get_cavity_size(context) result(nr_points) bind(C)
            import
            type(c_ptr), value :: context
            integer(c_size_t)  :: nr_points
        end function pcmsolver_get_cavity_size
    end interface pcmsolver_get_cavity_size

    interface pcmsolver_get_irreducible_cavity_size
        function pcmsolver_get_irreducible_cavity_size(context) result(nr_points_irr) bind(C)
            import
            type(c_ptr), value :: context
            integer(c_size_t)  :: nr_points_irr
        end function pcmsolver_get_irreducible_cavity_size
    end interface pcmsolver_get_irreducible_cavity_size

    interface pcmsolver_get_centers
        subroutine pcmsolver_get_centers(context, centers) bind(C)
            import
            type(c_ptr), value :: context
            real(c_double), intent(inout) :: centers(*)
        end subroutine pcmsolver_get_centers
    end interface pcmsolver_get_centers

    interface pcmsolver_get_center
        subroutine pcmsolver_get_center(context, its, center) bind(C)
            import
            type(c_ptr), value :: context
            integer(c_int), value, intent(in) :: its
            real(c_double), intent(inout) :: center(*)
        end subroutine pcmsolver_get_center
    end interface pcmsolver_get_center

    interface pcmsolver_compute_asc
        subroutine pcmsolver_compute_asc(context, mep_name, asc_name, irrep) bind(C)
            import
            type(c_ptr), value :: context
            character(c_char), intent(in) :: mep_name, asc_name
            integer(c_int), value, intent(in) :: irrep
        end subroutine pcmsolver_compute_asc
    end interface pcmsolver_compute_asc

    interface pcmsolver_compute_response_asc
        subroutine pcmsolver_compute_response_asc(context, mep_name, asc_name, irrep) bind(C)
            import
            type(c_ptr), value :: context
            character(c_char), intent(in) :: mep_name, asc_name
            integer(c_int), value, intent(in) :: irrep
        end subroutine pcmsolver_compute_response_asc
    end interface pcmsolver_compute_response_asc

    interface pcmsolver_compute_polarization_energy
        function pcmsolver_compute_polarization_energy(context, mep_name, asc_name) result(energy) bind(C)
            import
            type(c_ptr), value :: context
            character(c_char), intent(in) :: mep_name, asc_name
            real(c_double) :: energy
        end function pcmsolver_compute_polarization_energy
    end interface pcmsolver_compute_polarization_energy

    interface pcmsolver_get_surface_function
        subroutine pcmsolver_get_surface_function(context, f_size, values, name) bind(C)
            import
            type(c_ptr), value :: context
            integer(c_size_t), value, intent(in) :: f_size
            real(c_double), intent(inout) :: values(*)
            character(c_char), intent(in) :: name
        end subroutine pcmsolver_get_surface_function
    end interface pcmsolver_get_surface_function

    interface pcmsolver_set_surface_function
        subroutine pcmsolver_set_surface_function(context, f_size, values, name) bind(C)
            import
            type(c_ptr), value :: context
            integer(c_size_t), value, intent(in) :: f_size
            real(c_double), intent(in) :: values(*)
            character(c_char), intent(in) :: name
        end subroutine pcmsolver_set_surface_function
    end interface pcmsolver_set_surface_function

    interface pcmsolver_save_surface_functions
        subroutine pcmsolver_save_surface_functions(context) bind(C)
            import
            type(c_ptr), value :: context
        end subroutine pcmsolver_save_surface_functions
    end interface pcmsolver_save_surface_functions

    interface pcmsolver_save_surface_function
        subroutine pcmsolver_save_surface_function(context, name) bind(C)
            import
            type(c_ptr), value :: context
            character(c_char), intent(in) :: name
        end subroutine pcmsolver_save_surface_function
    end interface pcmsolver_save_surface_function

    interface pcmsolver_load_surface_function
        subroutine pcmsolver_load_surface_function(context, name) bind(C)
            import
            type(c_ptr), value :: context
            character(c_char), intent(in) :: name
        end subroutine pcmsolver_load_surface_function
    end interface pcmsolver_load_surface_function

    interface pcmsolver_write_timings
        subroutine pcmsolver_write_timings(context) bind(C)
            import
            type(c_ptr), value :: context
        end subroutine pcmsolver_write_timings
    end interface pcmsolver_write_timings

end module pcmsolver
