module pcmsolver

    implicit none

    private

    public pcmsolver_new
    public pcmsolver_delete
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

    interface pcmsolver_new
        function pcmsolver_new(F_collect_nctot, F_collect_atoms, F_host_writer, F_set_point_group) result(context) bind (C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_funptr
            type(c_funptr), value, intent(in) :: F_collect_nctot
            type(c_funptr), value, intent(in) :: F_collect_atoms
            type(c_funptr), value, intent(in) :: F_host_writer
            type(c_funptr), value, intent(in) :: F_set_point_group
            type(c_ptr) :: context
        end function
    end interface

    interface pcmsolver_delete
        subroutine pcmsolver_delete(context) bind (C)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr), value :: context
        end subroutine
    end interface

    interface pcmsolver_get_cavity_size
        function pcmsolver_get_cavity_size(context) result(nr_points) bind (C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t)  :: nr_points
        end function pcmsolver_get_cavity_size
    end interface pcmsolver_get_cavity_size

    interface pcmsolver_get_irreducible_cavity_size
        function pcmsolver_get_irreducible_cavity_size(context) result(nr_points_irr) bind (C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t)  :: nr_points_irr
        end function pcmsolver_get_irreducible_cavity_size
    end interface pcmsolver_get_irreducible_cavity_size

    interface pcmsolver_get_centers
        subroutine pcmsolver_get_centers(context, centers) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_double
            type(c_ptr), value :: context
            real(c_double), intent(inout) :: centers(*)
        end subroutine pcmsolver_get_centers
    end interface pcmsolver_get_centers

    interface pcmsolver_get_center
        subroutine pcmsolver_get_center(context, its, center) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
            type(c_ptr), value :: context
            real(c_double), value, intent(in) :: its
            real(c_double), intent(inout) :: center(*)
        end subroutine pcmsolver_get_center
    end interface pcmsolver_get_center

    interface pcmsolver_compute_asc
        subroutine pcmsolver_compute_asc(context, mep_name, asc_name, irrep) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
            type(c_ptr), value :: context
            character(c_char), intent(in) :: mep_name, asc_name
            integer(c_int), value, intent(in) :: irrep
        end subroutine pcmsolver_compute_asc
    end interface pcmsolver_compute_asc

    interface pcmsolver_compute_response_asc
        subroutine pcmsolver_compute_response_asc(context, mep_name, asc_name, irrep) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
            type(c_ptr), value :: context
            character(c_char), intent(in) :: mep_name, asc_name
            integer(c_int), value, intent(in) :: irrep
        end subroutine pcmsolver_compute_response_asc
    end interface pcmsolver_compute_response_asc

    interface pcmsolver_compute_polarization_energy
        subroutine pcmsolver_compute_polarization_energy(context, mep_name, asc_name) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_char
            type(c_ptr), value :: context
            character(c_char), intent(in) :: mep_name, asc_name
        end subroutine pcmsolver_compute_polarization_energy
    end interface pcmsolver_compute_polarization_energy

    interface pcmsolver_get_surface_function
        subroutine pcmsolver_get_surface_function(context, f_size, values, name) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_size_t, c_double
            type(c_ptr), value :: context
            integer(c_size_t), value, intent(in) :: f_size
            real(c_double), intent(inout) :: values
            character(c_char), intent(in) :: name
        end subroutine pcmsolver_get_surface_function
    end interface pcmsolver_get_surface_function

    interface pcmsolver_set_surface_function
        subroutine pcmsolver_set_surface_function(context, f_size, values, name) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_size_t, c_double
            type(c_ptr), value :: context
            integer(c_size_t), value, intent(in) :: f_size
            real(c_double), intent(in) :: values
            character(c_char), intent(in) :: name
        end subroutine pcmsolver_set_surface_function
    end interface pcmsolver_set_surface_function

    interface pcmsolver_save_surface_functions
        subroutine pcmsolver_save_surface_functions(context) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr), value :: context
        end subroutine pcmsolver_save_surface_functions
    end interface pcmsolver_save_surface_functions

    interface pcmsolver_save_surface_function
        subroutine pcmsolver_save_surface_function(context, name) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_char
            type(c_ptr), value :: context
            character(c_char), intent(in) :: name
        end subroutine pcmsolver_save_surface_function
    end interface pcmsolver_save_surface_function

    interface pcmsolver_load_surface_function
        subroutine pcmsolver_load_surface_function(context, name) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_char
            type(c_ptr), value :: context
            character(c_char), intent(in) :: name
        end subroutine pcmsolver_load_surface_function
    end interface pcmsolver_load_surface_function

    interface pcmsolver_write_timings
        subroutine pcmsolver_write_timings(context) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr), value :: context
        end subroutine pcmsolver_write_timings
    end interface pcmsolver_write_timings

end module pcmsolver
