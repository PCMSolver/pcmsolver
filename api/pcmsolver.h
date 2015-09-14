#ifndef PCMSOLVER_H_INCLUDED
#define PCMSOLVER_H_INCLUDED

#include <stddef.h>

#ifndef PCMSOLVER_API
#  ifdef _WIN32
#     if defined(PCMSOLVER_BUILD_SHARED) /* build dll */
#         define PCMSOLVER_API __declspec(dllexport)
#     elif !defined(PCMSOLVER_BUILD_STATIC) /* use dll */
#         define PCMSOLVER_API __declspec(dllimport)
#     else /* static library */
#         define PCMSOLVER_API
#     endif
#  else
#     if __GNUC__ >= 4
#         define PCMSOLVER_API __attribute__((visibility("default")))
#     else
#         define PCMSOLVER_API
#     endif
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct pcmsolver_context_s;
typedef struct pcmsolver_context_s pcmsolver_context_t;

typedef int(*collect_nctot)(void);
typedef void(*collect_atoms)(double[], double[]);
typedef void(*host_writer)(const char *, size_t);
typedef void(*set_point_group)(int *, int *, int *, int *);
//typedef void(*host_input)(cavityInput &, solverInput &, greenInput &);

typedef enum
{
    PCMSOLVER_READER_OWN,
    PCMSOLVER_READER_HOST
} pcmsolver_reader_t;

PCMSOLVER_API pcmsolver_context_t * pcmsolver_new(collect_nctot, collect_atoms, host_writer, set_point_group);

PCMSOLVER_API void pcmsolver_delete(pcmsolver_context_t * context);

PCMSOLVER_API void pcmsolver_print(pcmsolver_context_t * context);

PCMSOLVER_API size_t pcmsolver_get_cavity_size(pcmsolver_context_t * context);

PCMSOLVER_API size_t pcmsolver_get_irreducible_cavity_size(pcmsolver_context_t * context);

PCMSOLVER_API void pcmsolver_get_centers(pcmsolver_context_t * context, double centers[]);

PCMSOLVER_API void pcmsolver_get_center(pcmsolver_context_t * context, int its, double center[]);

PCMSOLVER_API void pcmsolver_compute_asc(pcmsolver_context_t * context,
                                        const char * mep_name,
                                        const char * asc_name,
                                        int irrep);

PCMSOLVER_API void pcmsolver_compute_response_asc(pcmsolver_context_t * context,
                                                 const char * mep_name,
                                                 const char * asc_name,
                                                 int irrep);

PCMSOLVER_API double pcmsolver_compute_polarization_energy(pcmsolver_context_t * context,
                                             const char * mep_name,
                                             const char * asc_name);

PCMSOLVER_API void pcmsolver_get_surface_function(pcmsolver_context_t * context,
                                                 size_t size,
                                                 double values[],
                                                 const char * name);

PCMSOLVER_API void pcmsolver_set_surface_function(pcmsolver_context_t * context,
                                                 size_t size,
                                                 double values[],
                                                 const char * name);

PCMSOLVER_API void pcmsolver_save_surface_functions(pcmsolver_context_t * context);

PCMSOLVER_API void pcmsolver_save_surface_function(pcmsolver_context_t * context,
                                                   const char * name);

PCMSOLVER_API void pcmsolver_load_surface_function(pcmsolver_context_t * context,
                                                   const char * name);

PCMSOLVER_API void pcmsolver_write_timings(pcmsolver_context_t * context);

#ifdef __cplusplus
}
#endif

#endif /* PCMSOLVER_H_INCLUDED */
