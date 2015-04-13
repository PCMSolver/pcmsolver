#ifndef MACROS_H
#define MACROS_H
/* This file: provides macros for the program */

#define       CRITICAL_ERROR                  {printf("%s(): Critical error, program exit!\n",__func__); exit(-1);}
#define       INVALID_ARG_EXIT                  {printf("Invalid argument passed: Program exit!\n"); exit(-1);}
#define       ALLOC_FAIL_EXIT                   {printf("Memory allocation failure: Program exit!\n"); exit(-1);}
#define       FATAL_EXIT                        {printf("Fatal error: Program exit!\n"); exit(-1);}
#define       PRINT_FUNC_NAME                   {printf("%s()\n",__func__);}
#define       DUMMY_ROUTINE_WARNING             {printf("%s(): Warning this is a dummy routine!\n",__func__);}
#define       CASE_NOT_IMP(X)                   {printf("%s():%s%d%s\n",__func__,"Case ",(X)," not yet implemented!");}
#define       FPRINT_HLINE(X)                   {fprintf((X),"--------------------------------------------------------------------------------\n");}
#define       my_floor(X)                       (floor(X))
#define       same_node_size(X,Y)               ( ( ((X)->dim == (Y)->dim) && ((X)->order == (Y)->order) ) ? (1):(-1) )
#define       max(X,Y)                          (((X) > (Y)) ? (X):(Y))
#define       min(X,Y)                          (((X) < (Y)) ? (X):(Y))
#define       square(X)                         ((X) * (X))
#define       cube(X)                           ((X) * (X) * (X))
#define       out_of_bounds(X,Y,Z)              ((X) < (Y) || (X) > (Z))
#define       alloc_check(X,Y)                  {if(X == NULL) {printf(Y); printf("Memory allocation failure: Program exit!\n"); exit(-1);}}

#ifndef DEBUG
#define PRINT(...)
#define FPRINT(...)
#else
#define PRINT(...) printf(__VA_ARGS__);
#define FPRINT(...) fprintf(__VA_ARGS__);
#endif

extern int DEBUG_LEVEL;
#define debug(level, ...) if (DEBUG_LEVEL < level) printf(__VA_ARGS__);
#define fdebug(level, ...) if (DEBUG_LEVEL < level) fprintf(__VA_ARGS__);

#endif

