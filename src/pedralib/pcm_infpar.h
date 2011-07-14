#ifdef COMMENT
! -- infpar.h --
!     my_MPI_INTEGER is used both in .c and .F routines in MPI calls
!        so we can handle "-i8" compilations on 32-bit machines,
!        using VAR_INT64 /Jan-2007 hjaaj
#endif
#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

#if defined(__CVERSION__)
#define MAXNOD 200
#define MAXCL2 10000
#define NPARI  ((MAXNOD + 1) + 6)
extern struct common_infpar {
#if defined (VAR_INT64)
    long iprpar, ntask, ncode, ndegdi, master, mynum, mytid;
    long nodtot, nodeid[MAXNOD+1], nfmat, mtottk, parher, debug, pario;
    long timing, slave;
#else
    int  iprpar, ntask, ncode, ndegdi, master, mynum, mytid;
    int  nodtot, nodeid[MAXNOD+1], nfmat, mtottk, parher, debug, pario;
    int  timing, slave;
#endif
    char nodnam[MAXNOD][20], myname[20];
} daltoninfpar_;
#else
C File: infpar.h for Dalton; special information for parallel calculations
C
C     Parameters NPARI must be updated after changes (for parallelization)
C
C     NOTE: Integers  (IPRPAR,...,MASTER,...,MYTID)
C           Logicals  (TIMING,SLAVE)
C           Character (NODNAM,MYNAME) should NOT be sent to slaves
C     THUS: NPARI is length from NODTOT,...,PARIO
C
      INTEGER MAXNOD, MAXCL2
      PARAMETER (MAXNOD = 200, MAXCL2 = 10000)
      PARAMETER (NPARI = (MAXNOD + 1) + 6)
      INTEGER IPRPAR, NTASK, NCODE, NDEGDI, MASTER, MYNUM, MYTID
      INTEGER NODTOT, NODEID(0:MAXNOD), NFMAT, MTOTTK
      LOGICAL PARHER, PARIO, DEBUG,     TIMING, SLAVE
      CHARACTER*20   NODNAM(0:MAXNOD), MYNAME
      COMMON /DALTONINFPAR/                                              &
     &        IPRPAR, NTASK, NCODE, NDEGDI, MASTER, MYNUM, MYTID         &
     &       ,NODTOT, NODEID, NFMAT, MTOTTK, PARHER, DEBUG, PARIO        &
     &       ,TIMING, SLAVE , NODNAM, MYNAME

#if defined (VAR_INT64)
!     integer array ISTAT contains MPI_SOURCE information.
!     Proper use of ISTAT on 64-bit machines in 
!     combination with VAR_INT64 requires explicit declaration 
!     as INTEGER*4 /March-2007 sk 
      INTEGER*4 ISTAT
#endif

C -- end of infpar.h --
#endif
