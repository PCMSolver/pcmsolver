      SUBROUTINE DZERO_(DX,LENGTH)
#include "pcm_implicit.h"
C                                                                               
C Last revision 5-May-1984 by Hans Jorgen Aa. Jensen                            
C                                                                               
C   Subroutine DZERO_ sets a real array of length *LENGTH*                       
C   to zero.                                                                    
C...................................................................            
      DIMENSION DX(*)
C                                                                               
      IF (LENGTH.LE.0) RETURN
C                                                                               
      DO I = 1,LENGTH
         DX(I) = 0.0D00
      END DO
C                                                                               
      RETURN
      END
      SUBROUTINE AROUND_(HEAD)
      CHARACTER HEAD*(*)
#include <pcm_priunit.h>
      LHEAD  = LNBLNK(HEAD)
      LNG    = LHEAD + 2
      IND = MAX(1,(80 - LNG)/2 + 1)
      WRITE (LVPRI,'(//150A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '! ', HEAD(1:LHEAD), ' !'
      WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
Cx    WRITE (LVPRI,'(//150A)') (' ',I=1,IND), '.', ('-',I=1,LNG), '.'           
Cx    WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '| ', HEAD(1:LHEAD), ' |'         
Cx    WRITE (LVPRI,'(150A)')   (' ',I=1,IND), '`', ('-',I=1,LNG), ''''          
      WRITE (LVPRI,'()')
      RETURN
      END
      SUBROUTINE AVERAG_ (VECS,NVDIM,NVSIM)
C                                                                               
C 7-Jul-1992 hjaaj+hh                                                           
C                                                                               
#include <pcm_implicit.h>
#include <pcm_priunit.h>
      DIMENSION VECS(NVDIM,NVSIM)
C                                                                               
C Used from common blocks:                                                      
C  INFINP : SUPSYM                                                              
C  INFVAR : NWOPT,JWOPSY                                                        
C  INFPRI : LVERR                                                               
C                                                                               
#include <pcm_maxorb.h>
#include <pcm_infinp.h>
#include <pcm_infvar.h>
#include <pcm_infpri.h>
C                                                                               
Ckm      CALL QENTER_('AVERAG_')
      IF (.NOT. SUPSYM) GO TO 9999
      IF (NWOPT  .EQ. 0) GO TO 9999
      IF (JWOPSY .NE. 1) GO TO 910
      IF (NVDIM  .LT. NWOPT) GO TO 920
      IF (SUPSYM) THEN
         CALL AVERSS_(VECS,NVDIM,NVSIM)
      END IF
Ckm 9999 CALL QEXIT('AVERAG_')
 9999 RETURN
C                                                                               
C *** Error sections                                                            
C                                                                               
 910   CONTINUE
      WRITE (LVERR,9010) JWOPSY
Ckm      CALL QTRACE_(LVERR)
      STOP
 9010  FORMAT(/' ERROR-AVERAG_, operator symmetry .ne. 1; JWOPSY =',I8)
C                                                                               
 920    CONTINUE
      WRITE (LVERR,9020) NVDIM,NWOPT
Ckm      CALL QTRACE_(LVERR)
      STOP
 9020  FORMAT(/' ERROR-AVERAG_, vector length',I8,' is less than NWOPT',
     *   I8)
C                                                                               
      END
      SUBROUTINE AVERSS_ (VECS,NVDIM,NVSIM)
C                                                                               
C 7-Jul-1992 Hans Joergen Aa. Jensen + Hinne Hettema                            
C                                                                               
C Purpose: Average degenerate super symmetries                                  
C          (which are different components of a degenerate irrep).              
C          This corresponds to doing a microcanonical average.                  
C                                                                               
#include <pcm_implicit.h>
      DIMENSION VECS(NVDIM,NVSIM)
C                                                                               
C Used from common blocks:                                                      
C   INFORB : NSSYM,NORBSS(),IORBSS(),NINFSS(*,3), NOCCT,NORBT                   
C   INFIND : ISW(*),ISSORD(*)                                                   
C   INFVAR : NWOPT                                                              
C   INFPRI : LVERR                                                              
C                                                                               
#include <pcm_maxash.h>
#include <pcm_maxorb.h>
#include <pcm_priunit.h>
#include <pcm_infind.h>
#include <pcm_inforb.h>
#include <pcm_infvar.h>
#include <pcm_infpri.h>
C                                                                               
      PARAMETER (MAXDEG = 20, D1 = 1.0D0)
      DIMENSION ISSDEG(MAXDEG),ISSOFF(MAXDEG),IDEGV(MAXDEG)
C                                                                               
Chjaaj-may2000:                                                                 
#ifdef NO_F90
      DIMENSION KLWOP(1,1)
      STOP
C if automatic arrays not allowed include #define NO_F90.                       
C This means you cannot use .SUPSYM. /hjaaj may2000                             
#else
C     Dynamic allocation of KLWOP with adjustable run-time dimensions           
      DIMENSION KLWOP(NOCCT,NORBT)
#endif
C                                                                               
Ckm      CALL QENTER_('AVERSS_')
      CALL MAKE_KLWOP_(KLWOP)
      IF (IPRAVE .GT. 100) THEN
         WRITE(LVPRI,'(/A)') ' AVERSS_, Incoming orbital vector(s):'
         CALL OUTPUT_(VECS,1,NWOPT,1,NVSIM,NVDIM,NVSIM,-1,LVPRI)
      END IF
C                                                                               
      DO 800 ISSYMR = 1,NSSYM
C        Skip this supsym if not degenerate or not "root supsym"                
         IF (NINFSS(ISSYMR,1) .EQ. 1) GO TO 800
         IF (NINFSS(ISSYMR,2) .NE. ISSYMR) GO TO 800
         NDEG = 0
         DO 100 ISSYM = ISSYMR,NSSYM
            IF (NINFSS(ISSYM,2) .EQ. ISSYMR) THEN
               NDEG = NDEG + 1
               ISSDEG(NDEG) = ISSYM
               ISSOFF(NDEG) = IORBSS(ISSYM) - IORBSS(ISSYMR)
            END IF
 100            CONTINUE
         IF (IPRAVE .GT. 5 .OR. NDEG .NE. NINFSS(ISSYMR,1)) THEN
            WRITE(LVPRI,'(/A,I3,A/A,(T35,8I5))')
     &         'AVERSS_ averaging ',NDEG,' degenerate super symmetries',
     &         ' namely sup.sym.s no.',(ISSDEG(IDEG),IDEG=1,NDEG)
         END IF
         IF (NDEG .NE. NINFSS(ISSYMR,1)) THEN
            WRITE(LVPRI,'(/A)')
     &      ' AVERSS_-software error: inconsistent degeneracy in NINFSS'
            CALL AVEDMP_(LVPRI)
            STOP
         END IF
C                                                                               
         AVFAC = NDEG
         AVFAC = D1 / AVFAC
C                                                                               
         IOFFSR = IORBSS(ISSYMR)
         DO 480 IMOSS = IOFFSR+1,IOFFSR+NORBSS(ISSYMR)
            IMO = ISSORD(IMOSS)
            IMOW = ISW(IMO)
            IF (IMOW .GT. NOCCT) GO TO 480
            DO 470 JMOSS = IMOSS+1,IOFFSR+NORBSS(ISSYMR)
               JMO = ISSORD(JMOSS)
               JMOW = ISW(JMO) - NISHT
               IF (JMOW .LE. 0) GO TO 470
Chjaaj         ... this is an inactive orbital                                  
               IF (KLWOP(IMOW,JMOW) .EQ. 0) GO TO 470
               IDEGV(1) = KLWOP(IMOW,JMOW)
               NERR = 0
               DO 410 IDEG = 2,NDEG
                  KMO = ISSORD(IMOSS + ISSOFF(IDEG))
                  LMO = ISSORD(JMOSS + ISSOFF(IDEG))
                  LMOW = ISW(LMO) - NISHT
                  IF (LMOW .LE. 0) THEN
                     IDEGV(IDEG) = 0
                     NERR = NERR + 1
                     WRITE(LVPRI,'(/A)')
     &                    'AVERSS_ error: degenerate orbitals not ',
     &                    'same orbital class.'
                  ELSE
                     IDEGV(IDEG) = KLWOP( ISW(KMO) , LMOW )
                     IF (IDEGV(IDEG) .EQ. 0) NERR = NERR + 1
                  END IF
 410                        CONTINUE
               IF (IPRAVE .GT. 10 .OR. NERR .GT. 0) THEN
                  IF (NERR .GT. 0) THEN
                     WRITE(LVPRI,'(/A/A/)')
     &         ' AVERSS_ error: not all degenerate rotations included',
     &         ' Probable error: some but not all'//
     &            ' of the degenerate orbitals have been frozen.'
                  END IF
                  WRITE(LVPRI,'(A,(T35,4I10))')
     &               ' AVERSS_ averaging orb.rot. no.s',
     &               (IDEGV(IDEG),IDEG=1,NDEG)
                  IF (NERR .GT. 0) CALL AVEDMP_(LVPRI)
               END IF
               IF (NERR .GT. 0) THEN
Ckm                  CALL QTRACE_(LVPRI)
                  STOP
               END IF
               DO 460 IVSIM = 1,NVSIM
                  XKAPIJ = VECS(IDEGV(1),IVSIM)
                  DO 420 IDEG = 2,NDEG
                     XKAPIJ = XKAPIJ + VECS(IDEGV(IDEG),IVSIM)
 420                              CONTINUE
                  XKAPIJ = AVFAC*XKAPIJ
                  DO 430 IDEG = 1,NDEG
                     VECS(IDEGV(IDEG),IVSIM) = XKAPIJ
 430                              CONTINUE
 460                                        CONTINUE
 470                                               CONTINUE
 480                                                   CONTINUE
 800                                                    CONTINUE
C                                                                               
      IF (IPRAVE .GT. 100) THEN
         WRITE(LVPRI,'(/A)') ' AVERSS_, Averaged orbital vector(s):'
         CALL OUTPUT_(VECS,1,NWOPT,1,NVSIM,NVDIM,NVSIM,-1,LVPRI)
      END IF
C                                                                               
Ckm      CALL QEXIT('AVERSS_')
      RETURN
      END
      SUBROUTINE MAKE_KLWOP_(KLWOP)
C                                                                               
C 1-May-2000 hjaaj                                                              
C                                                                               
C Purpose: Generate KLWOP array from JWOP                                       
C                                                                               
#include <pcm_implicit.h>
      DIMENSION   KLWOP(NOCCT,NORBT)
      CHARACTER*8 XLABEL
C                                                                               
C Used from common blocks:                                                      
C   INFORB: NOCCT,NORBT                                                         
C   INFIND: ISW(*)                                                              
C   INFVAR: NWOPT,JWOP(2,*)                                                     
C                                                                               
#include <pcm_priunit.h>
#include <pcm_maxash.h>
#include <pcm_maxorb.h>
#include <pcm_inforb.h>
#include <pcm_infind.h>
#include <pcm_infvar.h>
C                                                                               
      NKLWOP = NOCCT*NORBT
      CALL IZERO_(KLWOP,NKLWOP)
C                                                                               
      DO IG = 1,NWOPT
         KW = ISW(JWOP(1,IG))
         LW = ISW(JWOP(2,IG))-NISHT
         IF (KLWOP(KW,LW) .EQ. 0) THEN
            KLWOP(KW,LW)  = IG
         ELSE
            WRITE (LVPRI,'(/A,I10,A,I10//A)')
     &      'Fatal error, duplicate entry in JWOP found for elements',
     &      KLWOP(KW,LW),' and',IG,'Dump of JWOP:'
            WRITE (LVPRI,'(3I10)') (I,JWOP(1,I),JWOP(2,I),I=1,NWOPT)
            STOP
         END IF
      END DO
C                                                                               
      RETURN
      END
      SUBROUTINE AVEDMP_(LUWDMP)
C                                                                               
C     2-Jul-92 Hinne Hettema                                                    
C     Purpose : dump supersymmetry information to output for debugging          
C                                                                               
#include <pcm_implicit.h>
C                                                                               
#include <pcm_maxash.h>
#include <pcm_maxorb.h>
#include <pcm_inforb.h>
#include <pcm_infind.h>
C                                                                               
Ckm      CALL QENTER_('AVEDMP_')
      WRITE(LUWDMP,'(/A/A)')
     *   '    Dump of super symmetry information',
     *   ' ========================================'
      WRITE(LUWDMP,'(/A,I4)')
     *   ' Number of supersymmetries (NSSYM) is ',NSSYM
      WRITE(LUWDMP,'(/A/A)')
     *   '  ISS    NORBSS    IORBSS',
     *   ' =============================='
      DO 100 I = 1, NSSYM
         WRITE(LUWDMP,'(1X,I4,2I10)') I, NORBSS(I), IORBSS(I)
 100      CONTINUE
      WRITE(LUWDMP,'(A)') ' =============================='
      WRITE(LUWDMP,'(/A/A)')
     *   '    I     ISSMO    ISSORD      ISMO',
     *   ' ========================================'
      DO 200 I = 1, NORBT
         WRITE(LUWDMP,'(1X,I4,3I10)') I, ISSMO(I), ISSORD(I), ISMO(I)
 200      CONTINUE
      WRITE(LUWDMP,'(A/)') ' ========================================'
Ckm      CALL QEXIT('AVEDMP_')
      RETURN
      END
      SUBROUTINE IZERO_(INTVEC,LENGTH)
C...................................................................            
C Written 5-May-1984 by Hans Jorgen Aa. Jensen                                  
C                                                                               
C   Subroutine IZERO_ sets an integer array of length *LENGTH*                   
C   to zero.                                                                    
C...................................................................            
      INTEGER LEGNTH, INTVEC(*)
C                                                                               
      IF (LENGTH.LE.0) RETURN
C                                                                               
      DO I=1,LENGTH
         INTVEC(I) = 0
      END DO
C                                                                               
      RETURN
      END
      FUNCTION DNORM2__(N,DX,INCX)
C                                                                               
C     Forms the two-norm of a vector.                                           
C 19-Sep-1988 -- hjaaj -- based on DNRM2_ from LINPACK                           
C     This version does not use extended precision for intermediates            
C     as the LINPACK version does.C 1) DNORM2__ (emulate ESSL DNORM2: do not use \
Cextended precision for intermediates)                                           

C     Equivalent to DNORM2__ in IBM's ESSL library.                               
C                                                                               
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                          
C     DNRM2_: JACK DONGARRA, LINPACK, 3/11/78.                                   
C                                                                               
#include <pcm_implicit.h>
C                                                                               
      DIMENSION DX(*)
      PARAMETER ( ZERO = 0.0D0 )
C                                                                               
      IF (N.LE.0) THEN
         DNORM2__ = ZERO
         RETURN
      END IF
      DTEMP  = ZERO
      IF(INCX.EQ.1)GO TO 20
C                                                                               
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                        
C          NOT EQUAL TO 1                                                       
C                                                                               
      IX = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DX(IX)
        IX = IX + INCX
 10      CONTINUE
      DNORM2__ = SQRT(DTEMP)
      RETURN
C                                                                               
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                                    
C                                                                               
C                                                                               
C        CLEAN-UP LOOP                                                          
C                                       20 M = MOD(N,5)
 20   M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DX(I)
 30      CONTINUE
      IF( N .LT. 5 ) GO TO 60
 40    MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DX(I) + DX(I + 1)*DX(I + 1) +
     *   DX(I + 2)*DX(I + 2) + DX(I + 3)*DX(I + 3) + DX(I + 4)*DX(I + 4)
 50      CONTINUE
 60       DNORM2__ = SQRT(DTEMP)
      RETURN
      END
      SUBROUTINE GPOPEN_(LUNIT,FILEIN,STATIN,ACCEIN,FORMIN,NELEM,OLDDX)
C
C 23-Feb-2000 K.Ruud
C
C Purpose:
C   General purpose routine for opening files in the Dalton program.
C   The routine will dynamically allocate unit numbers that will
C   become available again when the file is closed using GPCLOSE_
C
C   Direct access files are provided with a UNIT-number, but the rest
C   of the file opening process is taken care of by a call to OPENDX_.
C
C   These files are strongly machine dependent, although care has been
C   taken to avoid using unit numbers that are illegal or reserved on a
C   particular architecture. However, the use of GPOPEN_ and GPCLSE should
C   remove much of the machine dependence in the rest of the Dalton
C   program
C
C Input:
C  LUNIT    Suggested unit number (OPTIONAL, but is mandatory if an unnamed
C           file is reopened after having been closed with STATUS='KEEP')
C           Otherwise it is recommended to not assign this.
C  FILEIN   Suggested name for the file (OPTIONAL, but strongly recommended)
C  STATIN   Suggested status of the file (OPTIONAL and maybe not recommended)
C  ACCEIN   Access method. 'DIRECT' or 'SEQUENTIAL'. Default is sequential.
C  FORMIN   Formatted or unformatted file format. Default is 'UNFORMATTED'
C
C Direct access-specific input:
C  NELEM    number of integer words per logical record
C  NREC     number of logical records
C
C
C Output:
C  LUNIT    Assigned file unit number
C
C Direct access-specific output:
C  OLDDX    logical, true if old LUDX file was opened
C
C
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_chrnos.h>
      CHARACTER*80 STATUF, ACCESF, FORMF
      CHARACTER*(*) FILEIN, STATIN, ACCEIN, FORMIN
      integer max_filename_length
      parameter (max_filename_length = 500)
      CHARACTER OUTFIL*(max_filename_length)
      CHARACTER WRKDIR*(max_filename_length)
      CHARACTER FILENM*(max_filename_length)
      INTEGER FILELN, STATLN, ACCELN, FORMLN
      LOGICAL OLDDX, FEXST, FIRST
      COMMON /PCM_UNITAR/ IUNTAB(1:99)
#include <pcm_gnrinf.h>
#include <pcm_infpar.h>
      DATA FIRST /.TRUE./
      SAVE FIRST
C
Ckm      CALL QENTER_('GPOPEN_')
C
C     Initialization
C     
      IF (FIRST) THEN
         CALL IZERO_(IUNTAB,99)
         FIRST = .FALSE.
      END IF
C
C
#ifdef SYS_LINUX
C     Just return if user has asked explicitly for stderr, stdinp, or stdout
C     /HJAaJ Sep 2007
      IF (LUNIT .EQ. 0) GO TO 9002
      IF (LUNIT .EQ. 0) GO TO 8000
      IF (LUNIT .EQ. 5) GO TO 8000
      IF (LUNIT .EQ. 6) GO TO 8000
#endif 
C
C     Copy input strings to local strings,
C     remove trailing blanks, check for use of defaults
C
      FILELN = LNBLNK(FILEIN)
      IF (FILELN .GT. 0) THEN
         FILENM(1:FILELN) = FILEIN(1:FILELN)
      ELSE
         FILENM = ' '
      END IF

      ACCELN = LNBLNK(ACCEIN)
      IF (ACCELN .GT. 0) THEN
         ACCESF(1:ACCELN) = ACCEIN(1:ACCELN)
         CALL UPCASE_(ACCESF)
      ELSE
         ACCESF = ' '
      END IF

      STATLN = LNBLNK(STATIN)
      IF (STATLN .GT. 0) THEN
         STATUF(1:STATLN) = STATIN(1:STATLN)
         CALL UPCASE_(STATUF)
      ELSE
         STATLN = 7
         STATUF(1:7) = 'UNKNOWN'
      END IF

      FORMLN = LNBLNK(FORMIN)
      IF (FORMLN .GT. 0) THEN
         FORMF (1:FORMLN) = FORMIN(1:FORMLN)
         CALL UPCASE_(FORMF)
      ELSE
         FORMLN = 11
         FORMF(1:11)='UNFORMATTED'
      END IF
C
C     We then deal with the unit number
C
      IF ((LUNIT .LT. 0) .OR. (LUNIT .GT. 99)) THEN
C
C     Unit number left unassigned, we get to decide!
C
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF (IUN .GT. 99) GOTO 9001
         IF ((IUN .EQ. 5) .OR. (IUN .EQ. 6))
     &        GOTO 10
         IF (IUNTAB(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
      ELSE
C
C     The user has requested a specific unit number. We don't 
C     quite trust the user, so we check that
C     1) It is not unit 5 or 6
C     2) The file either has been closed with status='KEEP' or is
C        not currently in use
C
         IF ((LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9002
#if defined (SYS_HPUX)
C
C     This avoids us writing to unit 7 of HP-systems (which corresponds to
C     standard error on this machine)
C
         IF ((LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6) .OR. (LUNIT .EQ. 7))
     &      GOTO 9002
#else
         IF ((LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9002
#endif
         IF (IUNTAB(LUNIT) .EQ. 1) GOTO 9003
      END IF
      IF (MYNUM .EQ. 0) THEN
      IF (FILENM(1:FILELN) .EQ. 'DALTON.OUT') THEN
         OUTFIL = ' '
#if defined (SYS_T3D) || defined (SYS_T90)
         CALL PXFGETENV ('OUTFIL',6,OUTFIL,LOUTFL,IERR)
#else 
         CALL GETENV ('OUTFIL',OUTFIL)
#endif
         LENOUT = LNBLNK(OUTFIL)
         if (LENOUT .EQ. max_filename_length) then
            write (lupri,*) 'max_filename_length is too short'
            write (lupri,*) OUTFIL
            STOP
         end if
         if (LENOUT .gt. 0) then
c        IF (OUTFIL(1:1) .NE. ' ') THEN
c           LENOUT = 0
c           DO I = 1, max_filename_length
c              IF (OUTFIL(I:I) .EQ. ' ') GO TO 43
c              LENOUT = LENOUT + 1
c           END DO
c43         CONTINUE
#if defined (SYS_T3D) || defined (SYS_T90)
            CALL PXFGETENV ('WRKDIR',6,WRKDIR,LWRKDR,IERR)
#else 
            CALL GETENV ('WRKDIR',WRKDIR)
#endif
            LENWRK = lnblnk(WRKDIR)
            if (LENWRK .EQ. max_filename_length) then
               write (lupri,*) 'max_filename_length is too short'
               write (lupri,*) WRKDIR
               STOP
            end if
c           LENWRK = 0
c           DO I = 1, max_filename_length
c              IF (WRKDIR(I:I) .EQ. ' ') GO TO 44
c              LENWRK = LENWRK + 1
c           END DO
c44         CONTINUE
            FILELN = LENWRK + LENOUT + 1
            if (FILELN .EQ. max_filename_length) then
               write (lupri,*) 'max_filename_length is too short'
               write (lupri,*) WRKDIR(1:LENWRK)//'/'//
     &                         OUTFIL(1:LENOUT)
               STOP
            end if
            FILENM(1:FILELN) = WRKDIR(1:LENWRK)//'/'//
     &                         OUTFIL(1:LENOUT)
         END IF
      END IF
      END IF
C
C     We've got a file number now, get rid of the direct access case
C
      IF (ACCESF(1:6) .EQ. 'DIRECT') THEN
         IF (FILELN .LE. 1) THEN
            IF (MYNUM .EQ. 0) THEN
               FILELN = 9
               FILENM(1:9) = 'UNIT'//CHRNOS(LUNIT/10)//
     &                        CHRNOS(MOD(LUNIT,10))//'.DA'
            ELSE
               FILELN = 13
               FILENM(1:9) = 'UNIT'//CHRNOS(LUNIT/10)//
     &                        CHRNOS(MOD(LUNIT,10))//'.DA'//
     &                        '#'//CHRNOS(MYNUM/100)//
     &                        CHRNOS(MOD(MYNUM,100)/10)//
     &                        CHRNOS(MOD(MOD(MYNUM,100),10))
            END IF
         ELSE
            IF (MYNUM .GT. 0) THEN
             IF (FILENM(FILELN-3:FILELN-3) .NE. '#') THEN
               FILENM(FILELN+1:FILELN+4) = '#'//CHRNOS(MYNUM/100)//
     &                        CHRNOS(MOD(MYNUM,100)/10)//
     &                        CHRNOS(MOD(MOD(MYNUM,100),10))
               FILELN = FILELN + 4
             END IF
            END IF
         END IF
         IF (IERPCM.GT.1) WRITE (LVERR,*) 'GPOPEN_ calling OPENDX_ ',
     &      LUNIT,' ',FILENM(1:FILELN),NELEM,' ',STATUF(1:STATLN)
         CALL OPENDX_(LUNIT,FILENM(1:FILELN),NELEM,STATUF(1:STATLN),
     &               OLDDX)
      ELSE
         IF (STATUF(1:3) .EQ. 'OLD' .AND. IUNTAB(LUNIT) .EQ. 0) THEN
C
C     This better be a file with a name, and it better exist
C
            IF (FILELN .LE. 0) GOTO 9005
C
            IF (MYNUM .GT. 0 .AND. .NOT. OLDDX) THEN
             IF (INDEX(FILENM,'#') .NE. FILELN-3 .AND.
     &           INDEX(FILENM,'NUMCAL') .EQ. 0) THEN
C
C     The last test will not be necessary anymore when we get
C     locking RMA operations in MPI2
C
               FILENM(FILELN+1:FILELN+4) = '#'//CHRNOS(MYNUM/100)//
     &                        CHRNOS(MOD(MYNUM,100)/10)//
     &                        CHRNOS(MOD(MOD(MYNUM,100),10))
               FILELN = FILELN + 4
             END IF
            END IF
C
            INQUIRE(FILE=FILENM(1:FILELN),EXIST=FEXST,IOSTAT=IOS)
            IF (.NOT. FEXST) GOTO 9007
            IF (IERPCM.gt.1) write (LVERR,*) 'GPOPEN_ seq OLD ',LUNIT,
     &         ' ',FILENM(1:FILELN),' ',FORMF(1:FORMLN)
            OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),STATUS='OLD',
     &           ACCESS='SEQUENTIAL',FORM=FORMF(1:FORMLN),ERR=9004,
     &           IOSTAT=IOS)
         ELSE
            IF (STATUF(1:3) .EQ. 'NEW' .AND. IUNTAB(LUNIT) .NE. 0)
     &           GOTO 9006
            IF (FILELN .LE. 0) THEN
               IF (MYNUM .EQ. 0) THEN
                  FILELN = 6
                  FILENM = 'UNIT'//CHRNOS(LUNIT/10)//
     &                     CHRNOS(MOD(LUNIT,10))
               ELSE
                  FILELN = 10
                  FILENM = 'UNIT'//CHRNOS(LUNIT/10)//
     &                      CHRNOS(MOD(LUNIT,10))//
     &                      '#'//CHRNOS(MYNUM/100)//
     &                      CHRNOS(MOD(MYNUM,100)/10)//
     &                      CHRNOS(MOD(MOD(MYNUM,100),10))
               END IF
            ELSE
               IF (MYNUM .GT. 0) THEN
                IF (INDEX(FILENM,'#') .NE. FILELN-3) THEN
                  FILENM(FILELN+1:FILELN+4) = '#'//CHRNOS(MYNUM/100)//
     &                                     CHRNOS(MOD(MYNUM,100)/10)//
     &                                 CHRNOS(MOD(MOD(MYNUM,100),10))
                  FILELN = FILELN + 4
                END IF
               END IF
            END IF
            IF (STATUF(1:7) .EQ. 'SCRATCH') THEN
               IF (IERPCM.gt.1) write (LVERR,*) 'GPOPEN_ seq scratch ',
     &            LUNIT,' ',FORMF(1:FORMLN)
               OPEN(UNIT=LUNIT,STATUS='SCRATCH',
     &              ACCESS='SEQUENTIAL',FORM=FORMF(1:FORMLN),
     &              ERR=9004,IOSTAT=IOS)
            ELSE
               IF (IERPCM.gt.1) write (LVERR,*) 'GPOPEN_ seq ',
     &            STATUF(1:STATLN),LUNIT,' ',FILENM(1:FILELN),
     &            ' ',FORMF(1:FORMLN)
               IF(STATUF(1:STATLN) .EQ. 'NEW') THEN
Chjaaj/may2000:... if filename already exists the program will abort
C                  thus we inquire first and delete if necessary.
C                  This will often be the case if we restart a
C                  calculation.
                  INQUIRE(FILE=FILENM(1:FILELN),EXIST=FEXST,IOSTAT=IOS)
                  IF (FEXST) THEN
                     IF (IERPCM.gt.1) write (LVERR,*)
     &                  'GPOPEN_ INFO: file ',FILENM(1:FILELN),
     &                  ' already exists and is deleted first.'
                     OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),
     &                   STATUS='OLD',ACCESS='SEQUENTIAL',
     &                   FORM=FORMF(1:FORMLN),ERR=9004,IOSTAT=IOS)
                     CLOSE(UNIT=LUNIT,STATUS='DELETE')
                  END IF
               END IF
               OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),
     &              STATUS=STATUF(1:STATLN),ACCESS='SEQUENTIAL',
     &              FORM=FORMF(1:FORMLN),ERR=9004,IOSTAT=IOS)
            END IF
         END IF
      END IF
 8000 IUNTAB(LUNIT) = 1
Ckm      CALL QEXIT('GPOPEN_')
      RETURN
C     
C error branches
C
 9001 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A/A//A)')
     &   '--> ERROR (GPOPEN_) NO MORE AVAILABLE FILENUMBERS!',
     &   '--> THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &   'FILES OR', 
     &   '--> SOMEBODY HAS FORGOTTEN TO CLOSE FILES IN THE SOURCE CODE',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
 9002 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A,I3/A//2A//A)')
     &   '--> ERROR (GPOPEN_) TRYING TO OPEN AN ILLEGAL FILE NUMBER',
     &   '--> SOMEBODY HAS TRIED TO OPEN UNIT',LUNIT,
     &   '--> THE PROGRAM DOES NOT ALLOW THE USE OF THIS RESERVED '//
     &   'UNIT NUMBER',
     &   '--> Name of offending file (if any): ',FILENM(1:FILELN),
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
 9003 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A//A,I5/2A//A)')
     &   '--> ERROR (GPOPEN_) TRYING TO OPEN A FILE ALREADY IN USE',
     &   '--> SOMEBODY IS TRYING TO USE A FILENUMBER THAT IS '//
     &   'ALREADY IN USE',
     &   '--> Offending UNIT number: ',LUNIT,
     &   '--> Name of offending file (if any): ', FILENM(1:FILELN),
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
 9004 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A,I3/2A/A,I7)')
     &   '--> ERROR (GPOPEN_) UPON TRYING TO OPEN FILE ON UNIT',LUNIT,
     &   '--> with filename ',FILENM(1:FILELN),
     &   '--> IOSTAT ERROR CODE RETURNED ',IOS
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
 9005 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A/A//A)')
     &   '--> ERROR (GPOPEN_) TRYING TO OPEN A NON-EXISTING OLD FILE',
     &   '--> A FILE HAS BEEN SPECIFIED TO BE OLD, BUT THE TABLE',
     &   '--> ENTRY INDICATES THAT IT DOES NOT EXIST',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
 9006 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A/A//A,I5/2A//A)')
     &   '--> ERROR (GPOPEN_) TRYING TO OPEN AN EXISTING NEW FILE',
     &   '--> A FILE HAS BEEN SPECIFIED TO BE NEW, BUT THE TABLE',
     &   '--> ENTRY INDICATES THAT IT ALREADY EXISTS',
     &   '--> Offending UNIT number: ',LUNIT,
     &   '--> Name of offending file (if any): ',FILENM(1:FILELN),
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
 9007 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A/A//A,I5/2A//A)')
     &   '--> ERROR (GPOPEN_) TRYING TO OPEN A NON-EXISTING OLD FILE',
     &   '--> A FILE HAS BEEN SPECIFIED TO BE OLD,',
     &   '--> BUT THE FILE DOES NOT EXIST',
     &   '--> Offending UNIT number: ',LUNIT,
     &   '--> Name of offending file (if any): ',FILENM(1:FILELN),
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
C end of GPOPEN_
C
      END
C  /* Deck gpclose */
      SUBROUTINE GPCLOSE_(LUNIT,DISP)
C
C 23-Feb-2000 K.Ruud
C
C Purpose:
C   General purpose routine for closing files in the Dalton program.
C   The routine will ensure that files that are closed and that will not
C   be reopened in a later part of the program will make their unit 
C   number available for reuse.
C
C   These files are strongly machine dependent, although care has been 
C   taken to avoid using unit numbers that are illegal or reserved on a
C   particular architecture. However, the use of GPOPEN_ and GPCLSE should
C   remove much of the machine dependence in the rest of the Dalton 
C   program
C
C Input:
C  LUNIT    Mandatory
C  STATUS   Indicate whether the file should be removed or kept
C           ('KEEP' or 'DELETE'). Note that the default is that the file
C           is to be deleted
C
#include <pcm_implicit.h>
#include <pcm_priunit.h>
      CHARACTER*(*) DISP
      CHARACTER*80 FNNAME, FNNM2
      CHARACTER*3 IODIR
      LOGICAL FEXIST, OPND, OLDDX
      COMMON /PCM_UNITAR/ IUNTAB(1:99)
#include <pcm_chrnos.h>
C
C     We first deal with the unit number
C
      IF ((LUNIT .LT. 1) .OR. (LUNIT .GT. 99) .OR.
     &    (LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9001
ckr     &    (LUNIT .EQ. 5)) GOTO 9001
ckrChj  &    (LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9001
ckr      IF (LUNIT .EQ. 6) GOTO 30
C
C     Check that the file actually have been opened
C
      IF (IUNTAB(LUNIT) .EQ. 0) GOTO 9002
C
      INQUIRE(UNIT=LUNIT,EXIST=FEXIST,OPENED=OPND)
      IF (.NOT. FEXIST .OR. .NOT. OPND) THEN
         WRITE (LVPRI,'(/A/A,I3)')
     &     ' WARNING: Tried to close a non-existent or already '//
     &     'closed file', ' Unit number was: ',LUNIT
         CALL FLSHFO_(LVPRI)
      END IF
      IF (DISP .EQ. 'KEEP') THEN
         IF (IERPCM.gt.1) write (LVERR,*) 'GPCLOSE_ keep',LUNIT
         CLOSE (LUNIT,STATUS='KEEP',ERR=9003)
         IUNTAB(LUNIT) = 0
         LUNIT = -30000-LUNIT
         GOTO 30
      ELSE
#if defined (VAR_SPLITFILES)
         INQUIRE(UNIT=LUNIT,DIRECT=IODIR)
         INQUIRE(UNIT=LUNIT,NAME=FNNAME)
         IF (IODIR .EQ. 'YES') THEN
            INQUIRE(UNIT=LUNIT,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
            LRECL = LRECL/8
#endif
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
            LRECL = LRECL/4
#endif
         END IF
         LN = 1
 15      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 15
         END IF
         LN = LN - 1
         IF (IODIR .NE. 'YES') THEN
            IF (FNNAME((LN-1):(LN-1)) .EQ. '-') THEN
               CLOSE(LUNIT,STATUS='KEEP')
               LN = LN - 2
               OPEN (UNIT=LUNIT,FILE=FNNAME(1:LN),STATUS='OLD')
            END IF
         END IF
         IF (IERPCM.gt.1) write (LVERR,*) 'GPCLOSE_ delete',LUNIT
         CLOSE (LUNIT,STATUS='DELETE',ERR=20)
         IUNIT = 0
 16      CONTINUE
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(IUNIT)
         LN2 = LN + 2
         INQUIRE(FILE=FNNM2(1:LN2),EXIST=FEXIST)
         IF (FEXIST) THEN
            IF (IODIR .EQ. 'YES') THEN
               CALL OPENDX_(LUNIT,FNNM2(1:LN2),LRECL,'OLD',OLDDX)
            ELSE
               OPEN(UNIT=LUNIT,FILE=FNNM2(1:LN2),STATUS='OLD')
            END IF
            CLOSE(LUNIT,STATUS='DELETE')
            IUNIT = IUNIT + 1
            GOTO 16
         END IF
#else
         CLOSE (LUNIT,STATUS='DELETE',ERR=20)
#endif
         IUNTAB(LUNIT) = 0
         LUNIT = -20000-LUNIT
         GOTO 30
      END IF
C
C     We just close it
C     
 20   CLOSE (LUNIT)
      IUNTAB(LUNIT) = 0
      LUNIT = -30000-LUNIT
C
C     We release that saved unit number by resetting it to -(code)-LUNIT
C
 30   CONTINUE
      RETURN
C     
C error branches
C
 9001 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A,I15/A//A)')
     &   '--> ERROR (GPCLOSE_) TRYING TO CLOSE AN ILLEGAL FILE NUMBER',
     &   '--> SOMEBODY HAS TRIED TO CLOSE UNIT',LUNIT,
     &   '--> THE PROGRAM DOES NOT ALLOW THE USE OF THIS UNIT NUMBER',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
 9002 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A,I15//A)')
     &   '--> ERROR (GPCLOSE_) TRYING TO CLOSE A FILE NOT IN USE',
     &   '--> SOMEBODY IS TRYING TO USE A FILENUMBER THAT HAS '//
     &   'NOT BEEN USED :',LUNIT,
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
 9003 CONTINUE
      IF (LVPRI .GE. 0) THEN
         LUGPPRI = LVPRI
      ELSE
         LUGPPRI = 6
      END IF
      WRITE (LUGPPRI,'(//A/A/A//A)')
     &   '--> ERROR (GPCLOSE_) TRYING TO KEEP A SCRATCH FILE',
     &   '--> A FILE HAS BEEN INDICATED TO BE CLOSED AND KEPT,',
     &   '--> BUT IT APPEARS THE FILE IS A SCRATCH FILE',
     &   '### Please report the problem to dalton-admin@kjemi.uio.no'
Ckm      CALL QTRACE_(LUGPPRI)
      STOP
C
C end of GPCLOSE_
C
      END
      SUBROUTINE UPCASE_(WORD)
#include "pcm_implicit.h"
C                                                                               
C     Routine that uppercases the word provided. Allows the input easily to     
C     be made case insensitive. F90 routine supplied by Anthony Stone and       
C     77-ified by K.Ruud                                                        
C                                                                               
      CHARACTER word*(*)
      CHARACTER*26 UPPER, LOWER
C                                                                               
      UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      LOWER = "abcdefghijklmnopqrstuvwxyz"
      DO I = 1, LEN(WORD)
         K = INDEX(LOWER,WORD(I:I))
         IF (K .NE. 0) WORD(I:I)=UPPER(K:K)
      END DO
C                                                                               
      RETURN
      END
      SUBROUTINE FLSHFO_ (IUNIT)
C                                                                               
C *** THIS SUBROUTINE IS SYSTEM DEPENDENT ***                                   
C                                                                               
C     Flush formatted output unit (empty buffers).                              
C     If no flush utility, this is achieved by                                  
C     CLose and reOPen Formatted Output                                         
C                                                                               
C Written 21-Nov-1983 by Hans Jorgen Aa. Jensen in Uppsala, Sweden.             
C Last revision 16-Jul-1984 hjaaj / 30-Oct-1984 hjaaj (extendsize)              
C 10-Feb-1989 hjaaj, renamed CLOPFO to FLSHFO_                                   
C                                                                               
C Calls to this subroutine makes it possible to read the output                 
C up to the moment of the last call while the program continues                 
C executing (provided the computer allows shared access).                       
C This subroutine may be a dummy routine.                                       
C                                                                               
      IF (IUNIT .GE. 0) THEN
C     ... do not try to flush unused units (e.g. LVPRI on a slave) /hjaaj       
#if  defined (SYS_AIX)  || defined (SYS_IRIX) || defined (SYS_NEC)  \
  || defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_LINUX) \
  || defined (SYS_SUN)  || defined (SYS_HAL) || defined (SYS_T90)   \
  || defined (SYS_HPUX) || defined (SYS_SX)
C                                                                               
C        Force transfer of all buffered output to the file or device            
C        associated with logical unit IUNIT.                                    
C                                                                               
         CALL FLUSH(IUNIT)
#endif
      END IF
      RETURN
      END
      SUBROUTINE OPENDX_ (LUDX,NAME,NELEM,STATUS,OLDDX)
C
C 15-Jun-1985 hjaaj
C
C Revisions :  9-Dec-1987 hjaaj (Alliant version)
C
C Purpose:
C   Open files for direct access through WRITDX and READDX routines.
C   The ....DX routines enables direct access, even when the number
C   of elements per record (the logical record length) is greater
C   than the maximum physical record length.
C   >>> THIS IS MACHINE DEPENDENT <<<
C
C Input:
C  LUDX     file unit number
C  NELEM    number of integer words per logical record
C  STATUS   file status: 'OLD', 'NEW', or 'UNKNOWN'
C
C Output:
C  OLDDX    logical, true if old LUDX file was opened
C
C
      CHARACTER*(*) NAME, STATUS
      LOGICAL OLDDX, FEXIST
#include <pcm_priunit.h>
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
C
C     CRAY has 8 byte integers.
C
      LRECL  = 8*NELEM
C
      IF (STATUS .EQ. 'NEW') THEN
         INQUIRE(FILE=NAME,EXIST=FEXIST)
         IF (FEXIST) THEN
            OPEN(UNIT=LUDX,STATUS='OLD',FILE=NAME)
            CLOSE(UNIT=LUDX,STATUS='DELETE')
         END IF
         GO TO 300
      END IF
C
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',ERR=301,IOSTAT=IOS,RECL=LRECL)
         OLDDX = .FALSE.
C
C     ERROR HANDLING (JAA 1997-06-06)
C
 301     IF (IOS .GT. 0) THEN
            WRITE(102,*) '--> ERROR (OPENDX_) NEW FILE OPEN FAILED'
            WRITE(102,*) '    FILE=',NAME,', UNIT=',LUDX,
     $           ', IOSTAT=',IOS
            IF (IOS .LT. 1000) THEN
               WRITE(102,'(A,I3.3,A)') '    (try explain sys-',IOS,')'
            ELSE IF (IOS .LT. 6000) THEN
               WRITE(102,'(A,I4,A)') '    (try explain lib-',IOS,')'
            ENDIF
            CALL FLUSH(102,IOS)
            CALL TRBK()
#if defined(SYS_T3D)
            CALL STOP_ALL()
#else
            STOP
#endif
         ENDIF
  600 CONTINUE
#endif
#if defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_LINUX) || defined (SYS_SUN) || defined (SYS_HAL)
C
C     Assume same for IBM-AIX--  1-Oct-1990 hjaaj
C     Assume same for HP-UX  -- 21-Aug-1991 hjaaj
C     Assume same for DEC    -- 21-May-1992 hjaaj
C     Assume same for IRIX   --  3-Feb-1994 hjaaj; except recl in *4 words
C     Assume same for PARAGON-- 13-Oct-1994 hjaaj
C
#if (defined (SYS_IRIX) && defined(NO_F90)) || defined (SYS_NEC)
      LRECL  = NELEM
#else
#ifdef VAR_INT64
      LRECL  = 8*NELEM
#else
      LRECL  = 4*NELEM
#endif
#endif
C
      IF (STATUS .EQ. 'NEW') THEN
C        some systems abort if we try to open as NEW a file,
C        which already exist
         INQUIRE(FILE=NAME,EXIST=FEXIST)
         IF (FEXIST) THEN
            OPEN(UNIT=LUDX,STATUS='OLD',FILE=NAME)
            CLOSE(UNIT=LUDX,STATUS='DELETE')
         END IF
         GO TO 300
      END IF
C
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .FALSE.
  600  CONTINUE
#endif
#if !defined (SYS_CRAY) && !defined (SYS_AIX) && !defined (SYS_HPUX) && !defined (SYS_DEC) && !defined (SYS_IRIX) && !defined (SYS_NEC) && !defined (SYS_T3D) && !defined (SYS_LINUX) && !defined (SYS_SUN) && !defined (SYS_HAL) && !defined (SYS_T90)
 >>>>> insert appropriate OPEN statements in OPENDX_.
#endif
      RETURN
C
C error branches
C
 9000 CONTINUE
      WRITE (LVPRI,'(//2A/2A/A,I5)')
     *   ' --> ERROR (OPENDX_) INVALID STATUS KEYWORD: ',STATUS,
     *   '                    FILE NAME   =',NAME,
     *   '                    FILE NUMBER =',LUDX
Ckm      CALL QTRACE_(LVPRI)
      STOP
C
 9100 CONTINUE
      WRITE (LVPRI,'(//A/2A/A,I5/A)')
     *   ' --> ERROR (OPENDX_) OLD FILE NOT FOUND',
     *   '                    FILE NAME   =',NAME,
     *   '                    FILE NUMBER =',LUDX,
     *   ' --- or wrong record length on old file.'
Ckm      CALL QTRACE_(LVPRI)
      STOP
C
C end of OPENDX_
C
      END
      SUBROUTINE STOPIT_(SUB,PLACE,INT1,INT2)
#include <pcm_implicit.h>
#include <pcm_priunit.h>
      CHARACTER*(*) SUB, PLACE
      WRITE (LVPRI,'(//3A)') ' Work space exceeded in subroutine ',
     *                         SUB,'.'
      IF (LNBLNK(PLACE) .GT. 0) THEN
         WRITE (LVPRI,'(/2A)') ' Location: ',PLACE
      END IF
      JNT1  = ABS(INT1)
      JNT2  = ABS(INT2)
      LWRKR = MAX(JNT1,JNT2)
      LWRKA = MIN(JNT1,JNT2)
      IF (INT1 .LT. 0 .OR. INT2 .LT. 0) THEN
         WRITE (LVPRI,'(/A,I10)  ') ' Space required  >',LWRKR
      ELSE
         WRITE (LVPRI,'(/A,I10)  ') ' Space required  =',LWRKR
      END IF
      WRITE (LVPRI,'( A,I10)') ' Space available =',LWRKA
      STOP
      END
      SUBROUTINE READT_ (IT,N,X)
#include "pcm_implicit.h"
#include "pcm_priunit.h"
      CHARACTER*60 FNAME
      DIMENSION X(N)
      IF (IT .LE. 0) GOTO 30
      READ (IT,END=10,ERR=20) X
      RETURN
 10      CONTINUE
      INQUIRE(UNIT=IT,NAME=FNAME)
      WRITE (LVPRI,*) 'READT_: END reading file ',FNAME,
     & ' UNIT ',IT, ' record length ',N
      STOP
 20      CONTINUE
      INQUIRE(UNIT=IT,NAME=FNAME)
      WRITE (LVPRI,*) 'READT_: ERROR reading file ',FNAME,
     & ' UNIT ',IT, ' record length ',N
      STOP
 30      CONTINUE
      WRITE (LVPRI,*) 'READT_ ERRROR: non-positive file unit number ',IT
      STOP
      END
      SUBROUTINE ERRWRK_ (STRING,LNEED,LAVAIL)
C                                                                               
C Version 6-Mar-1985 by hjaaj                                                   
C                                                                               
      CHARACTER*(*) STRING
C                                                                               
#include <pcm_priunit.h>
C                                                                               
Ckm      CALL QENTER_('ERRWRK_')
      IF (LNEED .GE. 0) THEN
         WRITE (LVPRI,1010) STRING,LNEED,LAVAIL
      ELSE
         WRITE (LVPRI,1020) STRING,-LNEED,LAVAIL
      END IF
Ckm      CALL QTRACE_(LVPRI)
      STOP
C                                                                               
 1010  FORMAT(/'  FATAL ERROR, insufficient core for ',A,
     *       /T16,'( Need:',I10,', available:',I10,' )')
 1020   FORMAT(/'  FATAL ERROR, insufficient core for ',A,
     *      //T16,'Need      :',I10,' + uncalculated amount',
     *       /T16,'Available :',I10)
      END
      SUBROUTINE GETTIM_(TIMCPU,TIMWAL)
C
C Written 18-Dec-1984 Hans Joergen Aa. Jensen
C
C Return elapsed CPU time and elapsed real time.
C
#include "pcm_implicit.h"
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
C...  Cray unicos timing code
      REAL    TCPU0, TWALL0
      LOGICAL FIRST
      PARAMETER ( DP001 = 0.001 )
      SAVE TCPU0, TWALL0, FIRST
      DATA FIRST /.TRUE./
C
      IF (FIRST) THEN
         FIRST  = .FALSE.
#if defined (SYS_CRAY) || defined (SYS_T90)
         TCPU0  = SECOND()
         TWALL0 = TIMEF()
      ENDIF
      TIMCPU = SECOND() - TCPU0
#else
         TCPU0 = TSECND()
         TWALL0 = TIMEF()
      ENDIF
      TIMCPU = TSECND() - TCPU0
#endif
      TIMWAL = DP001 * ( TIMEF() - TWALL0 )
#else
#if defined (SYS_NEC)
C...  NEC SX4     timing code
      REAL*8  TCPU0, TWALL0
      REAL*8  SECOND1,SECOND2
      LOGICAL FIRST
      PARAMETER ( DP001 = 0.001 )
      SAVE TCPU0, TWALL0, FIRST
      DATA FIRST /.TRUE./
C
      IF (FIRST) THEN
         FIRST  = .FALSE.
         call clock(second1)
         TCPU0  = SECOND1
         call etime(second2)
         TWALL0 = second2
      ENDIF

      call clock(second1)
      call etime(second2)
      TIMCPU = second1 - TCPU0
      TIMWAL = second2 - TWALL0
#else
#if defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_IRIX) \
 || defined (SYS_SUN) || defined (SYS_HAL) \
 || (defined (SYS_LINUX) && defined (VAR_PGF77)) || defined(VAR_ABSOFT) \
 || defined (SYS_DEC) || defined (SYS_SX)
#if defined (SYS_IRIX) || defined (SYS_SUN) || defined (SYS_DEC) \
 || defined (SYS_SX)  || defined (SYS_AIX)
C...  SGI IRIX  etc. code
      external time
C     IRIX: otherwise "subroutine time", an intrinsic fu. will be used
#endif
      integer time
      real*4  etime, tarray(2)
      logical first
      save    first,iwall0,TCPU0
      data    first /.true./
C
      if (first) then
         first = .false.
         TCPU0 = etime(tarray)
         iwall0 = time()
      end if
      TIMCPU = etime(tarray) - TCPU0
      itimwl = time() - iwall0
      TIMWAL = itimwl
#else
#if (defined (SYS_LINUX) && defined (VAR_IFC))
      integer time,darr(9)
      real*4  etime, tarray(2)
      logical first
      save    first,iwall0,TCPU0
      data    first /.true./
C
      if (first) then
         first = .false.
         TCPU0 = etime(tarray)
          call gettimeofday(darr,err)
          iwall0 = darr(1)
      end if
      TIMCPU = etime(tarray) - TCPU0
      call gettimeofday(darr,err)
      itimwl = darr(1) - iwall0
      TIMWAL = itimwl
#else
#if defined (VAR_G77) || defined (VAR_GFORTRAN)
C 990819-ekd: Use time8() and second() system routines for timings
C             ONLY available for GNU g77 compiler!!!!!!!!
      logical first
      save    TWALL0,TCPU0
      data    first /.true./
C
      if (first) then
         first = .false.
         TWALL0 = time8()
         TCPU0  = second()
      end if
      TIMWAL = time8() - TWALL0
      TIMCPU = second() - TCPU0
#else
      PARAMETER (D0 = 0.0D0)
      logical first
      save    TWALL0,TCPU0
      data    first /.true./
C
      if (first) then
         first = .false.
         TCPU0  = second()
      end if
      TIMCPU = SECOND() - TCPU0
C insert appropriate routine to get elapsed real time (wall time)
C here
      TIMWAL = D0
#endif  /* VAR_G77, VAR_GFORTRAN */
#endif  /* VAR_IFC */
#endif  /* SYS_AIX, SYS_IRIX, SYS_SUN, SYS_HPUX etc. */
#endif  /* SYS_NEC */
#endif  /* SYS_CRAY etc. */
C
      RETURN
      END
      SUBROUTINE DSPTSI_(N,ASP,ASI)
C                                                                               
C  8-Feb-1987 Hans Joergen Aa. Jensen                                           
C  900720-hjaaj: inserted NOSYNC directive.                                     
C                                                                               
C Purpose: Transform from SP format to SI format, that is:                      
C          unpack symmetric,   packed matrix ASP                                
C              to symmetric, unpacked matrix ASI.                               
C                                                                               
#include <pcm_implicit.h>
      DIMENSION ASP(*), ASI(N,*)
C                                                                               
      ENTRY      DSPTGE_(N,ASP,ASI)
C     ... equivalent subroutine name                                            
C                                                                               
#include <pcm_nosync.h>
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
#include <pcm_nosync.h>
         DO 100 I = 1,J-1
            ASI(I,J) = ASP(JOFF+I)
            ASI(J,I) = ASP(JOFF+I)
 100            CONTINUE
         ASI(J,J) = ASP(JOFF+J)
 200      CONTINUE
C                                                                               
      RETURN
      END
      SUBROUTINE DGETRN_(AGE,NROWA,NRDIMA)
C                                                                               
C Replace AGE by AGE(transposed)                                                
C                                                                               
C 3-Apr-1987 HJAaJ                                                              
C 900108-hjaaj: block with NBLK for reduced paging                              
C   when virtual memory                                                         
C new name 971103-hjaaj (old name DGETRS was same as                            
C   a routine in LAPACK for solving linear equations;                           
C   when linking with complib on SGI/IRIX the LAPACK routine                    
C   was loaded instead of this one).                                            
C                                                                               
#include <pcm_implicit.h>
      DIMENSION AGE(NRDIMA,*)
      PARAMETER (NBLK = 128)
      DO 400 JBLK = 1,NROWA,NBLK
         JEND = MIN(NROWA,JBLK-1+NBLK)
         DO 300 IBLK = 1,JBLK,NBLK
            IEND = MIN(NROWA,IBLK-1+NBLK)
            DO 200 J = JBLK,JEND
               DO 100 I = IBLK,MIN(J-1,IEND)
                  SWAP     = AGE(I,J)
                  AGE(I,J) = AGE(J,I)
                  AGE(J,I) = SWAP
 100                        CONTINUE
 200                               CONTINUE
 300                                   CONTINUE
 400                                    CONTINUE
      RETURN
      END
      SUBROUTINE DUNIT_(A,N)
C                                                                               
C  SUBROUTINE DUNIT_ SETS THE REAL SQUARE MATRIX A EQUAL                         
C  TO A UNIT MATRIX.                                                            
C  /VER 2/ 14-Sep-1983 hjaaj                                                    
C                                                                               
#include "pcm_implicit.h"
      DIMENSION A(*)
      PARAMETER (D1=1.0D00, D0=0.0D00)
C                                                                               
      NN = N*N
      DO 100 I = 1,NN
         A(I) = D0
 100      CONTINUE
      N1 = N + 1
      DO 200 I = 1,NN,N1
         A(I) = D1
 200      CONTINUE
      RETURN
      END
      SUBROUTINE DSPSOL_ (N,NSIM,AP,B,KPVT,INFO)
C                                                                               
C Written 8-Feb-1985 Hans Joergen Aa. Jensen                                    
C No revisions.                                                                 
C                                                                               
C Purpose:                                                                      
C  Solve the NSIM simultaneous eqautions:                                       
C                                                                               
C     B(n,nsim) := A(n,n) inverse * B(n,nsim)                                   
C                                                                               
C  AP is A in packed format.                                                    
C  KPVT is a scratch array of length N.                                         
C                                                                               
#include <pcm_implicit.h>
      DIMENSION AP(*),B(N,*),KPVT(*)
C                                                                               
      CALL DSPFA_ (AP,N,KPVT,INFO)
      IF (INFO.NE.0) RETURN
C                                                                               
      DO 100 J = 1,NSIM
        CALL DSPSL_ (AP,N,KPVT,B(1,J))
 100     CONTINUE
C                                                                               
      RETURN
      END
      SUBROUTINE DSPFA_(AP,N,KPVT,INFO)
      INTEGER N,KPVT(*),INFO
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL AP(*)
#else
      DOUBLE PRECISION AP(*)
#endif
C
C     DSPFA_ FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN
C     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW DSPFA_ BY DSPSL_.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPFA_ BY DSPSL_.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPFA_ BY DSPDI.
C     TO COMPUTE  INERTIA(A) , FOLLOW DSPFA_ BY DSPDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DSPFA_ BY DSPDI.
C
C     ON ENTRY
C
C        AP      DOUBLE PRECISION (N*(N+1)/2)
C                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER_ OF THE MATRIX  A .
C
C     OUTPUT_
C
C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
C                     BUT IT DOES INDICATE THAT DSPSL_ OR DSPDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A SYMMETRIC MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K)  = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY_,DSWAP_,IDAMAX_
C     FORTRAN DABS,DMAX1,DSQRT
C
C     INTER_NAL VARIABLES
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL ABSAKK,ALPHA,COLMAX,ROWMAX
#else
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
#endif
      INTEGER IDAMAX_,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
C
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
      IK = (N*(N - 1))/2
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (AP(1) .EQ. 0.0D0) INFO = 1
C     ......EXIT
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTER_CHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         KK = IK + K
         ABSAKK = DABS(AP(KK))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = IDAMAX_(K-1,AP(IK+1),1)
         IMK = IK + IMAX
         COLMAX = DABS(AP(IMK))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            IM = IMAX*(IMAX - 1)/2
            IMJ = IM + 2*IMAX
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,DABS(AP(IMJ)))
               IMJ = IMJ + J
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX_(IMAX-1,AP(IM+1),1)
               JMIM = JMAX + IM
               ROWMAX = DMAX1(ROWMAX,DABS(AP(JMIM)))
   50       CONTINUE
            IMIM = IMAX + IM
            IF (DABS(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTER_CHANGE.
C
               CALL DSWAP_(IMAX,AP(IM+1),1,AP(IK+1),1)
               IMJ = IK + IMAX
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  JK = IK + J
                  T = AP(JK)
                  AP(JK) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            IJ = IK - (K - 1)
            DO 130 JJ = 1, KM1
               J = K - JJ
               JK = IK + J
               MULK = -AP(JK)/AP(KK)
               T = MULK
               CALL DAXPY_(J,T,AP(IK+1),1,AP(IJ+1),1)
               IJJ = IJ + J
               AP(JK) = MULK
               IJ = IJ - (J - 1)
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            KM1K = IK + K - 1
            IKM1 = IK - (K - 1)
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTER_CHANGE.
C
               CALL DSWAP_(IMAX,AP(IM+1),1,AP(IKM1+1),1)
               IMJ = IKM1 + IMAX
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  JKM1 = IKM1 + J
                  T = AP(JKM1)
                  AP(JKM1) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  150          CONTINUE
               T = AP(KM1K)
               AP(KM1K) = AP(IMK)
               AP(IMK) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = AP(KK)/AP(KM1K)
               KM1KM1 = IKM1 + K - 1
               AKM1 = AP(KM1KM1)/AP(KM1K)
               DENOM = 1.0D0 - AK*AKM1
               IJ = IK - (K - 1) - (K - 2)
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  JK = IK + J
                  BK = AP(JK)/AP(KM1K)
                  JKM1 = IKM1 + J
                  BKM1 = AP(JKM1)/AP(KM1K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY_(J,T,AP(IK+1),1,AP(IJ+1),1)
                  T = MULKM1
                  CALL DAXPY_(J,T,AP(IKM1+1),1,AP(IJ+1),1)
                  AP(JK) = MULK
                  AP(JKM1) = MULKM1
                  IJJ = IJ + J
                  IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
      SUBROUTINE DSPSL_(AP,N,KPVT,B)
      INTEGER N,KPVT(*)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL AP(*),B(*)
#else
      DOUBLE PRECISION AP(*),B(*)
#endif
C
C     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY DSPFA_.
C
C     ON ENTRY
C
C        AP      DOUBLE PRECISION(N*(N+1)/2)
C                THE OUTPUT_ FROM DSPFA_.
C
C        N       INTEGER
C                THE ORDER_ OF THE MATRIX  A .
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DSPFA_.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF  DSPCO  HAS SET RCOND .EQ. 0.0
C        OR  DSPFA_  HAS SET INFO .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DSPFA_(AP,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DSPSL_(AP,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY_,DDOT_
C     FORTRAN IABS
C
C     INTER_NAL VARIABLES.
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL AK,AKM1,BK,BKM1,DDOT_,DENOM,TEMP
#else
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT_,DENOM,TEMP
#endif
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
      IK = (N*(N - 1))/2
   10 IF (K .EQ. 0) GO TO 80
         KK = IK + K
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTER_CHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY_(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/AP(KK)
            K = K - 1
            IK = IK - K
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IKM1 = IK - (K - 1)
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTER_CHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY_(K-2,B(K),AP(IK+1),1,B(1),1)
               CALL DAXPY_(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            KM1K = IK + K - 1
            KK = IK + K
            AK = AP(KK)/AP(KM1K)
            KM1KM1 = IKM1 + K - 1
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = B(K)/AP(KM1K)
            BKM1 = B(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            IK = IK - (K + 1) - K
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
      IK = 0
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT_(K-1,AP(IK+1),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTER_CHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            IK = IK + K
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT_(K-1,AP(IK+1),1,B(1),1)
               IKP1 = IK + K
               B(K+1) = B(K+1) + DDOT_(K-1,AP(IKP1+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTER_CHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            IK = IK + K + K + 1
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
      SUBROUTINE ORDER_(EVEC,EVAL,N,NEVEC)
C                                                                               
C Purpose: order the N values in EVAL and their associated vectors              
C          in EVEC so EVAL(i+1) .ge. EVAL(i)                                    
C                                                                               
C Revisions:                                                                    
C   29-Jul-1992 hjaaj (only dswap if nevec .gt. 0)                              
C    2-Nov-1984 hjaaj (new parameter NEVEC, EVEC(1:NEVEC,1:N))                  
C   27-Oct-1984 hjaaj (reduced number of swaps)                                 
C                                                                               
#include "pcm_implicit.h"
      DIMENSION EVEC(*),EVAL(*)
      IF (N.LE.1) RETURN
      IN = 1
      DO 10 I=1,N-1
        EMIN = EVAL(I)
        IMIN = I
        DO 20 J=I+1,N
          IF (EVAL(J) .LT. EMIN) THEN
            EMIN = EVAL(J)
            IMIN = J
          ENDIF
 20          CONTINUE
        IF (IMIN.NE.I) THEN
          EVAL(IMIN)=EVAL(I)
          EVAL(I)=EMIN
          IF (NEVEC .GT. 0) THEN
            CALL DSWAP_(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*NEVEC+1),1)
          ENDIF
        ENDIF
        IN = IN + NEVEC
 10      CONTINUE
      RETURN
      END
      SUBROUTINE JACO_ (F,V,NB,NMAX,NROWV,BIG,JBIG)
C
C     F is symmetric packed matrix of dimension NB.
C       The first block of size NMAX will be diagonalized.
C     V is for eigenvectors, only V(NROWV,NMAX) will be referenced.
C       On entry it must correspond the basis vectors corresponding to the
C       F matrix on entry, e.g. unit matrix or AO coefficients for each MO.
C
C Revisions:
C   2-Nov-1984 hjaaj (new parameter NROWV such that
C                     dim(V) = (NROWV,NMAX). This makes
C                     it possible to solve eigenproblem
C                     in a reduced basis but get the
C                     eigenvectors in the original full
C                     basis, e.g. less mo's than ao's)
C  23-Feb-1989 hjaaj  Note that if NROWV = 0 then only
C                     eigenvalues will be calculated,
C                     V matrix will not be referenced.
C  27-Jul-1990 hjaaj  Changed -CX,+SX transformation to +CX,-SX
C                     transformation; probably the -CX,+SX
C                     transformation was responsible for that
C                     the eigenvectors easily changed sign.
C                     Changed initial test on NB. Changed SD.
C                     Optimized IR loop.
C     Jun-1992 ov     Parameters for 0.5, 1.5, ... (for Cray)
C  20-Jul-1992 hjaaj  Changed C1,C2 to THRZER
C  30-oct-1992 hjaaj  zero f(ab) to avoid round-off errors
C                     absolute conv.threshold SD=C1
C  18-aug-2005 wmk    Changed C1 to 1.D-15
#include "pcm_implicit.h"
#include "pcm_priunit.h"
      DIMENSION F(*),V(NROWV,*)
      DIMENSION BIG(*) ,JBIG(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, ROOT2 = 0.707106781186548D0)
      PARAMETER(DP5 = 0.5D0, D1P5 = 1.5D0, D1P375 = 1.375D0,
     *          D3P875 = 3.875D0, DP25 = 0.25D0)
#include "pcm_thrzer.h"
      DATA C1,C2,C3,C4,C5,C6/1.D-15,THRZER,1.D-20,1.D-14,1.D-9,1.D-5/
Cwas: DATA C1,C2,C3,C4,C5,C6/THRZER,THRZER,1.D-20,1.D-14,1.D-9,1.D-5/
Cwas: DATA C1,C2,C3,C4,C5,C6/1.D-12,1.D-12,1.D-20,1.D-14,1.D-9,1.D-5/
      IF (NB.LE.1 .OR. NMAX.LE.0) RETURN
Cwas: IF (NB.EQ.1) RETURN !900727-hjaaj
Ckm      CALL QENTER_('JACO_')
      CALL GETTIM_(TSTRT, WSTRT)
      DO 190 I=1,NB
         JBIGI=0
         J=MIN(I-1,NMAX)
         IF (J .GT. 0) THEN
            II = (I*I-I)/2
            ABIGI=D0
            DO 18 K=1,J
            IF (ABIGI .GE. ABS(F(II+K))) GO TO  18
               ABIGI=ABS(F(II+K))
               JBIGI=K
   18       CONTINUE
         END IF
         IF (JBIGI .GT. 0) THEN
            JBIG(I) = JBIGI
            BIG(I)  = F(II+JBIGI)
         ELSE
            JBIG(I) = 0
            BIG(I)  = D0
         END IF
  190 CONTINUE
C
#if defined (VAR_OLDCODE)
C 900727-hjaaj:
C SD calculation was done in every Jacobi iteration.
C Now the largest absolute element in F is found once and
C the SD based on that value is used in every iteration.
  410 SD=1.05D 00
      DO 220 J=1,NMAX
         DAB=ABS(F(J*(J+1)/2))
CHJ-861103: commented out next line, it seems to make the loop
C           meaningless (setting SD equal to J=NMAX value always!)
C        IF (SD .GT. DAB) SD=DAB
  220    SD=MAX(SD,DAB)
      SD=MAX(C1,C2*SD)
#else
C 921030-hjaaj: SD = C1 now
      NNB = (NB*NB+NB)/2
C     SD = 1.05D0
C     DO 220 J = 1,NNB
C        SD = MAX(SD, ABS(F(J)) )
C 220 CONTINUE
C     SD=MAX(C1,C2*SD)
      SD=C1
C
      MXITJA = 50*NNB
      ITJACO_ = 0
  410 ITJACO_ = ITJACO + 1
      IF (ITJACO_ .GT. MXITJA) THEN
         STOP
      END IF
#endif
      T = D0
      DO 230 I=2,NB
      IF (T .GE.  ABS(BIG(I))) GO TO 230
         T = ABS(BIG(I))
         IB= I
  230 CONTINUE
      IF(T.LT.SD) GO TO 420
         IA =JBIG(IB)
         IAA=IA*(IA-1)/2
         IBB=IB*(IB-1)/2
         DIF=F(IAA+IA)-F(IBB+IB)
         IF( ABS(DIF) .GT. C3) GO TO 271
            SX=ROOT2
            CX=ROOT2
         GO TO 270
  271       T2X2 =BIG(IB)/DIF
            T2X25=T2X2*T2X2
         IF(T2X25 .GT. C4) GO TO 240
            CX=D1
            SX=T2X2
         GO TO 270
  240    IF(T2X25 .GT. C5) GO TO 250
            SX=T2X2*(D1 - D1P5*T2X25)
            CX=D1 - DP5*T2X25
         GO TO 270
  250    IF(T2X25 . GT . C6) GO TO 260
            CX=D1+T2X25*(T2X25*D1P375 - DP5 )
            SX= T2X2*(D1 + T2X25*(T2X25*D3P875 - D1P5))
         GO TO 270
  260       T = DP25  / SQRT(DP25   + T2X25)
            CX= SQRT(DP5   + T)
            SX= SIGN( SQRT(DP5 - T),T2X2)
  270    CONTINUE
#include "pcm_ivdep.h"
         DO 275 IR=1,IA
            T        = F(IAA+IR)*SX
            F(IAA+IR)= F(IAA+IR)*CX+F(IBB+IR)*SX
            F(IBB+IR)=-T           +F(IBB+IR)*CX
  275    CONTINUE
         IEAA=IAA+IA
         IEAB=IBB+IA
         TT  =F(IEAB)
         F(IEAB)=BIG(IB)
         IF (JBIG(IA) .EQ. 0) THEN
            IRST = IA   + 1
            IEAR = IEAA + IA
            IEBR = IEAB + 1
         ELSE
            IRST = IA
            IEAR = IEAA
            IEBR = IEAB
         END IF
         DO 390 IR = IRST,NB
#if !defined (VAR_OLDCODE)
            IF (IR .EQ. IA) GO TO 360
C              ... we have checked above that JBIG(IA) .ne. 0
#else
            IF (IR .EQ. IA) THEN
               GO TO 360
C              ... we have checked above that JBIG(IA) .ne. 0
C              IF(JBIG(IR)) 360,380,360
            END IF
#endif
            T      = F(IEAR)*SX
            F(IEAR)= F(IEAR)*CX+F(IEBR)*SX
            F(IEBR)=-T         +F(IEBR)*CX
            T   =F(IEAR)
            IT  =IA
  310          F(IEAA)=F(IEAA)*CX+F(IEAB)*SX
C              921030+hjaaj: zero f(ab) to avoid round-off errors
C              F(IEAB)=     TT*CX+F(IEBR)*SX
               F(IEAB)=     D0
               F(IEBR)=    -TT*SX+F(IEBR)*CX
            GO TO 360
  320       IF(ABS(T) .GE.  ABS(F(IEBR))) GO TO 340
               T   =F(IEBR)
               IT  =IB
  340       IF(ABS(T) .LT.  ABS(BIG(IR))) GO TO 350
               BIG(IR)  = T
               JBIG(IR) = IT
            GO TO 380
  350       IF(IA .NE. JBIG(IR) .AND. IB .NE. JBIG(IR))  GO TO 380
  360          K= IEAR - IA
               JBIGI = 0
               IR1=MIN (IR-1,NMAX)
               IF (IR1 .GT. 0) THEN
                  ABIGI = D0
                  DO 370 I=1,IR1
                  IF(ABIGI .GE. ABS(F(K+I)))  GO TO 370
                     ABIGI = ABS(F(K+I))
                     JBIGI =I
  370             CONTINUE
               END IF
               IF (JBIGI .GT. 0) THEN
                  JBIG(IR) = JBIGI
                  BIG(IR)  = F(K+JBIGI)
               ELSE
                  JBIG(IR) = 0
                  BIG(IR)  = D0
               END IF
  380       CONTINUE
               IEAR = IEAR + IR
               IF (IR .GE. IB) THEN
                  IEBR = IEBR + IR
               ELSE
                  IEBR = IEBR + 1
               END IF
  390       CONTINUE
C
#include "pcm_ivdep.h"
         DO I=1,NROWV
            VIIA = V(I,IA)
            VIIB = V(I,IB)
            V(I,IA) =  VIIA*CX + VIIB*SX
            V(I,IB) = -VIIA*SX + VIIB*CX
         END DO
      GO TO 410
  420 CONTINUE
c     CALL GETTIM_(TEND, WEND)
c     WRITE(LVPRI,'(/A,4I10/A,2F20.2)')
c    &    'JACO_ -- ITJACO, NB,NMAX,NROWV :',ITJACO, NB,NMAX,NROWV,
c    &    'JACO_ -- CPUTIME, WALLTIME     :',TEND-TSTRT, WEND-WSTRT
Ckm      CALL QEXIT('JACO_')
      RETURN
      END

