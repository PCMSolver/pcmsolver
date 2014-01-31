!
!  -- dalton/sirius/sircav.F --
!     (Luca Frediani)
! Extracted from DALTON
! RDR 300114 Wrap up in a F90 module. 
!
    module pedra_cavity

    implicit none

    public polyhedra_driver

    private

    ! Some print levels
    integer :: iprsol = 0

    contains

    subroutine polyhedra_driver(work,lwork)

    use pedra_utils, only : errwrk

#include <pcm_maxorb.h>
#include <pcm_maxaqn.h>
#include <pcm_priunit.h>
#include <pcm_iratdef.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_symmet.h>

    real(8) :: work(*)
    integer :: lwork

    logical :: some
    integer :: numts, numsph, natm, numver
    integer :: loadfm, lintsp, lvert, lcentr, lwrk
    integer :: lnewsp, licav1, licav2, lx, ly, lz
    integer :: ljtr, lcv, lnperm, last
    character(len=8) :: memory_used_format = '(A, I10)'
    integer, allocatable :: intsph(:, :), newsph(:, :)
    integer, allocatable :: icav1(:), icav2(:)
    integer, allocatable :: jtr(:, :), nperm(:, :)

#include <pcm_pcm.h>
#include <pcm_pcmlog.h>

!     ----- set memory pointers for polyhedra setup -----

    SOME = IPRPCM.NE.-5

!     maximum number of tessalations is not known, take worst case
!     maximum number of spheres is not known, take worst case

    NUMTS  = MXTS
    NUMSPH = MXSP
    NATM   = MXCENT
    NUMVER = MXVER

    LOADFM = 1
    LINTSP = LOADFM + 1
    LVERT  = LINTSP + (NUMTS*10 + 1)/IRAT
    LCENTR = LVERT  + NUMTS*10*3
    LNEWSP = LCENTR + NUMTS*10*3
    LICAV1 = LNEWSP + (NUMSPH*2 + 1)/IRAT
    LICAV2 = LICAV1 + MXCENT
    LX     = LICAV2 + MXCENT
    LY     = LX     + NUMTS
    LZ     = LY     + NUMTS
    LJTR   = LZ     + NUMTS
    LCV    = LJTR   + NUMTS*3
    LNPERM = LCV    + NUMVER*3
    LAST   = LNPERM + NUMSPH*8
    IF (LAST > LWORK) CALL ERRWRK('polyhedra_driver',LAST,LWORK, lvpri)
    LWRK   = LWORK - LAST + 1
    IF(SOME) then 
            wRITE(LVPRI, memory_used_format) ' MEMORY USED TO GENERATE CAVITY =', LAST
    end if

    allocate(intsph(numts, 10))
    allocate(newsph(numsph, 2))
    allocate(icav1(natm))
    allocate(icav2(natm))
    allocate(jtr(numts, 3))
    allocate(nperm(numsph, maxrep))

    call polyhedra(intsph, WORK(LVERT), WORK(LCENTR), newsph, &
    icav1, icav2, WORK(LX), WORK(LY), WORK(LZ), &
    jtr, WORK(LCV), NUMTS, NUMSPH, NUMVER, NATM, SOME, &
    WORK(LAST), LWRK, nperm)

    end subroutine polyhedra_driver

    subroutine polyhedra(INTSPH,VERT,CENTR,NEWSPH,ICAV1,ICAV2, XVAL,YVAL, &
    ZVAL,JTR,CV,NUMTS,NUMSPH,NUMVER,NATM,SOME,WORK,LWORK,NPERM)

    use pedra_dblas, only: dzero
    use pedra_print, only: output
    use pedra_cavity_derivatives, only: cavder

#include <pcm_maxorb.h>
#include <pcm_maxaqn.h>
#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_symmet.h>

    logical :: some

#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_nuclei.h>
#include <pcm_pcmnuclei.h>

    integer :: numts, natm, numsph, numver, lwork
    integer :: intsph(numts, 10), newsph(numsph, 2), icav1(natm), icav2(natm)
    real(8) :: vert(numts, 10, 3), centr(numts, 10, 3), cv(numver, *)
    real(8) :: xval(numts), yval(numts), zval(numts)
    real(8) :: pp(3), pp1(3), pts(3, 10), ccc(3, 10)
    real(8) :: work(lwork)
    integer :: jtr(numts, 3), nperm(numsph, *) 

    real(8), parameter :: d0 = 0.0d0
    real(8), parameter :: first = 0.0174533d0
    real(8) :: area, cosom2, fc, fc1, hh, omg, prod, r2gn
    real(8) :: reg, reg2, regd2, ren, rend2, reo, reo2
    real(8) :: rep, rep2, repd2, rgn, rij, rij2, rik, rik2
    real(8) :: rjd, rjk, rjk2, rtdd, rtdd2, senom, sp 
    real(8) :: test, test1, test2, test3, test7, test8
    real(8) :: xen, xi, xj, xn, yen, yi, yj, yn, zen, zi, zj, zn
    integer :: i, icoord, idisp, ii, ipflag, iprcav, iptype
    integer :: its, itseff, itsnum, itypc, iv, iver, j, jj, jcor
    integer :: k, katom, kg, idisrep, kidx, klast, kord
    integer :: kp, lend, lgeom, lmass, lwrk, n, n1, n2, n3
    integer :: natsph, ncav1, ncav2, ne, nes, net, nev, nn
    integer :: nsfe, nsfer, nv
    integer, allocatable :: idxpri(:)
    
    real(8), dimension(3, 3) ::  rotcav = reshape([1.0d0, 0.0d0, 0.0d0,  &
                                                   0.0d0, 1.0d0, 0.0d0,  &
                                                   0.0d0, 0.0d0, 1.0d0], [3, 3])
          

!     Se stiamo costruendo una nuova cavita' per il calcolo del
!     contributo dispersivo e repulsivo:

    IDISREP = 0
    IPRCAV = 0


! se icesph=0
!     legge i raggi dall'input e fa coincidere i centri con gli atomi
! se icesph=1
!     legge centri e raggi dall'input
! se icesph=2
!     legge i raggi dall'input e fa coincidere i centri con
!     alcuni atomi definiti dall'indice ina(k) con k=1,NESFP
!     es: xe(k)=c(1,ina(k))

! f      IF(ICESPH.LE.0) THEN
! f         DO J=1,NESFP
! f
! fc           XE(J)=CORD(1,J)
! fc           YE(J)=CORD(2,J)
! fc           ZE(J)=CORD(3,J)
! fClf            write(lvpri,*) 'icesph',icesph,nesfp
! f            XE(J)=PCMCORD(1,J)
! f            YE(J)=PCMCORD(2,J)
! f            ZE(J)=PCMCORD(3,J)
! f            INA(J)=J
! f         ENDDO
! fC
! f      END IF
! f      IF(ICESPH.EQ.2)THEN
! f         DO J=1,NESFP
! f           XE(J)=PCMCORD(1,INA(J))
! f           YE(J)=PCMCORD(2,INA(J))
! f           ZE(J)=PCMCORD(3,INA(J))
! f         ENDDO
! f      END IF
! f      NESF=NESFP
!      print * ,'at the beginning of pedra', nesf
!      write(lvpri,*) 'at the beginning of pedra', nesf
!      do i=1,nesf
!         write(lvpri,1111) xe(i), ye(i), ze(i), re(i)
!      enddo

!  PEDRA prevede che i dati geometrici siano espressi in ANGSTROM :
!  vengono trasformati, e solo alla fine i risultati tornano in bohr.

    90 CONTINUE
    DO I=1,NESFP
        XE(I)=XE(I)
        YE(I)=YE(I)
        ZE(I)=ZE(I)
        RE(I) = RIN(I) * ALPHA(I)
        write (lvpri,*) xe(i), ye(i), ze(i), re(i)
    ENDDO

    CALL DZERO(VERT,NUMTS*10*3)
    CALL DZERO(CENTR,NUMTS*10*3)

!                   creation of new spheres

    DO N = 1, NESF
        NEWSPH(N,1) = 0
        NEWSPH(N,2) = 0
    ENDDO

    ITYPC = 0
    OMG=OMEGA*FIRST
    SENOM=Sin(OMG)
    COSOM2=(cos(OMG))**2
    RTDD=RET+RSOLV
    RTDD2=RTDD*RTDD
    NET=NESF
    NN=2
    NE=NESF
    NEV=NESF
    GO TO 100
    110 NN=NE+1
    NE=NET
    100 CONTINUE

!     check on the number of spheres

    WRITE(LVPRI,*) 'Number of extra spheres = ', NE
    IF (NE > MXSP) THEN
        WRITE(LVPRI,*) NE
        WRITE(LVPRI,*) 'Too many spheres (Greater than MXSP=200).'
        STOP
    ENDIF


    DO 120 I=NN,NE
        NES=I-1
        DO 130 J=1,NES
            RIJ2=(XE(I)-XE(J))**2+ &
            (YE(I)-YE(J))**2+ &
            (ZE(I)-ZE(J))**2
            RIJ=SQRT(RIJ2)
            RJD=RE(J)+RSOLV
            TEST1=RE(I)+RJD+RSOLV
            IF(RIJ >= TEST1) GO TO 130
            REG=DMAX1(RE(I),RE(J))
            REP=DMIN1(RE(I),RE(J))
            REG2=REG*REG
            REP2=REP*REP
            TEST2=REP*SENOM+sqrt(REG2-REP2*COSOM2)
            IF(RIJ <= TEST2) GO TO 130
            REGD2=(REG+RSOLV)*(REG+RSOLV)
            TEST3=(REGD2+REG2-RTDD2)/REG
            IF(RIJ >= TEST3) GO TO 130
            DO 140 K=1,NEV
                IF(K == J .OR. K == I) GO TO 140
                RJK2=(XE(J)-XE(K))**2+ &
                (YE(J)-YE(K))**2+ &
                (ZE(J)-ZE(K))**2
                IF(RJK2 >= RIJ2) GO TO 140
                RIK2=(XE(I)-XE(K))**2+ &
                (YE(I)-YE(K))**2+ &
                (ZE(I)-ZE(K))**2
                IF(RIK2 >= RIJ2) GO TO 140
                RJK=sqrt(RJK2)
                RIK=sqrt(RIK2)
                SP=(RIJ+RJK+RIK)/2.0D0
                HH=4*(SP*(SP-RIJ)*(SP-RIK)*(SP-RJK))/RIJ2
                REO=RE(K)*FRO
                IF(K >= NE)REO=0.0002D0
                REO2=REO*REO
                IF(HH < REO2) GO TO 130
            140 END DO
            REPD2=(REP+RSOLV)**2
            TEST8=SQRT(REPD2-RTDD2)+sqrt(REGD2-RTDD2)
            IF(RIJ <= TEST8)GO TO 150
            REND2=REGD2+REG2-(REG/RIJ)*(REGD2+RIJ2-REPD2)
            IF(REND2 <= RTDD2) GO TO 130
            REN=sqrt(REND2)-RSOLV
            FC=REG/(RIJ-REG)
            TEST7=REG-RE(I)
            KG=I
            KP=J
            IF(TEST7 <= 0.000000001D0) GO TO 160
            KG=J
            KP=I
            160 FC1=FC+1.0
            XEN=(XE(KG)+FC*XE(KP))/FC1
            YEN=(YE(KG)+FC*YE(KP))/FC1
            ZEN=(ZE(KG)+FC*ZE(KP))/FC1
            ITYPC = 1
            GO TO 170
            150 R2GN=RIJ-REP+REG
            RGN=R2GN/2.0D0
            FC=R2GN/(RIJ+REP-REG)
            FC1=FC+1.0D0
            TEST7=REG-RE(I)
            KG=I
            KP=J
            IF(TEST7 <= 0.000000001D0) GO TO 180
            KG=J
            KP=I
            180 XEN=(XE(KG)+FC*XE(KP))/FC1
            YEN=(YE(KG)+FC*YE(KP))/FC1
            ZEN=(ZE(KG)+FC*ZE(KP))/FC1
            REN=sqrt(REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ))-RSOLV
            170 NET=NET+1
            XE(NET)=XEN
            YE(NET)=YEN
            ZE(NET)=ZEN
            RE(NET)=REN
        
        !     Nella matrice NEWSPH(NESF,2) sono memorizzati i numeri delle
        !     sfere "generatrici" della nuova sfera NET: se la nuova sfera e'
        !     del tipo A o B entrambi i numeri sono positivi, se e' di tipo
        !     C il numero della sfera "principale" e' negativo
        !     (per la definizione del tipo si veda JCC 11, 1047 (1990))
        
            IF(ITYPC == 0) THEN
                NEWSPH(NET,1) = KG
                NEWSPH(NET,2) = KP
            ELSEIF(ITYPC == 1) THEN
                NEWSPH(NET,1) = - KG
                NEWSPH(NET,2) = KP
            ENDIF
        
        130 END DO
        NEV=NET
    120 END DO
    IF(NET /= NE) GO TO 110
    NESF=NET

!     costruzione della tabella di permutazione delle sfere

! f      WRITE(LVPRI,*)'MAXREP', MAXREP

    CALL SPHPER(NESF,NESFP,NUMSPH,NPERM,NEWSPH,XE,YE,ZE,RE)

! f Determination of eigenvalues and eigenvectors of the inertia tensor

!      KATOM = NATOMS + NFLOAT

!     Modification by Ville & Chris for new cavity gen
    KATOM = NESFP
    LGEOM = 1
    LMASS = LGEOM + 3 * KATOM
    LEND  = LMASS + KATOM
    LWRK =  LWORK - LEND
! f      write(lvpri,*) 'katom',katom,lgeom,lmass,lend
    IF (LWRK < 0) THEN
        WRITE(*,*) '2'
        STOP
    ENDIF

    IF (KATOM > 1) THEN
        CALL PCMTNS(ROTCAV,WORK(LGEOM),WORK(LMASS),KATOM)
    END IF

!    Division of the surface into tesserae

    VOL=D0
    STOT=D0

!     Controlla se ciascuna tessera e' scoperta o va tagliata

    NN = 0
    DO 300 NSFE = 1, NESF
        XEN = XE(NSFE)
        YEN = YE(NSFE)
        ZEN = ZE(NSFE)
        REN = RE(NSFE)
        IF(IPRCAV >= 10) THEN
            WRITE(LVPRI,1242) NSFE
            WRITE(LVPRI,1243) XEN,YEN,ZEN,REN
        END IF
    !     Default options (same as traditional GEPOL)
        IPtype = 2
        IPFlag = 0
        ITSNUM = 60
    ! f Which type of tessellation?
    !         IF(IPOLYG.GT.0) THEN
    !            IPFlag = 0
    !            ITSNUM = IPolyg
    !         ELSEIF(IPolyg.lt.0) THEN
    !            IPFlag = 1
    !            TsAre = -1.0D-03*IPolyg
    ! f         ENDIF
        CALL POLYGEN(IPFLAG,AREATS,ITSNUM,XEN,YEN,ZEN,REN,ITSEFF, &
        CV,JTR,NPerm,NSFE,NUMTS,NUMSPH,NUMVER,ROTCAV)
        IF(IPRCAV >= 10) WRITE(LVPRI,*)'AFTER POLYGEN_. ITSEFF=',ITSEFF
        DO 310 ITS = 1, ITSEFF
            N1 = JTR(ITS,1)
            N2 = JTR(ITS,2)
            N3 = JTR(ITS,3)
            PTS(1,1)=CV(N1,1)
            PTS(2,1)=CV(N1,2)
            PTS(3,1)=CV(N1,3)
            PTS(1,2)=CV(N2,1)
            PTS(2,2)=CV(N2,2)
            PTS(3,2)=CV(N2,3)
            PTS(1,3)=CV(N3,1)
            PTS(2,3)=CV(N3,2)
            PTS(3,3)=CV(N3,3)
            NV=3
            IF(IPRCAV >= 10) THEN
                WRITE(LVPRI,1244) N1,N2,N3
                WRITE(LVPRI,*) 'VERTICES'' COORDINATES'
                CALL OUTPUT(PTS,1,3,1,3,3,3,1,LVPRI)
            END IF
            DO JJ = 1, 3
                PP(JJ) = D0
                PP1(JJ) = D0
            ENDDO
        
        !     Per ciascuna tessera, trova la porzione scoperta e ne
        !     calcola l'area con il teorema di Gauss-Bonnet; il punto rappresentativo
        !     e' definito come media dei vertici della porzione scoperta di tessera
        !     e passato in PP (mentre in PP1 ci sono le coordinate del punto sulla
        !     normale interna).
        !     I vertici di ciascuna tessera sono conservati in VERT(3,10,MxTs), il
        !     numero di vertici di ciascuna tessera e' in NVERT(MxTs), e i centri
        !     dei cerchi di ciascun lato sono in CENTR(3,10,MxTs).
        !     In INTSPH(MxTs,10) sono registrate le sfere a cui appartengono i lati
        !     delle tessere.
        
            IF (IPRCAV > 15) THEN
                write(lvpri,9100) NN + 1
                DO IVER = 1, 3
                    WRITE(LVPRI,9110) IVER, (PTS(JCOR,IVER),JCOR=1,3)
                END DO
            END IF
            CALL TESSERA(NSFE,NV,PTS,CCC,PP,PP1,AREA,INTSPH,NUMTS)
            IF(IPRCAV >= 10) THEN
                WRITE(LVPRI,*) 'AFTER TESSERA_ ROUTINE'
                WRITE(LVPRI,1245) NV
                WRITE(LVPRI,*) 'VERTICES'' COORDINATES'
                CALL OUTPUT(PTS,1,3,1,NV,3,NV,1,LVPRI)
            END IF
            IF(AREA == D0) THEN
                WRITE(LVPRI,*) 'ZERO AREA IN TESSERA_', NN + 1
                GOTO 310
            END IF
            IF (NN >= MXTS) THEN
                WRITE(LVPRI,*) 'TOO MANY TESSERA_E IN PEDRA'
                WRITE(LVPRI,*) 'NN=',NN,'  MXTS=',MXTS
                STOP
            END IF
            NN = NN + 1
            XTSCOR(NN) = PP(1)
            YTSCOR(NN) = PP(2)
            ZTSCOR(NN) = PP(3)
            XVAL(NN) = PP1(1)
            YVAL(NN) = PP1(2)
            ZVAL(NN) = PP1(3)
            AS(NN) = AREA
            ISPHE(NN) = NSFE
            NVERT(NN) = NV
            DO IV = 1, NV
                DO JJ = 1, 3
                    VERT(NN,IV,JJ) = PTS(JJ,IV)
                    CENTR(NN,IV,JJ) = CCC(JJ,IV)
                ENDDO
            ENDDO
            DO IV = 1, NV
                INTSPH(NN,IV) = INTSPH(NUMTS,IV)
            ENDDO
            IF(IPRCAV >= 10) THEN
                WRITE(LVPRI,1246) NN
                WRITE(LVPRI,1247) XTSCOR(NN),YTSCOR(NN),ZTSCOR(NN)
                WRITE(LVPRI,1248) XVAL(NN),YVAL(NN),ZVAL(NN)
                WRITE(LVPRI,1249) AS(NN),ISPHE(NN),NVERT(NN)
            END IF
        310 END DO
    ! f
        9100 format('Before Tessera n.',I4)
        9110 format('XYZ vert n.',I4,3f15.9)
        9120 format('After Tessera n.',I4)
        9130 format('XYZ centr n.',I4,3f15.9)
        9140 format('XYZA',4f15.9)
    ! f
    300 END DO
    NTS = NN

!     Verifica se due tessere sono troppo vicine
    TEST = 0.02D0
    TEST2 = TEST*TEST
    DO 400 I = 1, NTS-1
        IF(AS(I) == D0) GOTO 400
        XI = XTSCOR(I)
        YI = YTSCOR(I)
        ZI = ZTSCOR(I)
        II = I + 1
        DO 410 J = II , NTS
            IF(ISPHE(I) == ISPHE(J)) GOTO 410
            IF(AS(J) == D0) GOTO 410
            XJ = XTSCOR(J)
            YJ = YTSCOR(J)
            ZJ = ZTSCOR(J)
            RIJ = (XI-XJ)**2 + (YI-YJ)**2 + (ZI-ZJ)**2
            IF(RIJ > TEST2) GOTO 410
            WRITE(LVPRI,9010) I,J,RIJ,TEST2
        ! a vedere  IF(AS(I).LT.AS(J)) AS(I) = D0
        ! a vedere  IF(AS(I).GE.AS(J)) AS(J) = D0
        410 END DO
    400 END DO

! f Replication of geometrical parameters according to simmetry

    CALL REPCAV(VERT,CENTR,NPERM,NUMTS,NUMSPH)

! Prepare data for geomview

! f      write(lvpri,*) 'vertici delle tessere'
!      do i=1,nts
!         do j=1,nvert(i)
!            write(lvpri,1241) i,j,(vert(i,j,k),k=1,3)
!         enddo
!      enddo
!      write(lvpri,*) 'centers of tesserae'
!      do i=1,nts
!         write(lvpri,1240) i, xtscor(i),ytscor(i),ztscor(i),as(i)
!      enddo
!      call flshfo(lvpri)

! f      IF (IPRPCM .GT. 10) THEN
    IF ( .TRUE. ) THEN
        KORD  = 1
        KIDX  = KORD + 4 * NTS
        KLAST = KIDX + NTS
        allocate(idxpri(nts))
        IF (KLAST < LWORK) THEN
            CALL DZERO(WORK,KLAST)
            call ORDPCM(lvpri,nts,xtscor,ytscor,ztscor,as,work(kord),idxpri,work(klast),lwork)
        ELSE
            WRITE(LVPRI,*) &
            'WARNING: NOT ENOUGH MEM IN PEDRA TO PRINT THE CAVITY'
        END IF
    END IF
    1234 format(2i4,3f12.6)
    1235 format('XTSCOR(',i4,') = ',f18.12)
    1236 format('YTSCOR(',i4,') = ',f18.12)
    1237 format('ZTSCOR(',i4,') = ',f18.12)
    1238 format('AS(',i4,') = ',f18.12)
    1239 format('ISPHE(',i4,') = ',I3)
    1240 format(i4,4f12.6)
    1241 format(2i4,3f15.9)
    1242 FORMAT('Now tessellating sphere n.',i4)
    1243 FORMAT('Coordinates: X=',F12.8,'  Y=',F12.8,'  Z=',F12.8, &
    '  R=',F12.8)
    1244 FORMAT('VERTICES INDICES: ',3I6)
    1245 FORMAT('NO. OF VERTICES: ',I6)
    1246 FORMAT('ADDING TESSERA_ DATA TO THE LIST. NN=',I6)
    1247 FORMAT('CENTER: ',3F15.11)
    1248 FORMAT('BOH???: ',3F15.11)
    1249 FORMAT('AREA=',F10.8,'ISPHE=',I4,'NVERT=',I3)
    CALL PLOTCAV(VERT,NUMTS)

!***********************************************************
!     Calcola il volume della cavita' con la formula (t. di Gauss):
!                V=SOMMAsulleTESSERE{A r*n}/3
!     dove r e' la distanza del punto rappresentativo dall'origine,
!     n e' il versore normale alla tessera, A l'area della tessera,
!     e * indica il prodotto scalare.
!***********************************************************
    VOL = D0
    DO ITS = 1, NTS
        NSFE = ISPHE(ITS)
    !     Trova il versore normale
        XN = (XTSCOR(ITS) - XE(NSFE)) / RE(NSFE)
        YN = (YTSCOR(ITS) - YE(NSFE)) / RE(NSFE)
        ZN = (ZTSCOR(ITS) - ZE(NSFE)) / RE(NSFE)
    !     Trova il prodotto scalare
        PROD = XTSCOR(ITS)*XN + YTSCOR(ITS)*YN + ZTSCOR(ITS)*ZN
        VOL = VOL + AS(ITS) * PROD / 3.D0
    ENDDO
!***********************************************************
!     Stampa la geometria della cavita'
!***********************************************************
    STOT=D0
    CALL DZERO(SSFE,NESF)
    DO I=1,NTS
        K=ISPHE(I)
        SSFE(K)=SSFE(K)+AS(I)
    ENDDO

    IF(SOME) THEN
        WRITE(LVPRI,9020) NESF
        DO 520 I=1,NESF
            WRITE(LVPRI,9030) I,XE(I),YE(I),ZE(I),RE(I),SSFE(I)
            STOT=STOT+SSFE(I)
        520 END DO
        WRITE(LVPRI,9040) NTS,STOT,VOL
    END IF

!     ----- set up for possible gradient calculation -----
!           DONE ONLY IF WE HAVE SPHERES ON ATOMS

    IF(ICESPH /= 1) THEN
    
        IF(NESFP > NUCDEP) THEN
            WRITE(LVPRI,9050) NESFP,NUCDEP
            WRITE(*,*) '3'
            STOP
        END IF
    
        CALL DZERO(DERCEN,MXSP*MXCENT*3*3)
        CALL DZERO(DERRAD,MXSP*MXCENT*3)
        DO NSFE = 1,NUCDEP
            NATSPH=0
            NSFER=NSFE
            IF(ICESPH == 2)THEN
                NATSPH=1
                DO JJ=1,NESFP
                    IF(INA(JJ) == NSFE)THEN
                        NATSPH=0
                        NSFER=JJ
                    ENDIF
                ENDDO
            ENDIF
            IF(NATSPH == 0) then
                DO ICOORD = 1, 3
                    CALL CAVDER(NSFE,NSFER,ICOORD,INTSPH,NEWSPH)
                ENDDO
            ENDIF
        ENDDO
    END IF

    IF(IPRSOL > 5) THEN
        WRITE(LVPRI,9070)
        WRITE(LVPRI,9090) (I,ISPHE(I),AS(I),XTSCOR(I),YTSCOR(I), &
        ZTSCOR(I),XTSCOR(NTS+I),YTSCOR(NTS+I), &
        ZTSCOR(NTS+I), I=1,NTS)
    END IF

!     Call the routine which checks if the cavity is single or divided.

    call cavspl(icav1,icav2,ncav1,ncav2,npcmn,some)

!     The dispersion calculation is allowed only in the case of
!     single cavity.

    IF(NCAV2 /= 0) IDISP=0

    IDISREP = 0
    IRETCAV = 0
    IF (IDISP == 2) IDISP = 1

    9000 FORMAT(10X,'-- CENTER OF CHARGE --'/ &
    ' X =',F8.4,' AU  Y =',F8.4,' AU  Z =',F8.4,' AU')
    9010 FORMAT(/' WARNING: The distance between center of tessera',I4, &
    'and',I4,' is ',F8.6,', less than ',F8.6,' AU'/)
    9020 FORMAT(/' Total number of spheres =',I5/ &
    ' Sphere             Center  (X,Y,Z) (AU)            ', &
    '   Radius (AU)      Area (AU^2)')
    9030 FORMAT(I4,4F15.9,F15.9)
    9040 FORMAT(/' Total number of tesserae =',I8 &
    /' Surface area =',F14.8,' (AU^2)    Cavity volume =', &
    F14.8,' (AU^3)')
    9050 FORMAT( ' PEDRA: CONFUSION ABOUT SPHERE COUNTS. NESFP, NATM=',2I6)
    9060 FORMAT(/' ADDITIONAL MEMORY NEEDED TO SETUP GRADIENT RUN=',I10)
    9061 FORMAT(/' ADDITIONAL MEMORY NEEDED TO SETUP IEF RUN=',I10)
    9070 FORMAT(/' ***  PARTITION OF THE SURFACE  ***' &
    //' TESSERA_  SPHERE   AREA   X Y Z TESSERA CENTER  ', &
    'X Y Z NORMAL VECTOR')
    9090 FORMAT(2I4,7F12.7)

    end subroutine polyhedra
    
    subroutine sphper(nesf,nesf0,numsph,nperm,newsph,xe,ye,ze,re)
  
    use pedra_ibtfun

#include <pcm_maxorb.h>
#include <pcm_maxaqn.h>
#include <pcm_mxcent.h>
#include <pcm_priunit.h>
#include <pcm_symmet.h>

! f#include <pcm_trans.h>
    
    integer :: nesf, nesf0, numsph
    real(8) :: xe(*), ye(*), ze(*), re(*), v1(3), v2(3)
    integer :: nperm(numsph,*), newsph(numsph,2)

    real(8) :: diff1, r1
    integer :: i, j, k, l, n1, n2, nesf1

!  crea la tabella di permutazione delle sfere
!  e completa per simmetria le sfere aggiunte

          
    NESF1=NESF

    DO I=1,NESF
        V1(1)=XE(I)
        V1(2)=YE(I)
        V1(3)=ZE(I)
        R1=RE(I)
    
    !  ciclo sulle operazioni di simmetria
    
                
        DO J=1,MAXREP
            DO K=1,NESF
                DO L=1,3
                    V2(L)=PT(IBTAND(ISYMAX(L,1),J))*V1(L)
                ENDDO
            ! f
            !            DIFF1=
            !     &        SQRT((XE(K)-V2(1))**2+(YE(K)-V2(2))**2+(ZE(K)-V2(3))**2)
            !            DIFF2=ABS(R1-RE(K))
            !            IF ((DIFF1.LT.1.0d-3).AND.(DIFF2.LT.1.0d-3)) THEN
            ! f
                DIFF1 = SQRT((XE(K)-V2(1))**2 + (YE(K)-V2(2))**2 + &
                (ZE(K)-V2(3))**2 + (RE(K)-R1)**2)
                IF ((DIFF1 < 1.0d-3)) THEN
                    NPERM(I,J+1)=K
                    GOTO 10
                ENDIF
            ENDDO
        
        !  caso in cui occorre completare la creazione di sfere aggiunte
        !  oppure le sfere iniziali sono messe male
        
            If (I < NESF0) THEN
                WRITE(LVPRI,*) 'CAVITIY IS NOT CONSISTENT WITH SYMMETRY'
                STOP
            ELSE
                NESF1=NESF1+1
                NPerm(i,j+1)=NESF1
                XE(NESF1)=v2(1)
                YE(NESF1)=v2(2)
                ZE(NESF1)=v2(3)
                RE(NESF1)=r1
                N1=NPERM(IAbs(NewSph(i,1)),j+1)
                N2=NPERM(IAbs(NEwSph(i,2)),j+1)
                IF (NewSph(i,1) < 0) THEN
                    N1=-N1
                    N2=-N2
                ENDIF
                NewSph(NESF1,1)=N1
                NewSph(NESF1,2)=N2
                WRITE(LVPRI,*) 'new!',NESF1,N1,N2
            ENDIF
            10 CONTINUE
        ENDDO
    ENDDO

!  eventuale completamento della tabella delle permutazioni per le
!  nuove sfere

    DO I=NESF+1,NESF1
        V1(1)=XE(I)
        V1(2)=YE(I)
        V1(3)=ZE(I)
        R1=RE(I)
                
    !  ciclo sulle operazioni di simmetria
                
        DO J=1,MAXREP
            DO K=1,NESF1
                DO L=1,3
                    V2(L)=PT(IBTAND(ISYMAX(L,1),J))*V1(L)
                ENDDO
                DIFF1 = SQRT((XE(K)-V2(1))**2 + (YE(K)-V2(2))**2 + &
                (ZE(K)-V2(3))**2 + (RE(K)-R1)**2)
                IF ((DIFF1 < 1.0d-3)) THEN
                    NPERM(I,J+1)=K
                    GOTO 20
                ENDIF
            ENDDO
            WRITE(*,*) '4'
            STOP

            20 CONTINUE
        ENDDO
    ENDDO
          
    NESF=NESF1

!  vengono segnate le sfere da tassellare: la prima in ordine di ogni
!  classe di sfere equivalenti per simmetria

    DO I=1,NESF
        NPERM(I,1)=1
        DO K=1,MAXREP
            IF (NPERM(I,K+1) < I) THEN
                NPERM(I,1)=0
                GOTO 30
            ENDIF
        ENDDO
        30 CONTINUE
    ENDDO
    
    END SUBROUTINE SPHPER
    
    SUBROUTINE POLYGEN(IPFLAG,TSARE,ITSNUM,XEN,YEN,ZEN,REN,ITSEFF,CV,JTR,NPERM,NSFE,NUMTS,NUMSPH,NUMVER,ROTCAV)

    use pedra_ibtfun

#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_symmet.h>

! Polygen: a program to generate spherical polyhedra with triangular
! faces. An equilater division algorithm is used.

    integer :: ipflag, itsnum, itseff, nsfe, numts, numsph, numver
    integer :: nperm(numsph, *), jtr(numts, *)
    real(8) :: cv(numver, *)
    real(8) :: tsare, xen, yen, zen, ren
    
    real(8) :: v1(3), v2(3), v3(3), rotcav(3, 3)
    integer :: itrvo(60,3) = 0
    integer :: itreo(60,3) = 0 
    integer :: iedo(90,2)  = 0
    integer :: oldtr(100,100), ednew(90,100), trnew(60,100,100)

    real(8), parameter :: d0 = 0.0d0
    real(8), parameter :: d1 = 1.0d0
    real(8), parameter :: pi = acos(-1.0d0)
    real(8) :: alpha, beta, cos1, cos2, costheta, dl, dm, dn, dnf, dnorm
    real(8) :: sintheta, theta
    integer :: i, ii, isymop, j, jj, jsymop, k, l, m, n, ne0, nf
    integer :: noppt, nt, nt0, ntpt, ntra, nv, nvpt


!  Se la sfera non deve essere tassellata ritorna a pedra con 0 tessere

! f      IF (NPerm(NSFE,1).EQ.0.or.nsfe.ne.2) THEN
    IF (NPerm(NSFE,1) == 0) THEN
        ITsEff=0
    ! f         print *, 'return!', ITsEff,NPerm(NSFE,1)
        RETURN
    ENDIF

!  A seconda delle due opzioni di funzionamento (Area ottimale o
!  numero di tessere ottimale ) vengono stabilite il tipo di poliedro
!  e la  frequenza di divisione (NF)

!  IPflag = 0 Numero di tessere: si genera il poliedro disponibile
!           con il numero di tessere piu' prossimo a ITSNUM
!  IPflag = 1  Area: si genera il poliedro con il numero di tessere tali
!           da avere un area il piu' simile possibile

    IPFLAG = 1
    IF (IPflag == 1) THEN
        ITSNUM = INT (4.0d0 * pi * REN**2 / TsAre + 0.5D0 )
    ENDIF

    IF (ITSNUM > MxTs) THEN
        WRITE(*,*) '5', ITSNUM, MxTs
        STOP
    ENDIF

!  costruisce la tassellazione iniziale che dipende dal gruppo di
!  simmetria e (per C1) dal numero di tessere richieste
    nt0 = 1
    nv = 3
    ne0 = 3
    itrvo(1,1)=1
    itrvo(1,2)=2
    itrvo(1,3)=3
    iedo(1,1)=1
    iedo(1,2)=2
    iedo(2,1)=2
    iedo(2,2)=3
    iedo(3,1)=1
    iedo(3,2)=3
    itreo(1,1)=1
    itreo(1,2)=2
    itreo(1,3)=3

    CV(1,1)=D1
    CV(1,2)=D0
    CV(1,3)=D0

    CV(2,1)=D0
    CV(2,2)=D1
    CV(2,3)=D0

    CV(3,1)=D0
    CV(3,2)=D0
    CV(3,3)=D1

!  determination of NF

          
    NF = INT(SQRT(D1 * ITSNUM / ( D1 * NT0  * 8 )) + 0.5d0)
    IF (NF <= 1) NF=2
    IF (NF <= 2) WRITE(LVPRI,1000)
    IF (NF >= 8) WRITE(LVPRI,1100)
    DNF = DFLOAT(NF)
! eck WRITE(LVPRI,*) 'dopo polydata',NT0,NE0,NV,NF

    ITsEff=NT0*NF**2

!  -nuovi vertici posti lungo i vecchi spigoli
!  -le regole di calcolo derivano da  calcoli di algebra
!   lineare e trigonometria
!  -i vertici vengono memorizzati in ordine progressivo ed riferiti
!  tramite ednew al vertice di appartenenza

    NVPT=NV+1
    DO J = 1,NE0
        DO K = 1,3
            v1(k)=CV(IEDO(j,1),k)
            v2(k)=CV(IEDO(j,2),k)
        ENDDO
        costheta= &
        ( v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)) / &
        (sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2) * &
        sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2))
        theta=acos(costheta)
        sintheta=sin(theta)
        DO l=1,NF-1
            DL = DFLOAT(L)
        ! F          m=NF-l
            cos1=cos(theta*Dl/DNF)
            cos2=cos(theta*(DNF-Dl)/DNF)
            alpha=(cos1-costheta*cos2)/sintheta**2
            beta=(cos2-costheta*cos1)/sintheta**2
            DM = 0.0D0
            ALPHA = DFLOAT(NF-L)
            BETA  = DFLOAT(L)
            DO k=1,3
                v3(k)=alpha*v1(k)+beta*v2(k)
            ENDDO
            dnorm = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
            DO k=1,3
                v3(k) = v3(k) /dnorm
            ENDDO
            DO K = 1,3
                CV(NVPT,K)=V3(K)
            ENDDO
            ednew(j,l+1)=NVPT
            NVPT=NVPT+1
        ENDDO
    ENDDO
          

!  -nuovi vertici non posti lungo i vecchi spigoli
!  -a partire dai vertici in ednew secondo regole analoghe alle
!  precedenti, vengono memorizzati a seconda del triangolo in trnew
!  trnew(triangolo,fila,n ordine)=nvertice

! f      DO J=1,NT0
    j=1
    ii=1
    jj=3
    DO L =3,NF
        DL = DFLOAT(L)
        DO N=1,l-2
            DM = DFLOAT(NF - L + 1)
            DL = DFLOAT(L - 1 - N)
            DN = DFLOAT(N)
        !            CV(NVPT,1) = SQRT(DM/DNF)
        !            CV(NVPT,2) = SQRT(DL/DNF)
        !            CV(NVPT,3) = SQRT(DN/DNF)
            CV(NVPT,1) = DM
            CV(NVPT,2) = DL
            CV(NVPT,3) = DN
            DNORM = SQRT(CV(NVPT,1)**2 + CV(NVPT,2)**2 + CV(NVPT,3)**2)
            CV(NVPT,1) = CV(NVPT,1) / DNORM
            CV(NVPT,2) = CV(NVPT,2) / DNORM
            CV(NVPT,3) = CV(NVPT,3) / DNORM
            trnew(j,l,n+1)=NVPT
            NVPT=NVPT+1
        ENDDO
    ENDDO
! f       ENDDO
    NV=NVPT-1

!  -ora per ogni triangolo originario vengono posti nella matrice oldtr
!  i numeri di ordine dei vertici originali,creati lungo gli spigoli e
!  no dei vecchi triangoli secondo lo schema:

!          11                    La disposizione e' secondo
!          | \                   la posizione geometrica del
!          |  \                  vertice.
!          21--22
!          | \  |\               i vertici (1,1)-(NF+1,1)-(NF+1,NF+1) sono
!          |  \ | \              quelli originari.
!          31--32--33
!          | \  |\ | \           i nuovi triangoli sono:
!          |  \ | \|  \          (i,j)-(i+1,j)-(i+1,j+1)
!          41--42--43--44        i=1,...,NF j=1,...,i

!                                (i,j)--(i,j+1),(i+1,j+1)
!                                i=2,...,NF j=1,...,i-1

    NTPT=1
    do 2310 n=1,NT0
    
    !  -1 vecchi spigoli
    
        oldtr(1,1)=ITRVO(n,1)
        oldtr(NF+1,1)=ITRVO(n,2)
        oldtr(NF+1,NF+1)=ITRVO(n,3)
    
    !  -2 nuovi vertici lungo i vecchi spigoli
    
        DO l=2,NF
            oldtr(l,1)=ednew(ITREO(n,1),l)
            oldtr(NF+1,l)=ednew(ITREO(n,2),l)
            oldtr(l,l)=ednew(ITREO(n,3),l)
        ENDDO
    
    !  -3 nuovi vertici non lungo i vecchi spigoli
    
        DO l=3,NF
            DO m=2,l-1
                oldtr(l,m)=trnew(n,l,m)
            ENDDO
        ENDDO
    
    !  -ora si creano i nuovi triangoli
    
        DO i=1,NF
            DO j=1,i
                JTR(NTPT,1)=oldtr(i,j)
                JTR(NTPT,2)=oldtr(i+1,j)
                JTR(NTPT,3)=oldtr(i+1,j+1)
                NTPT=NTPT+1
            ENDDO
        ENDDO
        DO i=2,NF
            DO j=1,i-1
                JTR(NTPT,1)=oldtr(i,j)
                JTR(NTPT,2)=oldtr(i,j+1)
                JTR(NTPT,3)=oldtr(i+1,j+1)
                NTPT=NTPT+1
            ENDDO
        ENDDO
    2310 END DO
    NV=NVPT-1
    NT=NTPT-1

! f Replication generate the right irreducible part of the cavity
!   for the selected group

    CALL PREREP(NV,NT,ITSEFF,CV,JTR,NUMVER,NUMTS)
    DO i=1,nv
        V1(1) = 0.0D0
        V1(2) = 0.0D0
        V1(3) = 0.0D0
        DO J = 1,3
            DO K = 1,3
                V1(J) = V1(J) + ROTCAV(K,J) * CV(I,K)
            ENDDO
        ENDDO
        CV(I,1) = V1(1)
        CV(I,2) = V1(2)
        CV(I,3) = V1(3)
    ! f         WRITE(LVPRI,*)'CV',I
    ! f         WRITE(LVPRI,'(3F10.4)')  (CV(i,II),II=1,3)
    ENDDO


!  replicazione con operatori locali
!        nei casi in cui l'o.l. trasforma la sfera in un altra sfera

    NOpPt=0
    DO ISYMOP = 1, MAXREP
        NTRA=NPerm(NSFE,ISYMOP + 1)
    
    !     the If statment is entered only if the sphere which is equivalent
    !     to the one we are tesselating now under the operation ISYMOP, is
    !     not itself.
    
        IF (NTRA /= NSFE) THEN
        
        !     If the replication of this part of the cavity has already been
        !     done we start checking a new simmetry operation
        
            DO JSYMOP = 1 , ISYMOP-1
                IF(Nperm(NSFE, JSYMOP+1) == NTRA) Goto 100
            ENDDO
            NOpPt=NOpPt+1
        
        !     riproduzione vertici
        
            Do I=1,NV
                II=I+NOpPt*NV
                Do K=1,3
                    CV(II,K) = PT(IBTAND(ISYMAX(K,1),ISYMOP)) * CV(I,K)
                ENDDO
            ENDDO
        
        !  riproduzione topologia
        
            Do i=1,NT
                ii=i+NOpPt*NT
                jj=NOpPt*NV
                Do k=1,3
                    JTR(ii,k)=JTR(i,k)+jj
                ENDDO
            ENDDO
            100 CONTINUE
        ENDIF
    ENDDO

!  aggiornamento indici

    NV=NV*(NOpPt+1)
    ITsEff=ITsEff*(NOpPt+1)

! scrittura della CV

! heckWRITE(LVPRI,*) 'OFF'
! heckWRITE(LVPRI,*) NV,ITsEff,NV
    DO i=1,NV
        CV(I,1) = CV(I,1) * REN + XEN
        CV(I,2) = CV(I,2) * REN + YEN
        CV(I,3) = CV(I,3) * REN + ZEN
    ! f         WRITE(LVPRI,*)'CV',I
    ! f         WRITE(LVPRI,'(3F10.4)')  (CV(i,II),II=1,3)
    ENDDO

    1000 FORMAT(/' ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN', &
    /' IT IS VALUABLE ALMOST ONLY FOR TESTING')
    1100 FORMAT(/' ** WARNING ** A VERY EXPENSIVE TESSELATION ', &
    'HAS BEEN CHOSEN', &
    /' IT WILL PROBABLY PRODUCE A HUGE AMOUNT OF TESSERA_E')

    END SUBROUTINE POLYGEN
    
    SUBROUTINE REPCAV(VERT,CENTR,NPERM,NUMTS,NUMSPH)

    use pedra_ibtfun

#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_nuclei.h>
#include <pcm_symmet.h>
! f#include <pcm_trans.h>

!  reproduce the irreducibile part of the cavity
    
    integer :: numts, numsph
    real(8) :: vert(numts,10,3), centr(numts,10,3)
    integer :: nperm(numsph,*)

    real(8), parameter :: d0 = 0.0d0
    integer :: i, ii, ii2, isymop, k, l

!  replication of the cavity

! heckWRITE(LVPRI,*)'IN REPCAV_',NTS
    NTSIRR = NTS
    NTS = NTS * (MAXREP + 1)

    IF (NTS > MXTS) THEN
        WRITE(*,*) 'Errornumber 6 - check in the code'
        STOP
    END IF
    DO ISYMOP=1,MAXREP
        DO I=1,NTSIRR
            II=I+Ntsirr*ISYMOP
            AS(II)=AS(I)
            NVERT(II)=NVERT(I)
            ISPHE(II)=NPERM(ISPHE(I),ISYMOP+1)
            DO K=1,NVERT(I)
                DO L=1,3
                    VERT(II,K,L)  = &
                    PT(IBTAND(ISYMAX(L,1),ISYMOP)) * VERT(I,K,L)
                    CENTR(II,K,L) = &
                    PT(IBTAND(ISYMAX(L,1),ISYMOP)) * CENTR(I,K,L)
                ENDDO
            ENDDO
        ENDDO
    ENDDO

! f Moving the normal points in a safe location where they will not be overwritten

    DO I=1,NTSIRR
        II=I+NTSIRR
        II2=I+(MAXREP+1)*NTSIRR
        XTSCOR(II2)=XTSCOR(II)
        YTSCOR(II2)=YTSCOR(II)
        ZTSCOR(II2)=ZTSCOR(II)
    ENDDO

! f      write(lvpri,*) 'see if the pt are correct',nts
    DO ISYMOP=1,MAXREP
        DO I=1,NTSIRR
            II=I+NTSIRR*ISYMOP
            XTSCOR(II)=PT(IBTAND(ISYMAX(1,1),ISYMOP)) * XTSCOR(I)
            YTSCOR(II)=PT(IBTAND(ISYMAX(2,1),ISYMOP)) * YTSCOR(I)
            ZTSCOR(II)=PT(IBTAND(ISYMAX(3,1),ISYMOP)) * ZTSCOR(I)
        !            II=I+NTS*(ISYMOP+MAXREP+1)
        !            II2=I+NTS
        !            XTSCOR(II)=PT(IBTAND(ISYMAX(1,1),ISYMOP)) * XTSCOR(II2)
        !            YTSCOR(II)=PT(IBTAND(ISYMAX(2,1),ISYMOP)) * YTSCOR(II2)
        !            ZTSCOR(II)=PT(IBTAND(ISYMAX(3,1),ISYMOP)) * ZTSCOR(II2)
        ENDDO
    ENDDO
! f
!      do i=1,nts
!         do j=1,nts
!            distance=sqrt((xtscor(j)-xtscor(i))**2
!     $                   +(ytscor(j)-ytscor(i))**2
!     $                   +(ztscor(j)-ztscor(i))**2)
!            if(i.ne.j .and. distance.lt.1.0d-8) then
!               write(lvpri,*) 'tesserae too close:',i,j
!            endif
!         enddo
!      enddo
! f
    
    END SUBROUTINE REPCAV

    SUBROUTINE TESSERA(NS,NV,PTS,CCC,PP,PP1,AREA,INTSPH,NUMTS)

    use pedra_utils, only: around
    use pedra_print, only: output

#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>

    integer :: ns, nv, numts
    real(8) :: area
    real(8) :: pts(3,10), ccc(3,10), pp(3), pp1(3)
    integer :: intsph(numts,10)
    real(8) :: p1(3),p2(3),p3(3),p4(3),p4bis(3),point(3),p5(3),p6(3)
    real(8) :: pscr(3,10),cccp(3,10),pointl(3,10), distck(10,10)
    integer :: ind(10),ltyp(10),intscr(10),ntrhso(10)
    logical :: lan
    character(16) :: typlab_
    
    real(8) :: dcheck, de2, delr, delr2, diffdr, dist, dist1, dist2
    real(8) :: dnorm, rc, rc2, tol
    real(8) :: x1, x2, y1, y2, z1, z2
    integer :: i, j, ic, icop, icut, idx, idx2, ii, intcas, iprcav
    integer :: iv1, iv2, ivnew, ivold, jj, k, l, n, nsfe1, nvleft, nvnegl

#include <pcm_pcm.h>
#include <pcm_pcmlog.h>

!     Coord. del centro che sottende l`arco tra i vertici
!     n e n+1 (per i primi tre vertici e' sicuramente il centro della
!     sfera) e sfera alla cui intersezione con NS appartiene l'arco (se
!     appartiene alla sfera originaria INTSPH(numts,N)=NS)


    IPRCAV = 0
    LAN = .FALSE.
    AREA = 0.0D+00
    DO J=1, 3
        CCC(1,J) = XE(NS)
        CCC(2,J) = YE(NS)
        CCC(3,J) = ZE(NS)
    ENDDO
    IF(IPRCAV >= 10) THEN
        CALL AROUND('INPUT DATA IN TESSERA_', lvpri)
        WRITE(LVPRI,1000) NS,NV,NUMTS
        WRITE(LVPRI,*) '=======PTS======='
        CALL OUTPUT(PTS,1,3,1,3,3,3,1,LVPRI)
        WRITE(LVPRI,*) '=======CCC======='
        CALL OUTPUT(CCC,1,3,1,3,3,3,1,LVPRI)
    END IF

!     INTSPH viene riferito alla tessera -numts-, e in seguito riceve il
!     numero corretto.

    DO N = 1, 3
        INTSPH(NUMTS,N) = NS
    ENDDO

!     Loop sulle altre sfere

    DO 150 NSFE1=1,NESF
        IF(NSFE1 == NS) GO TO 150
        IF(IPRCAV >= 10) THEN
            WRITE(LVPRI,1005) NSFE1
            WRITE(LVPRI,1010) XE(NSFE1),YE(NSFE1),ZE(NSFE1),RE(NSFE1)
        END IF
    
    !     Memorizza i vertici e i centri che sottendono gli archi
    
        DO J =1, NV
            INTSCR(J) = INTSPH(NUMTS,J)
            DO I = 1,3
                PSCR(I,J) = PTS(I,J)
                CCCP(I,J) = CCC(I,J)
            ENDDO
        ENDDO
        IF(IPRCAV >= 10) THEN
            CALL AROUND('ACTUAL TESSERA_ STATUS', lvpri)
            WRITE(LVPRI,1000) NS,NV,NUMTS
            WRITE(LVPRI,*) '=======PSCR======='
            CALL OUTPUT(PSCR,1,3,1,NV,3,NV,1,LVPRI)
            WRITE(LVPRI,*) '=======CCCP======='
            CALL OUTPUT(CCCP,1,3,1,NV,3,NV,1,LVPRI)
        END IF
                 
    
        ICOP = 0
        DO J =1, 10
            IND(J) = 0
            LTYP(J) = 0
        ENDDO
    
    !     Loop sui vertici della tessera considerata
    
        DO 100 I=1,NV
            DELR2=(PTS(1,I)-XE(NSFE1))**2+(PTS(2,I)-YE(NSFE1))**2+ &
            (PTS(3,I)-ZE(NSFE1))**2
            DELR=SQRT(DELR2)
            IF (IPRPCM >= 10) THEN
                WRITE(LVPRI,1015) DELR,RE(NSFE1)
            ENDIF
            DIFFDR = DELR - RE(NSFE1)
            IF(DIFFDR < 0.0D0) THEN
                IND(I) = 1
                ICOP = ICOP+1
            END IF
        100 END DO
    !     Se la tessera e' completamente coperta, la trascura
        IF(ICOP == NV) THEN
            IF(IPRCAV >= 10) WRITE(LVPRI,1020) NSFE1
            RETURN
        !           ******
        END IF

    
    !     Controlla e classifica i lati della tessera: LTYP = 0 (coperto),
    !     1 (tagliato con il II vertice coperto), 2 (tagliato con il I
    !     vertice coperto), 3 (bitagliato), 4 (libero)
    !     Loop sui lati
    
        DO L = 1, NV
            IV1 = L
            IV2 = L+1
            IF(L == NV) IV2 = 1
            IF(IND(IV1) == 1 .AND. IND(IV2) == 1) THEN
                LTYP(L) = 0
            ELSE IF(IND(IV1) == 0 .AND. IND(IV2) == 1) THEN
                LTYP(L) = 1
            ELSE IF(IND(IV1) == 1 .AND. IND(IV2) == 0) THEN
                LTYP(L) = 2
            ELSE IF(IND(IV1) == 0 .AND. IND(IV2) == 0) THEN
                LTYP(L) = 4
                RC2 = (CCC(1,L)-PTS(1,L))**2 + (CCC(2,L)-PTS(2,L))**2 + &
                (CCC(3,L)-PTS(3,L))**2
                RC = SQRT(RC2)
            
            !     Su ogni lato si definiscono 11 punti equispaziati, che vengono
            !     controllati
            
                TOL = - 1.0D-10
                DO II = 1, 11
                    POINT(1) = PTS(1,IV1) + &
                    II * (PTS(1,IV2)-PTS(1,IV1)) / 11
                    POINT(2) = PTS(2,IV1) + &
                    II * (PTS(2,IV2)-PTS(2,IV1)) / 11
                    POINT(3) = PTS(3,IV1) + &
                    II * (PTS(3,IV2)-PTS(3,IV1)) / 11
                    POINT(1) = POINT(1) - CCC(1,L)
                    POINT(2) = POINT(2) - CCC(2,L)
                    POINT(3) = POINT(3) - CCC(3,L)
                    DNORM = SQRT(POINT(1)**2 + POINT(2)**2 + POINT(3)**2)
                    POINT(1) = POINT(1) * RC / DNORM + CCC(1,L)
                    POINT(2) = POINT(2) * RC / DNORM + CCC(2,L)
                    POINT(3) = POINT(3) * RC / DNORM + CCC(3,L)
                    DIST = SQRT( (POINT(1)-XE(NSFE1))**2 + &
                    (POINT(2)-YE(NSFE1))**2 + &
                    (POINT(3)-ZE(NSFE1))**2 )
                    IF((DIST - RE(NSFE1)) < TOL) THEN
                    !     IF(DIST.LT.RE(NSFE1)) then
                        LTYP(L) = 3
                        DO JJ = 1, 3
                            POINTL(JJ,L) = POINT(JJ)
                        ENDDO
                        GO TO 160
                    END IF
                ENDDO
            END IF
            160 CONTINUE
            IF(IPRCAV >= 10) THEN
                WRITE(LVPRI,1025) L,TYPLAB(LTYP(L))
            END IF
        ENDDO
    
    !     Se la tessera e' spezzata in due o piu' tronconi, la trascura
    
        ICUT = 0
        DO L = 1, NV
            IF(LTYP(L) == 1 .OR. LTYP(L) == 2) ICUT = ICUT + 1
            IF(LTYP(L) == 3) ICUT = ICUT + 2
        ENDDO
        ICUT = ICUT / 2
        IF(ICUT > 1) THEN
            WRITE(LVPRI,*) 'Tessera cut in pieces and removed.'
            RETURN
        END IF
    
    !     Creazione dei nuovi vertici e lati della tessera
    !     Loop sui lati
    
        N = 1
        IF(IPRCAV >= 10) WRITE(LVPRI,*) &
        'NOW CREATING NEW VERTICES AND EDGES IF NEED BE....'
        DO 300 L = 1, NV
        !     Se il lato L e' coperto:
            IF(LTYP(L) == 0) GO TO 300
            IV1 = L
            IV2 = L+1
            IF(L == NV) IV2 = 1
            IF (IPRCAV >= 10) THEN
                WRITE(LVPRI,1030) L,TYPLAB(LTYP(L)),LTYP(L)
                WRITE(LVPRI,1035) IV1,IV2
            ENDIF
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     Se il lato L e' tagliato (con il I vertice scoperto):
            IF(LTYP(L) == 1) THEN
            ! f store first point in the final set
                DO JJ = 1, 3
                    PTS(JJ,N) = PSCR(JJ,IV1)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                ENDDO
                INTSPH(NUMTS,N) = INTSCR(IV1)
                N = N+1
            
            !     Trova l'intersezione tra i due vertici del lato L
            
            !     P1 = coord. del primo vertice
            !     P2 = coord. del secondo vertice
            !     P3 = coord. del centro dell`arco sotteso
            !     P4 = coord. dell'intersezione
            
                DO JJ = 1, 3
                    P1(JJ) = PSCR(JJ,IV1)
                    P2(JJ) = PSCR(JJ,IV2)
                    P3(JJ) = CCCP(JJ,IV1)
                    P4(JJ) = 0.0D0
                ENDDO
                INTCAS=0
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1040) INTCAS,NSFE1
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                END IF
                CALL INTER(P1,P2,P3,P4,NSFE1,INTCAS)
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1050)
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                END IF
                DIST1 = SQRT((P4(1)-P1(1))**2 + (P4(2)-P1(2))**2 + &
                (P4(3)-P1(3))**2)
                DIST2 = SQRT((P4(1)-P2(1))**2 + (P4(2)-P2(2))**2 + &
                (P4(3)-P2(3))**2)
            !     Aggiorna i vertici della tessera e il centro dell'arco
                DO JJ = 1,3
                    PTS(JJ,N) = P4(JJ)
                ENDDO
            
            !     Il nuovo arco sara' sotteso tra questo e il prossimo punto
            !     di intersezione: il centro che lo sottende
            !     sara' il centro del cerchio di intersezione tra la sfera NS
            !     e la sfera NSFE1.
            
                DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+ &
                (ZE(NSFE1)-ZE(NS))**2
                CCC(1,N)=XE(NS)+(XE(NSFE1)-XE(NS))* &
                (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
                CCC(2,N)=YE(NS)+(YE(NSFE1)-YE(NS))* &
                (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
                CCC(3,N)=ZE(NS)+(ZE(NSFE1)-ZE(NS))* &
                (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
                INTSPH(NUMTS,N) = NSFE1
                N = N+1
            END IF
        
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     Se il lato L e' tagliato (con il II vertice scoperto):
            IF(LTYP(L) == 2) THEN
            !     Trova l'intersezione tra i due vertici del lato L
            
            !     P1 = coord. del primo vertice
            !     P2 = coord. del secondo vertice
            !     P3 = coord. del centro dell`arco sotteso
            !     P4 = coord. dell'intersezione
            
                DO JJ = 1, 3
                    P1(JJ) = PSCR(JJ,IV1)
                    P2(JJ) = PSCR(JJ,IV2)
                    P3(JJ) = CCCP(JJ,IV1)
                    P4(JJ) = 0.0D0
                ENDDO
                INTCAS=1
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1040) INTCAS,NSFE1
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                END IF
                CALL INTER(P1,P2,P3,P4,NSFE1,INTCAS)
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1050)
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                END IF
            !     Aggiorna i vertici della tessera e il centro dell'arco
                DO JJ = 1,3
                    PTS(JJ,N) = P4(JJ)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                ENDDO
                INTSPH(NUMTS,N) = INTSCR(IV1)
                N = N+1
            END IF
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     Se il lato e' intersecato due volte:
            IF(LTYP(L) == 3) THEN
                DO JJ = 1, 3
                    PTS(JJ,N) = PSCR(JJ,IV1)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                ENDDO
                INTSPH(NUMTS,N) = INTSCR(IV1)
                N = N+1
            
            !     Trova l'intersezione tra il primo vertice e un punto intermedio
            !     coperto
            
            !     P1 = coord. del primo vertice
            !     P2 = coord. del secondo vertice
            !     P3 = coord. del centro dell`arco sotteso
            !     P4 = coord. dell'intersezione
            
                DO JJ = 1, 3
                    P1(JJ) = PSCR(JJ,IV1)
                    P2(JJ) = POINTL(JJ,L)
                    P3(JJ) = CCCP(JJ,IV1)
                    P4(JJ) = 0.0D0
                ENDDO
                INTCAS=0
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1040) INTCAS,NSFE1
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                END IF
                CALL INTER(P1,P2,P3,P4,NSFE1,INTCAS)
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1050)
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                END IF
            !     Aggiorna i vertici della tessera e il centro dell'arco
                DO JJ = 1,3
                    PTS(JJ,N) = P4(JJ)
                ENDDO
            
            !     Il nuovo arco sara' sotteso tra questo e il prossimo punto
            !     di intersezione: il centro che lo sottende
            !     sara' il centro del cerchio di intersezione tra la sfera NS
            !     e la sfera NSFE1.
            
                DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+ &
                (ZE(NSFE1)-ZE(NS))**2
                CCC(1,N)=XE(NS)+(XE(NSFE1)-XE(NS))* &
                (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
                CCC(2,N)=YE(NS)+(YE(NSFE1)-YE(NS))* &
                (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
                CCC(3,N)=ZE(NS)+(ZE(NSFE1)-ZE(NS))* &
                (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
                INTSPH(NUMTS,N) = NSFE1
                N = N+1
            
            !     Trova l'intersezione tra un punto intermedio coperto e il
            !     secondo vertice
            
            !     P1 = coord. del primo vertice
            !     P2 = coord. del secondo vertice
            !     P3 = coord. del centro dell`arco sotteso
            !     P4 = coord. dell'intersezione
            
                DO JJ = 1, 3
                    P1(JJ) = POINTL(JJ,L)
                    P2(JJ) = PSCR(JJ,IV2)
                    P3(JJ) = CCCP(JJ,IV1)
                    P4(JJ) = 0.0D0
                ENDDO
                INTCAS=1
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1040) INTCAS,NSFE1
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                END IF
                CALL INTER(P1,P2,P3,P4,NSFE1,INTCAS)
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1050)
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                END IF
                DIST1 = SQRT((P4(1)-P1(1))**2 + (P4(2)-P1(2))**2 + &
                (P4(3)-P1(3))**2)
                DIST2 = SQRT((P4(1)-P2(1))**2 + (P4(2)-P2(2))**2 + &
                (P4(3)-P2(3))**2)
            !     Aggiorna il vertice e il centro dell'arco
                DO JJ = 1,3
                    PTS(JJ,N) = P4(JJ)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                ENDDO
                INTSPH(NUMTS,N) = INTSCR(IV1)
                N = N + 1
            END IF
        
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     Se il lato e' scoperto:
            IF(LTYP(L) == 4) THEN
                DO JJ = 1, 3
                    PTS(JJ,N) = PSCR(JJ,IV1)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                ENDDO
                INTSPH(NUMTS,N) = INTSCR(IV1)
                N = N+1
            END IF
        300 END DO
        NV = N - 1
        IF(IPRCAV >= 10) THEN
            CALL AROUND('AFTER INTER_SECTION TESSERA_ STATUS', lvpri)
            WRITE(LVPRI,1000) NS,NV,NUMTS
            WRITE(LVPRI,*) '=======PTS======='
            CALL OUTPUT(PTS,1,3,1,NV,3,NV,1,LVPRI)
            WRITE(LVPRI,*) '=======CCC======='
            CALL OUTPUT(CCC,1,3,1,NV,3,NV,1,LVPRI)
            DO IDX=1,NV
                IDX2=IDX+1
                IF (IDX2 > NV) IDX2=1
                DIST1=0.0D0
                DIST2=0.0D0
                DO IC=1,3
                    DIST1 = DIST1 + (PTS(IC,IDX ) - CCC(IC,IDX))**2
                    DIST2 = DIST2 + (PTS(IC,IDX2) - CCC(IC,IDX))**2
                END DO
                DCHECK=DIST1-DIST2
                WRITE(LVPRI,1055) IDX,DCHECK
            END DO
        END IF
    !     Controlla che il numero di vertici creati non sia eccessivo
        IF(NV > 10) THEN
            WRITE(LVPRI,*)'TOO MANY VERTICES IN TESSERA_: BYE BYE...'
            STOP
        END IF
    150 END DO
    IF(IPRCAV >= 10) THEN
        CALL AROUND('FINAL TESSERA_ STATUS', lvpri)
        WRITE(LVPRI,1000) NS,NV,NUMTS
        WRITE(LVPRI,*) '=======PSCR======='
        CALL OUTPUT(PSCR,1,3,1,NV,3,NV,1,LVPRI)
        WRITE(LVPRI,*) '=======CCCP======='
        CALL OUTPUT(CCCP,1,3,1,NV,3,NV,1,LVPRI)
    END IF

!     Se la tessera non e' stata scartata, a questo punto ne troviamo
!     l'area e il punto rappresentativo


! f We label the edges whose lenght is lower than a certain threshold

    NVNEGL = 0
    DO I=1,NV
        NTRHSO(I) = 0
        II=I+1
        IF (I == NV) II = 1
        X1 = PTS(1,I)
        Y1 = PTS(2,I)
        Z1 = PTS(3,I)
        X2 = PTS(1,II)
        Y2 = PTS(2,II)
        Z2 = PTS(3,II)
        DIST = SQRT((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)
        IF (DIST < 1.0D-6) THEN
            NTRHSO(I) = 1
            NVNEGL = NVNEGL + 1
        ! f            NTRHSO(I) = 0
        END IF
    ENDDO
    NVLEFT = NV - NVNEGL

! If we end up with a tessera with less than three long enough edges
! we discard it.

    IF (NVLEFT < 3) THEN
        WRITE(LVPRI,*) '*******WARNING*******\\', &
        ' NEGLECTED TESSERA_: TOO FEW ACCEPTABLE EDGES'
        RETURN
    
    ! Removal of too small edges to avoid numerical problems
    ! in area calculation (done only if one ore more edges are
    ! to be neglected).
    
    ELSE IF (NVNEGL > 0) THEN
        IVNEW = 0
        DO IVOLD = 1,NV
            IF(NTRHSO(IVOLD) == 0) THEN
                IVNEW = IVNEW + 1
                INTSPH(NUMTS,IVNEW) = INTSPH(NUMTS,IVOLD)
                DO K =1,3
                    CCC(K,IVNEW) = CCC(K,IVOLD)
                    PTS(K,IVNEW) = PTS(K,IVOLD)
                ENDDO
            END IF
        END DO
        IF (IVNEW /= NVLEFT) THEN
            WRITE(LVPRI,*) 'SIRCAV: BADLY MADE ALGORITHM!'
            STOP
        ENDIF
        NV = NVLEFT
    END IF

! Finally we calculate area and center of the tessera!

    CALL GAUBON(NV,NS,PTS,CCC,PP,PP1,AREA,INTSPH,NUMTS)
    RETURN
    1000 FORMAT('NS=',I4,' NV=',I4,' NUMTS=', I4)
    1005 FORMAT('CHECKING INTER_SECTIONS WITH SPHERE N.',I4)
    1010 FORMAT('Coordinates: X=',F12.8,'  Y=',F12.8,'  Z=',F12.8, &
    ' R=',F12.8)
    1015 FORMAT('VERT-OTHER CENTER DIST:',F10.8,' SPHERE RADIUS',F10.8)
    1020 FORMAT('THIS TESSERA_ IS COMPLETELY COVERED BY SPHERE N.',I4)
    1025 FORMAT('EDGE N.',I4,' STATUS: ',A18)
    1030 FORMAT('EDGE N.',I2,' IS IN STATUS ',A16,' CASE N.',I2)
    1035 FORMAT('VERTICES NUMBERS:',I3,' AND',I3)
    1040 FORMAT('CALLING INTER_, CASE N.',I2,' INTERSECTING SPHERE:',I4)
    1045 FORMAT('P',I1,' X = ',F12.9,' Y = ',F12.9,' Z = ',F12.9)
    1050 FORMAT('AFTER INTER_')
    1055 FORMAT('DIST CHECK EDGE N.',I2,' DIFFERENCE:',F12.10)

    END SUBROUTINE TESSERA
    
    SUBROUTINE INTER(P1,P2,P3,P4,NS,I)

#include <pcm_priunit.h>
#include <pcm_mxcent.h>
#include <pcm_pcmdef.h>

    real(8) :: p1(3), p2(3), p3(3), p4(3)
    integer :: ns, i
    logical :: lalow,lblow

    real(8) :: alphat, delta, diff, diff2, diff2a, diff2b
    real(8) :: diffa, diffb, dnorm, p1p3, p2p3, r, r2, tol
    integer :: j, jj, m
    integer :: iprcav = 0

#include <pcm_pcm.h>
#include <pcm_pcmlog.h>

!     Trova il punto P4, sull`arco P1-P2 sotteso dal centro P3, che
!     si trova sulla superficie della sfera NS
!     P4 e' definito come combinazione lineare di P1 e P2, con
!     il parametro ALPHA ottimizzato per tentativi.

    R2 = (P1(1)-P3(1))**2+(P1(2)-P3(2))**2+(P1(3)-P3(3))**2
    R = SQRT(R2)
! f Don't change this threshold unless you know exactly what you are doing.
    TOL = 1.0D-12
    ALPHAT = 0.5D+00
    DELTA = 0.0D+00
! f Check distances of P1 and P2 from P3
    P1P3 = SQRT((P1(1)-P3(1))**2+(P1(2)-P3(2))**2+(P1(3)-P3(3))**2)
    P2P3 = SQRT((P2(1)-P3(1))**2+(P2(2)-P3(2))**2+(P2(3)-P3(3))**2)
    IF(IPRCAV >= 10) THEN
        WRITE(LVPRI,1000) P1P3-P2P3
    END IF
! f Check if one or both starting points are below the threshold
    LALOW=.FALSE.
    LBLOW=.FALSE.
    DIFF2A= (P1(1)-XE(NS))**2 + (P1(2)-YE(NS))**2 + (P1(3)-ZE(NS))**2
    DIFFA = ABS(SQRT(DIFF2A) - RE(NS))
    DIFF2B= (P2(1)-XE(NS))**2 + (P2(2)-YE(NS))**2 + (P2(3)-ZE(NS))**2
    DIFFB = ABS(SQRT(DIFF2B) - RE(NS))
    LALOW=(DIFFA.LT.TOL)
    LBLOW=(DIFFB.LT.TOL)
    IF(LALOW .AND. .NOT. LBLOW) THEN
        IF(IPRCAV >= 10) WRITE (LVPRI,*) &
        'INTER_: TAKEN FIRST POINT'
        DO J=1,3
            P4(J)=P1(J)
        END DO
        GOTO 100
    END IF
    IF(LBLOW .AND. .NOT. LALOW) THEN
        IF(IPRCAV >= 10) WRITE (LVPRI,*) &
        'INTER_: TAKEN SECOND POINT'
        DO J=1,3
            P4(J)=P2(J)
        END DO
        GOTO 100
    END IF
    IF(LALOW .AND. LBLOW) THEN
        IF(DIFFA <= DIFFB) THEN
            IF(IPRCAV >= 10) WRITE (LVPRI,*) &
            'INTER_: TAKEN FIRST POINT'
            DO J=1,3
                P4(J)=P1(J)
            END DO
        ELSE
            IF(IPRCAV >= 10) WRITE (LVPRI,*) &
            'INTER_: TAKEN SECOND POINT'
            DO J=1,3
                P4(J)=P2(J)
            END DO
        END IF
        GOTO 100
    END IF

!     Start iterations

    M = 1
    10 CONTINUE
    ALPHAT = ALPHAT + DELTA
    DNORM = 0.0D+00
    DO JJ = 1,3
        P4(JJ)=P1(JJ)+ALPHAT*(P2(JJ)-P1(JJ))-P3(JJ)
        DNORM = DNORM + P4(JJ)**2
    ENDDO
    DNORM = SQRT(DNORM)
    DO JJ = 1,3
        P4(JJ)= P4(JJ)*R/DNORM + P3(JJ)
    ENDDO
    DIFF2=(P4(1)-XE(NS))**2 + (P4(2)-YE(NS))**2 + (P4(3)-ZE(NS))**2
    DIFF = SQRT(DIFF2) - RE(NS)
    IF(ABS(DIFF) < TOL) GOTO 100
    IF(I == 0) THEN
        IF(DIFF > 0.0D+00) DELTA = 1.0D+00/(2.0D+00**DBLE(M+1))
        IF(DIFF < 0.0D+00) DELTA = - 1.0D+00/(2.0D+00**DBLE(M+1))
        M = M + 1
    ELSE IF(I == 1) THEN
        IF(DIFF > 0.0D+00) DELTA = - 1.0D+00/(2.0D+00**DBLE(M+1))
        IF(DIFF < 0.0D+00) DELTA = 1.0D+00/(2.0D+00**DBLE(M+1))
        M = M + 1
    ELSE
        WRITE(*,*) '7'
        STOP
    END IF
    IF (M > 300)THEN
        WRITE(LVPRI,*)'Too many iterations in INTER_! BYE BYE ...'
        WRITE(LVPRI,*)'P1',P1(1),P1(2),P1(3),DIFFA,LALOW
        WRITE(LVPRI,*)'P2',P2(1),P2(2),P2(3),DIFFB,LBLOW
        WRITE(LVPRI,*)'P3',P3(1),P3(2),P3(3)
        WRITE(LVPRI,*)'P4',P4(1),P4(2),P4(3),DIFF,I
        WRITE(LVPRI,*)'SPHERE',XE(NS),YE(NS),ZE(NS),RE(NS)
        WRITE(*,*) '8'
        STOP
    END IF
    GO TO 10

! Final printing and return

    100 CONTINUE
    IF(IPRPCM >= 10) THEN
        WRITE(LVPRI,1005)
        WRITE(LVPRI,1010) 1,P1(1),P1(2),P1(3)
        WRITE(LVPRI,1010) 2,P2(1),P2(2),P2(3)
        WRITE(LVPRI,1010) 3,P3(1),P3(2),P3(3)
        WRITE(LVPRI,1010) 4,P4(1),P4(2),P4(3)
        WRITE(LVPRI,1015) XE(NS),YE(NS),ZE(NS),RE(NS)
    END IF
    RETURN

    1000 FORMAT('INTER_, distance consistency check: ',F14.12)
    1005 FORMAT(/'Final result from INTER_!')
    1010 FORMAT('P',I1,' X = ',F12.9,' Y = ',F12.9,' Z = ',F12.9)
    1015 FORMAT('SPHERE',' X = ',F12.9,' Y = ',F12.9,' Z = ',F12.9, &
    ' R = ',F12.9)

    END SUBROUTINE INTER
    
    SUBROUTINE GAUBON(NV,NS,PTS,CCC,PP,PP1,AREA,INTSPH,NUMTS)

    use pedra_dblas, only: vector_product

#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>

    integer :: nv, ns, numts
    real(8) :: area
    real(8) :: pts(3, 10), ccc(3, 10), pp(3), pp1(3), beta(10)
    integer :: intsph(numts, 10), ntrhso(10)
    real(8) :: p1(3),p2(3),p3(3),u1(3),u2(3),phin(10),weight(0:10)
    
    real(8), parameter :: d0 = 0.0d0
    real(8), parameter :: dp5 = 0.5d0
    real(8), parameter :: d1 = 1.0d0
    real(8), parameter :: pi = acos(-1.0d0)

    real(8) :: cosphin, costn, sum1, sum2, sumphi, dnorm1, dnorm2, tpi
    real(8) :: scal, dnorm, dnorm3
    real(8) :: x1, x2, y1, y2, z1, z2
    integer :: i, jj, n, n0, n1, n2, nsfe1

#include <pcm_pcm.h>
#include <pcm_pcmlog.h>

!     Sfrutta il teorema di Gauss-Bonnet per calcolare l'area
!     della tessera con vertici PTS(3,NV). Consideriamo sempre
!     che il lato N della tessera e' quello compreso tra i vertici
!     N e N+1 (oppure NV e 1). In CCC(3,NV) sono le posizioni dei
!     centri degli archi che sottendono i vari lati della tessera.
!     La formula di Gauss-Bonet per le sfere e':
!            Area=R^2[2pi+S(Phi(N)cosT(N))-S(Beta(N))]
!     dove Phi(N) e' la lunghezza d'arco (in radianti) del lato N,
!     T(N) e' l'angolo polare del lato N, Beta(N) l'angolo esterno
!     relativo al vertice N.

    TPI=2*PI

!     Calcola la prima sommatoria
    SUM1 = D0
    SUMPHI = D0
    DO 100 N = 1, NV
        PHIN(N) = 0.0D0
    
    !         IF (NTRHSO(N) .EQ. 1) GOTO 100
    
        X1 = PTS(1,N) - CCC(1,N)
        Y1 = PTS(2,N) - CCC(2,N)
        Z1 = PTS(3,N) - CCC(3,N)
        IF(N < NV) THEN
            X2 = PTS(1,N+1) - CCC(1,N)
            Y2 = PTS(2,N+1) - CCC(2,N)
            Z2 = PTS(3,N+1) - CCC(3,N)
        ELSE
            X2 = PTS(1,1) - CCC(1,N)
            Y2 = PTS(2,1) - CCC(2,N)
            Z2 = PTS(3,1) - CCC(3,N)
        END IF
        DNORM1 = X1*X1 + Y1*Y1 + Z1*Z1
        DNORM2 = X2*X2 + Y2*Y2 + Z2*Z2
        SCAL = X1*X2 + Y1*Y2 + Z1*Z2
        COSPHIN = SCAL / (SQRT(DNORM1*DNORM2))
        IF(COSPHIN > 1.0D+00) COSPHIN = 1.0D+00
        PHIN(N) = ACOS(COSPHIN)
        SUMPHI = SUMPHI + PHIN(N)
    
    !     NSFE1 e' la sfera con cui la sfera NS si interseca (eventualmente)
        NSFE1 = INTSPH(NUMTS,N)
        X1 = XE(NSFE1) - XE(NS)
        Y1 = YE(NSFE1) - YE(NS)
        Z1 = ZE(NSFE1) - ZE(NS)
        DNORM1 = SQRT(X1*X1 + Y1*Y1 + Z1*Z1)
        IF(DNORM1 == D0) DNORM1 = 1.0D+00
        X2 = PTS(1,N) - XE(NS)
        Y2 = PTS(2,N) - YE(NS)
        Z2 = PTS(3,N) - ZE(NS)
        DNORM2 = SQRT(X2*X2 + Y2*Y2 + Z2*Z2)
        COSTN = (X1*X2+Y1*Y2+Z1*Z2)/(DNORM1*DNORM2)
        SUM1 = SUM1 + PHIN(N) * COSTN
    100 END DO

! WEIGHTS GENERATION FOR CENTER DEFINITION

! the use of old weights is only for backward compatibility
! the new weights give a smooth variation of the
! tessera change during geometry optimizations especially when
! the number of vertices of a tessera changes.

! the difference in the sum has no importance since at the end
! the point is renormalized in order to be on the sphere

    IF (OLDCEN) THEN
        DO I = 1,NV
            WEIGHT(I) = 1.0D0
        ENDDO
    ELSE
        DO I = 1, NV
            WEIGHT(I) = PHIN(I)
        END DO
        WEIGHT(0) = WEIGHT(NV)
        DO I = NV,1,-1
            WEIGHT(I) = WEIGHT(I) + WEIGHT(I-1)
        END DO
    END IF

!     Calcola la seconda sommatoria: l'angolo esterno Beta(N) e'
!     definito usando i versori (u(N-1),u(N)) tangenti alla sfera
!     nel vertice N lungo le direzioni dei lati N-1 e N:
!                cos( Pi-Beta(N) )=u(N-1)*u(N)
!            u(N-1) = [V(N) x (V(N) x V(N-1))]/NORM
!            u(N) = [V(N) x (V(N) x V(N+1))]/NORM
!     dove V(I) e' il vettore posizione del vertice I RISPETTO AL
!     CENTRO DELL'ARCO CHE SI STA CONSIDERANDO.

    SUM2 = 0.0D+00
!     Loop sui vertici
    DO 200 N = 1, NV
        DO JJ = 1, 3
            P1(JJ) = 0.0D+00
            P2(JJ) = 0.0D+00
            P3(JJ) = 0.0D+00
        ENDDO
        N1 = N
    ! f
    !         IF(N.GT.1) N0 = N - 1
    !         IF(N.EQ.1) N0 = NV
    !         IF(N.LT.NV) N2 = N + 1
    !         IF(N.EQ.NV) N2 = 1
    ! f
    
    ! If one or both the edges are labelled "small we jump to the preavious and/or next and so on
    ! until we reach a long enough edge.
    
        N0 = N
    ! f 220     CONTINUE
        N0 = MOD(NV+N0-1,NV)
        IF(N0 == 0) N0 = NV
    ! f         write(lvpri,*) "N0 is", N0, NTRHSO(N0)
    ! f         IF(NTRHSO(N0).EQ.1) GOTO 220
        N2 = N
    ! f         NCONT = 0
    ! f 230     CONTINUE
        N2 = MOD(N2+1,NV)
        IF(N2 == 0) N2 = NV
    ! f         write(lvpri,*) "N2 is", N2, NTRHSO(N1)
    ! f         IF(NTRHSO(N1 + NCONT).EQ.1) THEN
    ! f            NCONT = NCONT + 1
    ! f            GOTO 230
    ! F         END IF
    
    !     Trova i vettori posizione rispetto ai centri corrispondenti
    !     e i versori tangenti
    
    !     Lato N0-N1:
        DO JJ = 1, 3
            P1(JJ) = PTS(JJ,N1) - CCC(JJ,N0)
            P2(JJ) = PTS(JJ,N0) - CCC(JJ,N0)
        ENDDO
    
        CALL vector_product(P1,P2,P3,DNORM3)
        DO JJ = 1, 3
            P2(JJ) = P3(JJ)
        ENDDO
        CALL vector_product(P1,P2,P3,DNORM3)
        DO JJ = 1, 3
            U1(JJ) = P3(JJ)/DNORM3
        ENDDO
    
    !     Lato N1-N2:
        DO JJ = 1, 3
            P1(JJ) = PTS(JJ,N1) - CCC(JJ,N1)
            P2(JJ) = PTS(JJ,N2) - CCC(JJ,N1)
        ENDDO
    
        CALL vector_product(P1,P2,P3,DNORM3)
        DO JJ = 1, 3
            P2(JJ) = P3(JJ)
        ENDDO
        CALL vector_product(P1,P2,P3,DNORM3)
        DO JJ = 1, 3
            U2(JJ) = P3(JJ)/DNORM3
        ENDDO
    
        BETA(N) = ACOS(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
    ! F         SUM2 = SUM2 + (PI - BETAN)
    200 END DO
    do I=1,NV
    ! F         II = I - 1
    ! f         IF (II .EQ. 0) II = NV
    ! f         SUM2 = SUM2 - BETA(I) * (D1 - DP5 * DBLE(NTRHSO(I)+NTRHSO(II)))
    ! f     $               + PI * (D1 - DBLE(NTRHSO(I)))
        SUM2 = SUM2 - BETA(I)
    enddo
!     Calcola l'area della tessera
    AREA = RE(NS)*RE(NS)*(DBLE(2-NV) * PI + SUM1 - SUM2)
!     Trova il punto rappresentativo (come media dei vertici)
    DO JJ = 1, 3
        PP(JJ) = 0.0D+00
    ENDDO
    DO I = 1, NV
        PP(1) = PP(1) + (PTS(1,I)-XE(NS)) * WEIGHT(I)
        PP(2) = PP(2) + (PTS(2,I)-YE(NS)) * WEIGHT(I)
        PP(3) = PP(3) + (PTS(3,I)-ZE(NS)) * WEIGHT(I)
    END DO
    DNORM = 0.0D+00
    DO JJ = 1, 3
        DNORM = DNORM + PP(JJ)*PP(JJ)
    ENDDO
    PP(1) = XE(NS) + PP(1) * RE(NS) / SQRT(DNORM)
    PP(2) = YE(NS) + PP(2) * RE(NS) / SQRT(DNORM)
    PP(3) = ZE(NS) + PP(3) * RE(NS) / SQRT(DNORM)
!     Trova il punto sulla normale (interna!) distante DR dal punto
!     rappresentativo
    PP1(1) = XE(NS) + (PP(1) - XE(NS)) * (RE(NS) - DR) / RE(NS)
    PP1(2) = YE(NS) + (PP(2) - YE(NS)) * (RE(NS) - DR) / RE(NS)
    PP1(3) = ZE(NS) + (PP(3) - ZE(NS)) * (RE(NS) - DR) / RE(NS)

!     A causa delle approssimazioni numeriche, l'area di alcune piccole
!     tessere puo' risultare negativa, e viene in questo caso trascurata
    IF(AREA < 0.0D+00) THEN
        WRITE(LVPRI,1000) NS,area
        AREA = 0.0D+00
        1000 FORMAT(/,'WARNING: THE AEREA OF ONE TESSERA_ ON SPHERE ',I3, &
        ' IS IGNORED',F16.10)
    END IF
    
    END SUBROUTINE GAUBON
    
    SUBROUTINE CAVSPL(ICAV1,ICAV2,NCAV1,NCAV2,NATM,SOME)

#include <pcm_priunit.h>
#include <pcm_mxcent.h>
#include <pcm_pcmdef.h>

    integer :: icav1(mxcent), icav2(mxcent)
    integer :: ncav1, ncav2, natm
    logical :: some

    integer :: i, n1, n2, nn, icen, n
    real(8) :: r, rr, sum, x, y, z, xx, yy, zz

#include <pcm_pcm.h>
#include <pcm_pcmlog.h>

!  QUESTA ROUTINE CONTROLLA SE IL SOLUTO E' CONTENUTO IN UNA
!  UNICA CAVITA' O IN CAVITA' DISTINTE: IN TAL CASO I NUMERI
!  D'ORDINE DELLE SFERE CHE COSTITUISCONO LA PRIMA E LA SECONDA
!  CAVITA' VENGONO MEMORIZZATI RISPETTIVAMENTE IN ICAV1 E ICAV2

    DO I=1,MXCENT
        ICAV1(I)=0
        ICAV2(I)=0
    ENDDO
    NCAV1=0
    NCAV2=0
    N1=1
    N2=1
    NN=1
    ICAV1(N1)=1
    N1=N1+1
    50 DO I=1,MXCENT
        ICAV2(I)=0
    ENDDO
    N2=1
    N=ICAV1(NN)
    DO 100 ICEN=1,NESF
        IF(ICEN == N) GO TO 100
        DO I=1,NESF
            IF(ICEN == ICAV1(I)) GO TO 100
        ENDDO
        X=XE(N)-XE(ICEN)
        XX=X*X
        Y=YE(N)-YE(ICEN)
        YY=Y*Y
        Z=ZE(N)-ZE(ICEN)
        ZZ=Z*Z
        RR=XX+YY+ZZ
        R=SQRT(RR)
        SUM=RE(N)+RE(ICEN)
        IF(SUM > R) THEN
            ICAV1(N1)=ICEN
            N1=N1+1
        ELSE
            ICAV2(N2)=ICEN
            N2=N2+1
        END IF
    100 END DO
    NN=NN+1
    IF(ICAV1(NN) /= 0) GO TO 50
    NCAV1=NN-1
    NCAV2=NESF-NCAV1
    IF(SOME) THEN
        IF(NCAV2 == 0) THEN
            WRITE(LVPRI,200)
        ELSE
            WRITE(LVPRI,300) NCAV1,NCAV2
            WRITE(LVPRI,400)
            WRITE(LVPRI,*) (ICAV1(I),I=1,NCAV1)
            WRITE(LVPRI,500)
            WRITE(LVPRI,*) (ICAV2(I),I=1,NCAV2)
        END IF
    END IF
    RETURN

    200 FORMAT(/10X,'THE SOLUTE IS ENCLOSED IN ONE CAVITY')
    300 FORMAT(/10X,'THE SOLUTE IS ENCLOSED IN TWO DISTINCT CAVITIES'/ &
    10X,'OF',I3,3X,'E',I3,3X,'SPHERE(S),  RESPECTIVELY')
    400 FORMAT(/10X,'THE FIRST CAVITY IS FORMED BY SPHERE(S) :'/)
    500 FORMAT(/10X,'THE SECOND CAVITY IS FORMED BY SPHERE(S) :'/)

    END SUBROUTINE CAVSPL
    
    SUBROUTINE PLOTCAV(Vert,NUMTS)

#include <pcm_priunit.h>
#include <pcm_iratdef.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>

    integer :: numts
    real(8) :: VERT(NUMTS,10,*)
    integer :: IVTS(MXTS,10)
    logical :: cavity_file_exists

    real(8) :: c1, c2, c3
    integer :: n, numv, i, j, k, last, lucav, index
    integer :: jcord

!     Prepare the input file for GeomView (coloured polyhedra)

    IF (ICESPH == 0 .OR. ICESPH == 2 .OR. ICESPH == 3) INDEX = 1
    LUCAV = 12121201
    inquire(file = 'cavity.off', exist = cavity_file_exists)
    if (cavity_file_exists) then
        open(lucav, &
        file = 'cavity.off', &
        status = 'old', &
        form = 'formatted', &
        access = 'sequential')
        close(lucav, status = 'delete')
    end if
    open(lucav, &
    file = 'cavity.off', &
    status = 'new', &
    form = 'formatted', &
    access = 'sequential')
    rewind(lucav)
    Numv = 0
    DO 10 I = 1, NTs
        NUMV = NUMV + NVert(i)
    10 END DO
    WRITE(LUCAV,500)'COFF'
    WRITE(LUCAV,1000)NUMV, NTs, NumV
    K = 0
    LAST = 0
    DO 2010 I = 1, NTs
        N=ISphe(I)
        IF(N /= Last) WRITE(LUCAV,1500)N
        Last = N
        IF(Index == 1) THEN
            Call COLOUR(N,C1,C2,C3)
        ELSE
            C1=1.0D0
            C2=1.0D0
            C3=1.0D0
        END IF
        Do 2020 j = 1, NVert(i)
            IVTS(I,J) = K
            K = K + 1
            WRITE(LUCAV,2001)(VERT(I,J,JCORD),JCORD=1,3), &
            C1,C2,C3,0.75,I
        2020 END DO
    2010 END DO
    Do 20 I = 1, NTs
        WRITE(LUCAV,3000)NVert(I),(IVTS(I,J),J=1,NVERT(I))
    20 END DO
    500 FORMAT(1x,a)
    1000 FORMAT(3i10)
    1500 FORMAT('# Sphere number ',i4)
    2001 FORMAT('  ',3f16.9,4f5.2,' # Tess. ',i4)
    3000 FORMAT('  ',14i10)
    
    END SUBROUTINE PLOTCAV
    
    SUBROUTINE COLOUR(N,C1,C2,C3)
#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_nuclei.h>
#include <pcm_pcmnuclei.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
    
    integer :: n
    real(8) :: c1, c2, c3

    Character(20) :: Col
    Data Col /'                    '/
    real(8) :: delta, diff, ddiff
    integer :: num, i

!     Assign tesserae colours for GeomView:
!     Carbon: green, nitrogen: blue, oxygen: red, hydrogen: light blue,
!     others: violet, added spheres: gray.

    Delta=1.d-03
    IF(N > NESFP)THEN
        Col = 'Gray'
        Call COLTSS(Col,C1,C2,C3)
        Return
    END IF
    NUM = 0
    DO I = 1, NUCDEP
        DIFF=( PCMCORD(1,I) - XE(N) )**2 &
        +( PCMCORD(2,I) - YE(N) )**2 &
        +( PCMCORD(3,I) - ZE(N) )**2
        DDIFF=Sqrt(DIFF)
        IF(DDIFF <= Delta) THEN
            NUM=INT(PCMCHG(I))
            GOTO 11
        END IF
    END DO
    11 CONTINUE
    IF(NUM == 6)THEN
        COL = 'Green'
    ELSEIF(NUM == 7)THEN
        COL = 'Blue'
    ELSEIF(NUM == 8)THEN
        COL = 'Red'
    ELSEIF(NUM == 1)THEN
        COL = 'Light Blue'
    ELSE
        COL = 'Fuchsia'
    ENDIF
    CALL COLTSS(Col,C1,C2,C3)
    
    end subroutine colour
    
    subroutine coltss(colour_,c1,c2,c3)

!     Assigne tesserae colours for GeomView:

    character(20) :: COLOUR_
    real(8) :: c1, c2, c3

    IF(COLOUR_ == 'White') THEN
        C1 = 1.0
        C2 = 1.0
        C3 = 1.0
    ELSE IF(COLOUR_ == 'Gray') THEN
        C1 = 0.66
        C2 = 0.66
        C3 = 0.66
    ELSE IF(COLOUR_ == 'Blue' .OR. COLOUR_ == 'Dark Blue') THEN
        C1 = 0.0
        C2 = 0.0
        C3 = 1.0
    ELSE IF(COLOUR_ == 'Light Blue') THEN
        C1 = 0.0
        C2 = 1.0
        C3 = 1.0
    ELSE IF(COLOUR_ == 'Green') THEN
        C1 = 0.0
        C2 = 1.0
        C3 = 0.0
    ELSE IF(COLOUR_ == 'Yellow') THEN
        C1 = 1.0
        C2 = 1.0
        C3 = 0.0
    ELSE IF(COLOUR_ == 'Orange') THEN
        C1 = 1.0
        C2 = 0.5
        C3 = 0.0
    ELSE IF(COLOUR_ == 'Violet') THEN
        C1 = 0.6
        C2 = 0.0
        C3 = 1.0
    ELSE IF(COLOUR_ == 'Pink' .OR. COLOUR_ == 'Light Red') THEN
        C1 = 1.0
        C2 = 0.5
        C3 = 1.0
    ELSE IF(COLOUR_ == 'Fuchsia') THEN
        C1 = 1.0
        C2 = 0.0
        C3 = 1.0
    ELSE IF(COLOUR_ == 'Red' .OR. COLOUR_ == 'Dark Red') THEN
        C1 = 1.0
        C2 = 0.0
        C3 = 0.0
    ELSE IF(COLOUR_ == 'Black') THEN
        C1 = 0.0
        C2 = 0.0
        C3 = 0.0
    ELSE
        WRITE(*,*) '8'
        STOP
    ENDIF
    
    end subroutine coltss
    
    subroutine prerep(nv,nt,its,cv,jtr,nvert,numts)

    use pedra_ibtfun

#include <pcm_priunit.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_mxcent.h>
#include <pcm_symmet.h>
#include <pcm_pgroup.h>

    integer :: nv, nt, its, nvert, numts
    integer :: jtr(numts, *)
    real(8) :: cv(nvert, *)
    logical :: lsymop(0:7)

    integer :: i, j, isymop, ii, jj, k

!      SYMOP(0) = ' E '
!      SYMOP(1) = 'Oyz'
!      SYMOP(2) = 'Oxz'
!      SYMOP(3) = 'C2z'
!      SYMOP(4) = 'Oxy'
!      SYMOP(5) = 'C2y'
!      SYMOP(6) = 'C2x'
!      SYMOP(7) = ' i '
    DO I = 0, 7
        LSYMOP(I) = .FALSE.
    END DO
    LSYMOP(0) = .TRUE.
    IF(GROUP == 'C1 ') THEN
        LSYMOP(1) = .TRUE.
        LSYMOP(2) = .TRUE.
        LSYMOP(4) = .TRUE.
    ELSEIF(GROUP == 'C2 ' .OR. GROUP == 'Ci ') THEN
        LSYMOP(1) = .TRUE.
        LSYMOP(2) = .TRUE.
    ELSEIF(GROUP == 'Cs ') THEN
        LSYMOP(2) = .TRUE.
        LSYMOP(4) = .TRUE.
    ELSEIF(GROUP == 'D2 ' .OR. GROUP == 'C2v') THEN
        LSYMOP(1) = .TRUE.
    ELSEIF(GROUP == 'C2h') THEN
        LSYMOP(2) = .TRUE.
    ELSEIF(GROUP /= 'D2h') THEN
        WRITE(*,*) 'Check symmetry group.'
        STOP
    ENDIF
    DO ISYMOP = 1,7
        IF(LSYMOP(ISYMOP)) THEN
        
        !     riproduzione vertici
        
            DO I = 1,NV
                II = I + NV
                DO K = 1,3
                    CV(II,K) = &
                    PT(IBTAND(ISYMOP,2**(K-1))) * CV(I,K)
                ENDDO
            ENDDO
        ! f            WRITE(LVPRI,*) 'I,II,JJ,K,JTR(II,K),JTR(I,K)'
        
        !     riproduzione topologia
        
            DO I = 1,NT
                II = I + NT
                JJ = NV
                Do K = 1,3
                    JTR(II,K) = JTR(I,K) + JJ
                ! f                  WRITE(LVPRI,*) I,II,JJ,K,JTR(II,K),JTR(I,K)
                ENDDO
            ENDDO
        
        !  aggiornamento indici
        
            NT = NT * 2
            NV = NV * 2
            ITS = ITS * 2
        ! f            write(lvpri,*) 'nt,nv,its',nt,nv,its
        ENDIF
    END DO
    DO I = 1,ITS
        DO J = 1,3
        ! f            WRITE(LVPRI,*) I,JTR(I,J),(CV(JTR(I,J),K),K=1,3)
        ENDDO
    ENDDO
    
    end subroutine prerep
    
    subroutine pcmtns(vmat, geom, amass, katom)

    use pedra_utils, only: wlkdin                
    use pedra_ibtfun

#include <pcm_priunit.h>
#include <pcm_mxcent.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_iratdef.h>
#include <pcm_pcmdef.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_symmet.h>
#include <pcm_nuclei.h>
#include <pcm_orgcom.h>

    integer :: katom
    real(8) :: vmat(3, 3), geom(katom, *), amass(*)

    real(8) :: eigval(3), eigvec(3, 3), tinert(3, 3)
    real(8) :: angmom(3), omegad(3), eiginv(3, 3), scal(3)
    integer :: iax(6), norder(3)
    logical :: planar, linear
    integer :: jatom, i, j, k, jax, nmax
    integer :: nopax, nshift

    real(8) :: dij
    
    iax = (/1,2,3,3,2,1/)

    jatom = 1

! As correctly pointed out by Kenneth, we want two molecules
! like HCCD and HCCH to have the same cavity so NUMIS is always 1

!      NUMIS = 1
!      DO 100 IATOM = 1, NUCIND
!         NATTYP = NINT(CHARGE(IATOM))
! f         NUMIS  = ISOTOP(IATOM)
!         DO 110 ISYMOP = 0, MAXOPR
!            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
!               AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')
!               DO 120 ICOOR = 1, 3
!                  GEOM(JATOM,ICOOR) = PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
!     &                               *CORD(ICOOR,IATOM) - CMXYZ(ICOOR)
! 120           CONTINUE
! f               WRITE(LVPRI,*) 'NATTYP, NUMIS,JATOM,AMASS',
! f     $              NATTYP, NUMIS,JATOM,AMASS(JATOM)
!               JATOM = JATOM + 1
!            END IF
! 110     CONTINUE
! 100  CONTINUE
! f
!      DO I=1,JATOM-1
!         WRITE(LVPRI,*) 'ATOM ',I,' COORD ',(GEOM(I,J),J=1,3),AMASS(I)
!      ENDDO
! f
    do i = 1, nesfp
        amass(i) = rin(i)
        geom(i,1)  = xe(i)
        geom(i,2)  = ye(i)
        geom(i,3)  = ze(i)
    enddo

    angmom(1) = 1.0d0
    angmom(2) = 1.0d0
    angmom(3) = 1.0d0
    call wlkdin(geom,amass,nesfp,angmom,tinert,omegad,eigval,eigvec,.true.,planar,linear)
! f      WRITE(LVPRI,*) 'INERTIA TENSOR EIGENVECTORS AND EIGENVALUES'
    do i = 1, 3
        do j = 1, 3
            eiginv(i,j) = eigvec(j,i)
        enddo
    ! f         WRITE(LVPRI,*) (EIGINV(I,J),J=1,3), EIGVAL(I), OMEGA(I), QUAD
    enddo
!      do i=1,3
!         write(lvpri,*) 'tinert',(tinert(i,j),j=1,3)
!      enddo
    do i = 1,3
        do j= 1,3
            vmat(j,i) = 0.0d0
        enddo
    enddo
    nopax = nrots + nrefl
    norder(1) = 1
    norder(2) = 2
    norder(3) = 3
    if (nopax >= 3) then
        do i = 1,3
        !            write(lvpri,*) 'iax',i,isop(i),iax(isop(i))
            do j = 1,3
                dij = 0.0d0
                if (i == j) dij = 1.0d0
                vmat(j,iax(isop(i))) = dij
            enddo
        enddo
    elseif (nopax >= 1) then
        jax = iax(isop(1))
        scal(1) = abs(eiginv(1,jax))
        scal(2) = abs(eiginv(2,jax))
        scal(3) = abs(eiginv(3,jax))
        nmax = 1
        do j = 2,3
            if (scal(j) > scal(nmax)) nmax = j
        enddo
        nshift = mod(nmax-1,3)
        do i = 0,2
            k = mod(i + nshift,3) + 1
            do j =1,3
                vmat(i+1,j) = eiginv(k,j)
            enddo
        enddo
    elseif (nopax == 0) then
        do i = 1,3
            do j = 1,3
                vmat(i,j) = eiginv(i,j)
            enddo
        enddo
    else
        write(*,*) '9'
        stop
    endif

    end subroutine pcmtns
    
    subroutine updcav(coord)

    ! This subroutine is used to update the cavity in a gradient calculation.
    ! It is called in abacus/abaopt.F
    ! It might be useless within the module as the host program always
    ! calls some other higher-level function to perform the update...

    use pedra_ibtfun

#include <pcm_mxcent.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_pcmdef.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_nuclei.h>
#include <pcm_pcmnuclei.h>
#include <pcm_symmet.h>

    real(8) :: coord(3, *)

    integer :: jatom, i, isym, jcord, mulcnt

    if (icesph == 0 .or. icesph == 2) then
        jatom = 0
        do i = 1, nucind
            mulcnt = istbnu(i)
            do isym = 0, maxrep
                if(ibtand(isym,mulcnt) == 0) then
                    jatom = jatom + 1
                    do jcord = 1,3
                        pcmcord(jcord,jatom) = pt(ibtand(isymax(jcord,1),isym))*coord(jcord,i)
                    end do
                end if
            end do
        end do
    end if

    end subroutine updcav
    
    subroutine ordpcm(lvpri, nts, xtscor, ytscor, ztscor, as, privec, idxpri, work, lwork)

    integer :: lvpri, lwork, nts
    real(8) :: work(*)
    real(8) :: privec(4, nts), xtscor(nts), ytscor(nts), ztscor(nts), as(nts)
    integer :: idxpri(nts)
    
    logical :: lchk, lswtch
    real(8) :: xbak, ybak, zbak, abak
    integer :: ibak, i, j, ii

    lswtch = .false.
    
    do i = 1, nts
        privec(1,i) = xtscor(i)
        privec(2,i) = ytscor(i)
        privec(3,i) = ztscor(i)
        privec(4,i) = as(i)
        idxpri(i) = i
    enddo

    if(lswtch) then
        lswtch = .false.
        do i = 1, nts - 1                                                     
            ii = i + 1
            lchk = chktss(privec(1,i),privec(2,i),privec(3,i),privec(4,i),privec(1,ii),privec(2,ii),privec(3,ii),privec(4,ii))
            lswtch = (lswtch .or. lchk)
            if (lchk) then
                xbak         = privec(1,i)
                ybak         = privec(2,i)
                zbak         = privec(3,i)
                abak         = privec(4,i)
                ibak         = idxpri(i)
                privec(1,i)  = privec(1,ii)
                privec(2,i)  = privec(2,ii)
                privec(3,i)  = privec(3,ii)
                privec(4,i)  = privec(4,ii)
                idxpri(i)    = idxpri(ii)
                privec(1,ii) = xbak
                privec(2,ii) = ybak
                privec(3,ii) = zbak
                privec(4,ii) = abak
                idxpri(ii)   = ibak
            end if
        enddo
    endif
    
    do i = 1, nts
        write(lvpri, '(4f15.9)') (abs(privec(j,i)),j=1,4)
    enddo

    end subroutine ordpcm

    logical function chktss(x1,y1,z1,a1,x2,y2,z2,a2)

    real(8) :: x1, y1, z1, a1, x2, y2, z2, a2

    if(abs(x1) < abs(x2)) then
        chktss = .true.
    else if(abs(x1) > abs(x2)) then
        chktss = .false.
    else if(abs(y1) < abs(y2)) then
        chktss = .true.
    else if(abs(y1) > abs(y2)) then
        chktss = .false.
    else if(abs(z1) < abs(z2)) then
        chktss = .true.
    else if(abs(z1) > abs(z2)) then
        chktss = .false.
    else if(abs(a1) < abs(a2)) then
        chktss = .true.
    else
        chktss = .false.
    endif

    end function chktss
          
    function typlab(i)
    
    character(16)       :: typlab            
    integer, intent(in) :: i
    
    if(i == 4) then
        typlab='ALL EDGE IS FREE'
    else if(i == 1) then
        typlab='2ND VERT COVERED'
    else if(i == 2) then
        typlab='1ST VERT COVERED'
    else if(i == 3) then
        typlab='EDGE  CUT  TWICE'
    else if(i == 0) then
        typlab='ALL EDGE COVERED'
    else
        typlab='UNDEFINED CASE!!'
    end if
    
    end function typlab 

    end module pedra_cavity
