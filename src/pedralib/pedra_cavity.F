C
C  -- dalton/sirius/sircav.F --
C     (Luca Frediani)
C
C/* Deck pedram */
      SUBROUTINE PEDRA_M_(WORK,LWORK)
C
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_iratdef.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_infpri.h>
C
      LOGICAL SOME
      DIMENSION WORK(*)
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
C
C     ----- set memory pointers for polyhedra setup -----
C
      SOME = IPRPCM.NE.-5
C
C     maximum number of tessalations is not known, take worst case
C     maximum number of spheres is not known, take worst case
C
      NUMTS  = MXTS
      NUMSPH = MXSP
      NATM   = MXCENT
      NUMVER = MXVER
C
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
      IF (LAST .GT. LWORK) CALL ERRWRK_('PEDRA_M_',LAST,LWORK)
      LWRK   = LWORK - LAST + 1
      IF(SOME) WRITE(LVPRI,910) LAST
C
      CALL PEDRA_(WORK(LINTSP),WORK(LVERT),WORK(LCENTR),WORK(LNEWSP),
     $     WORK(LICAV1),WORK(LICAV2),WORK(LX),WORK(LY),WORK(LZ),
     $     WORK(LJTR),WORK(LCV),NUMTS,NUMSPH,NUMVER,NATM,SOME,
     $     WORK(LAST),LWRK,WORK(LNPERM))
C
Ckm      CALL QEXIT('PEDRA_')
      RETURN
C
  910 FORMAT(/' MEMORY USED TO GENERATE CAVITY =',I10/)
      END
c/* Deck pedra*/
      SUBROUTINE PEDRA_(INTSPH,VERT,CENTR,NEWSPH,ICAV1,ICAV2, XVAL,YVAL,
     $      ZVAL,JTR,CV,NUMTS,NUMSPH,NUMVER,NATM,SOME,WORK,LWORK,NPERM)
C
#include <pcm_implicit.h>
#include <pcm_maxorb.h>
#include <pcm_maxaqn.h>
#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_codata.h>
#include <pcm_mxcent.h>
#include <pcm_infpri.h>
#include <pcm_symmet.h>
C
      LOGICAL SOME
C
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_nuclei.h>
#include <pcm_pcmnuclei.h>
C
      DIMENSION INTSPH(NUMTS,10),VERT(NUMTS,10,3),CENTR(NUMTS,10,3),
     *          NEWSPH(NUMSPH,2),ICAV1(NATM),ICAV2(NATM),
     *          XVAL(NUMTS),YVAL(NUMTS),ZVAL(NUMTS)
      DIMENSION PP(3),PP1(3),PTS(3,10),CCC(3,10),ROTCAV(3,3)
      DIMENSION WORK(LWORK)
      DIMENSION JTR(NUMTS,*),CV(NUMVER,*)
      DIMENSION NPERM(NUMSPH,*)

      Save D0,First
      DATA D0,FIRST/0.0D0,0.0174533D0/
      DATA ROTCAV/1.0D0, 0.0D0, 0.0D0,
     &            0.0D0, 1.0D0, 0.0D0, 
     &            0.0D0, 0.0D0, 1.0D0/ 
      
C
C     Se stiamo costruendo una nuova cavita' per il calcolo del
C     contributo dispersivo e repulsivo:
C
      IDISREP = 0
      IPRCAV = 0
C     
C
C se icesph=0
C     legge i raggi dall'input e fa coincidere i centri con gli atomi
C se icesph=1
C     legge centri e raggi dall'input
C se icesph=2
C     legge i raggi dall'input e fa coincidere i centri con
C     alcuni atomi definiti dall'indice ina(k) con k=1,NESFP
C     es: xe(k)=c(1,ina(k))
C
Clf      IF(ICESPH.LE.0) THEN
Clf         DO J=1,NESFP
Clf
Clfc           XE(J)=CORD(1,J)
Clfc           YE(J)=CORD(2,J)
Clfc           ZE(J)=CORD(3,J)
ClfClf            write(lvpri,*) 'icesph',icesph,nesfp
Clf            XE(J)=PCMCORD(1,J)
Clf            YE(J)=PCMCORD(2,J)
Clf            ZE(J)=PCMCORD(3,J)
Clf            INA(J)=J
Clf         ENDDO
ClfC
Clf      END IF
Clf      IF(ICESPH.EQ.2)THEN
Clf         DO J=1,NESFP
Clf           XE(J)=PCMCORD(1,INA(J))
Clf           YE(J)=PCMCORD(2,INA(J))
Clf           ZE(J)=PCMCORD(3,INA(J))
Clf         ENDDO
Clf      END IF
Clf      NESF=NESFP
C      print * ,'at the beginning of pedra', nesf
C      write(lvpri,*) 'at the beginning of pedra', nesf
C      do i=1,nesf
C         write(lvpri,1111) xe(i), ye(i), ze(i), re(i)
C      enddo 
C
C  PEDRA_ prevede che i dati geometrici siano espressi in ANGSTROM :
C  vengono trasformati, e solo alla fine i risultati tornano in bohr.
C
   90 CONTINUE
      DO I=1,NESFP
         XE(I)=XE(I)
         YE(I)=YE(I)
         ZE(I)=ZE(I)
         RE(I) = RIN(I) * ALPHA(I)
         write (lvpri,*) xe(i), ye(i), ze(i), re(i)
      ENDDO
C
      CALL DZERO_(VERT,NUMTS*10*3)
      CALL DZERO_(CENTR,NUMTS*10*3)
C
C                   creation of new spheres
C
      DO N = 1, NESF
         NEWSPH(N,1) = 0
         NEWSPH(N,2) = 0
      ENDDO
C
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
C
C     check on the number of spheres
C
      WRITE(LVPRI,*) 'Number of extra spheres = ', NE
      IF (NE .GT. MXSP) THEN
        WRITE(LVPRI,*) NE
        WRITE(LVPRI,*) 'Too many spheres (Greater than MXSP=200).'
        STOP
      ENDIF
C

      DO 120 I=NN,NE
         NES=I-1
         DO 130 J=1,NES
            RIJ2=(XE(I)-XE(J))**2+
     $           (YE(I)-YE(J))**2+
     $           (ZE(I)-ZE(J))**2
            RIJ=SQRT(RIJ2)
            RJD=RE(J)+RSOLV
            TEST1=RE(I)+RJD+RSOLV
            IF(RIJ.GE.TEST1) GO TO 130
            REG=DMAX1(RE(I),RE(J))
            REP=DMIN1(RE(I),RE(J))
            REG2=REG*REG
            REP2=REP*REP
            TEST2=REP*SENOM+sqrt(REG2-REP2*COSOM2)
            IF(RIJ.LE.TEST2) GO TO 130
            REGD2=(REG+RSOLV)*(REG+RSOLV)
            TEST3=(REGD2+REG2-RTDD2)/REG
            IF(RIJ.GE.TEST3) GO TO 130
            DO 140 K=1,NEV
               IF(K.EQ.J .OR. K.EQ.I) GO TO 140
               RJK2=(XE(J)-XE(K))**2+
     $              (YE(J)-YE(K))**2+
     $              (ZE(J)-ZE(K))**2
               IF(RJK2.GE.RIJ2) GO TO 140
               RIK2=(XE(I)-XE(K))**2+
     $              (YE(I)-YE(K))**2+
     $              (ZE(I)-ZE(K))**2
               IF(RIK2.GE.RIJ2) GO TO 140
               RJK=sqrt(RJK2)
               RIK=sqrt(RIK2)
               SP=(RIJ+RJK+RIK)/2.0D0
               HH=4*(SP*(SP-RIJ)*(SP-RIK)*(SP-RJK))/RIJ2
               REO=RE(K)*FRO
               IF(K.GE.NE)REO=0.0002D0
               REO2=REO*REO
               IF(HH.LT.REO2) GO TO 130
  140       CONTINUE
            REPD2=(REP+RSOLV)**2
            TEST8=SQRT(REPD2-RTDD2)+sqrt(REGD2-RTDD2)
            IF(RIJ.LE.TEST8)GO TO 150
            REND2=REGD2+REG2-(REG/RIJ)*(REGD2+RIJ2-REPD2)
            IF(REND2.LE.RTDD2) GO TO 130
            REN=sqrt(REND2)-RSOLV
            FC=REG/(RIJ-REG)
            TEST7=REG-RE(I)
            KG=I
            KP=J
            IF(TEST7.LE.0.000000001D0) GO TO 160
            KG=J
            KP=I
  160       FC1=FC+1.0
            XEN=(XE(KG)+FC*XE(KP))/FC1
            YEN=(YE(KG)+FC*YE(KP))/FC1
            ZEN=(ZE(KG)+FC*ZE(KP))/FC1
            ITYPC = 1
            GO TO 170
  150       R2GN=RIJ-REP+REG
            RGN=R2GN/2.0D0
            FC=R2GN/(RIJ+REP-REG)
            FC1=FC+1.0D0
            TEST7=REG-RE(I)
            KG=I
            KP=J
            IF(TEST7.LE.0.000000001D0) GO TO 180
            KG=J
            KP=I
  180       XEN=(XE(KG)+FC*XE(KP))/FC1
            YEN=(YE(KG)+FC*YE(KP))/FC1
            ZEN=(ZE(KG)+FC*ZE(KP))/FC1
            REN=sqrt(REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ))-RSOLV
  170       NET=NET+1
            XE(NET)=XEN
            YE(NET)=YEN
            ZE(NET)=ZEN
            RE(NET)=REN
C
C     Nella matrice NEWSPH(NESF,2) sono memorizzati i numeri delle
C     sfere "generatrici" della nuova sfera NET: se la nuova sfera e'
C     del tipo A o B entrambi i numeri sono positivi, se e' di tipo
C     C il numero della sfera "principale" e' negativo
C     (per la definizione del tipo si veda JCC 11, 1047 (1990))
C
            IF(ITYPC.EQ.0) THEN
               NEWSPH(NET,1) = KG
               NEWSPH(NET,2) = KP
            ELSEIF(ITYPC.EQ.1) THEN
               NEWSPH(NET,1) = - KG
               NEWSPH(NET,2) = KP
            ENDIF
C
  130    CONTINUE
      NEV=NET
  120 CONTINUE
      IF(NET.NE.NE) GO TO 110
      NESF=NET
C
c     costruzione della tabella di permutazione delle sfere
C
Clf      WRITE(LVPRI,*)'MAXREP', MAXREP
C
      CALL SPHPER_(NESF,NESFP,NUMSPH,NPERM,NEWSPH,XE,YE,ZE,RE)
C
Clf Determination of eigenvalues and eigenvectors of the inertia tensor
C
C      KATOM = NATOMS + NFLOAT

C     Modification by Ville & Chris for new cavity gen
      KATOM = NESFP
      LGEOM = 1
      LMASS = LGEOM + 3 * KATOM
      LEND  = LMASS + KATOM
      LWRK =  LWORK - LEND
Clf      write(lvpri,*) 'katom',katom,lgeom,lmass,lend
      IF (LWRK .LT. 0) THEN
      WRITE(*,*) '2'
      STOP
      ENDIF

      IF (KATOM .GT. 1) THEN
         CALL PCMTNS_(ROTCAV,WORK(LGEOM),WORK(LMASS),KATOM)
      END IF
C     
C    Division of the surface into tesserae
C
      VOL=D0
      STOT=D0
C
C     Controlla se ciascuna tessera e' scoperta o va tagliata
C
      NN = 0
      DO 300 NSFE = 1, NESF
         XEN = XE(NSFE)
         YEN = YE(NSFE)
         ZEN = ZE(NSFE)
         REN = RE(NSFE)
         IF(IPRCAV.GE.10) THEN
            WRITE(LVPRI,1242) NSFE
            WRITE(LVPRI,1243) XEN,YEN,ZEN,REN
         END IF
C     Default options (same as traditional GEPOL)
         IPtype = 2
         IPFlag = 0
         ITSNUM = 60
Clf Which type of tessellation?
C         IF(IPOLYG.GT.0) THEN
C            IPFlag = 0
C            ITSNUM = IPolyg
C         ELSEIF(IPolyg.lt.0) THEN
C            IPFlag = 1
C            TsAre = -1.0D-03*IPolyg
Clf         ENDIF
         Call POLYGEN_(IPFLAG,AREATS,ITSNUM,XEN,YEN,ZEN,REN,ITSEFF,
     $        CV,JTR,NPerm,NSFE,NUMTS,NUMSPH,NUMVER,ROTCAV)
         IF(IPRCAV.GE.10) WRITE(LVPRI,*)'AFTER POLYGEN_. ITSEFF=',ITSEFF
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
            IF(IPRCAV.GE.10) THEN
               WRITE(LVPRI,1244) N1,N2,N3
               WRITE(LVPRI,*) 'VERTICES'' COORDINATES'
               CALL OUTPUT_(PTS,1,3,1,3,3,3,1,LVPRI)
            END IF
            DO JJ = 1, 3
               PP(JJ) = D0
               PP1(JJ) = D0
            ENDDO
C
C     Per ciascuna tessera, trova la porzione scoperta e ne
C     calcola l'area con il teorema di Gauss-Bonnet; il punto rappresentativo
C     e' definito come media dei vertici della porzione scoperta di tessera
C     e passato in PP (mentre in PP1 ci sono le coordinate del punto sulla
C     normale interna).
C     I vertici di ciascuna tessera sono conservati in VERT(3,10,MxTs), il
C     numero di vertici di ciascuna tessera e' in NVERT(MxTs), e i centri
C     dei cerchi di ciascun lato sono in CENTR(3,10,MxTs).
C     In INTSPH(MxTs,10) sono registrate le sfere a cui appartengono i lati
C     delle tessere.
C
            IF (IPRCAV.GT.15) THEN
               write(lvpri,9100) NN + 1
               DO IVER = 1, 3
                  WRITE(LVPRI,9110) IVER, (PTS(JCOR,IVER),JCOR=1,3)
               END DO
            END IF
            CALL TESSERA_(NSFE,NV,PTS,CCC,PP,PP1,AREA,INTSPH,NUMTS)
            IF(IPRCAV.GE.10) THEN
               WRITE(LVPRI,*) 'AFTER TESSERA_ ROUTINE'
               WRITE(LVPRI,1245) NV
               WRITE(LVPRI,*) 'VERTICES'' COORDINATES'
               CALL OUTPUT_(PTS,1,3,1,NV,3,NV,1,LVPRI)
            END IF
            IF(AREA.EQ.D0) THEN 
               WRITE(LVPRI,*) 'ZERO AREA IN TESSERA_', NN + 1
               GOTO 310
            END IF
            IF (NN.GE.MXTS) THEN
               WRITE(LVPRI,*) 'TOO MANY TESSERA_E IN PEDRA_'
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
            IF(IPRCAV.GE.10) THEN
               WRITE(LVPRI,1246) NN
               WRITE(LVPRI,1247) XTSCOR(NN),YTSCOR(NN),ZTSCOR(NN)
               WRITE(LVPRI,1248) XVAL(NN),YVAL(NN),ZVAL(NN)
               WRITE(LVPRI,1249) AS(NN),ISPHE(NN),NVERT(NN)
            END IF
 310        CONTINUE
Clf
 9100    format('Before Tessera n.',I4)
 9110    format('XYZ vert n.',I4,3f15.9)
 9120    format('After Tessera n.',I4)
 9130    format('XYZ centr n.',I4,3f15.9)
 9140    format('XYZA',4f15.9)
Clf
 300  CONTINUE
      NTS = NN
C
C     Verifica se due tessere sono troppo vicine
      TEST = 0.02D0
      TEST2 = TEST*TEST
      DO 400 I = 1, NTS-1
         IF(AS(I).EQ.D0) GOTO 400
         XI = XTSCOR(I)
         YI = YTSCOR(I)
         ZI = ZTSCOR(I)
         II = I + 1
         DO 410 J = II , NTS
            IF(ISPHE(I).EQ.ISPHE(J)) GOTO 410
            IF(AS(J).EQ.D0) GOTO 410
            XJ = XTSCOR(J)
            YJ = YTSCOR(J)
            ZJ = ZTSCOR(J)
            RIJ = (XI-XJ)**2 + (YI-YJ)**2 + (ZI-ZJ)**2
            IF(RIJ.GT.TEST2) GOTO 410
            WRITE(LVPRI,9010) I,J,RIJ,TEST2
Cda vedere  IF(AS(I).LT.AS(J)) AS(I) = D0
Cda vedere  IF(AS(I).GE.AS(J)) AS(J) = D0
 410     CONTINUE
 400  CONTINUE
C
Clf Replication of geometrical parameters according to simmetry
C
      CALL REPCAV_(VERT,CENTR,NPERM,NUMTS,NUMSPH)
C
C Prepare data for geomview
C
Clf      write(lvpri,*) 'vertici delle tessere'
C      do i=1,nts
C         do j=1,nvert(i)
C            write(lvpri,1241) i,j,(vert(i,j,k),k=1,3)
C         enddo
C      enddo
C      write(lvpri,*) 'centers of tesserae'
C      do i=1,nts
C         write(lvpri,1240) i, xtscor(i),ytscor(i),ztscor(i),as(i)
C      enddo
C      call flshfo(lvpri)

Clf      IF (IPRPCM .GT. 10) THEN
      IF (.true.) THEN
         KORD  = 1
         KIDX  = KORD + 4 * NTS
         KLAST = KIDX + NTS
         IF (KLAST.LT.LWORK) THEN
            CALL DZERO_(WORK,KLAST)
            call ORDPCM_(lvpri,nts,xtscor,ytscor,ztscor,as,work(kord),
     $           work(kidx),work(klast),lwork)
         ELSE
            WRITE(LVPRI,*)
     $           'WARNING: NOT ENOUGH MEM IN PEDRA_ TO PRINT THE CAVITY'
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
 1243 FORMAT('Coordinates: X=',F12.8,'  Y=',F12.8,'  Z=',F12.8,
     $     '  R=',F12.8)
 1244 FORMAT('VERTICES INDICES: ',3I6)
 1245 FORMAT('NO. OF VERTICES: ',I6)
 1246 FORMAT('ADDING TESSERA_ DATA TO THE LIST. NN=',I6)
 1247 FORMAT('CENTER: ',3F15.11)
 1248 FORMAT('BOH???: ',3F15.11)
 1249 FORMAT('AREA=',F10.8,'ISPHE=',I4,'NVERT=',I3)
      CALL PLOTCAV_(VERT,NUMTS)
C
C***********************************************************
C     Calcola il volume della cavita' con la formula (t. di Gauss):
C                V=SOMMAsulleTESSERE{A r*n}/3
C     dove r e' la distanza del punto rappresentativo dall'origine,
C     n e' il versore normale alla tessera, A l'area della tessera,
C     e * indica il prodotto scalare.
C***********************************************************
      VOL = D0
      DO ITS = 1, NTS
         NSFE = ISPHE(ITS)
C     Trova il versore normale
         XN = (XTSCOR(ITS) - XE(NSFE)) / RE(NSFE)
         YN = (YTSCOR(ITS) - YE(NSFE)) / RE(NSFE)
         ZN = (ZTSCOR(ITS) - ZE(NSFE)) / RE(NSFE)
C     Trova il prodotto scalare
         PROD = XTSCOR(ITS)*XN + YTSCOR(ITS)*YN + ZTSCOR(ITS)*ZN
         VOL = VOL + AS(ITS) * PROD / 3.D0
      ENDDO
C***********************************************************
C     Stampa la geometria della cavita'
C***********************************************************
      STOT=D0
      CALL DZERO_(SSFE,NESF)
      DO I=1,NTS
         K=ISPHE(I)
         SSFE(K)=SSFE(K)+AS(I)
      ENDDO
C
      IF(SOME) THEN
         WRITE(LVPRI,9020) NESF
          DO 520 I=1,NESF
            WRITE(LVPRI,9030) I,XE(I),YE(I),ZE(I),RE(I),SSFE(I)
            STOT=STOT+SSFE(I)
  520     CONTINUE
         WRITE(LVPRI,9040) NTS,STOT,VOL
      END IF
C
C     ----- set up for possible gradient calculation -----
C           DONE ONLY IF WE HAVE SPHERES ON ATOMS
C
      IF(ICESPH.NE.1) THEN
C
        IF(NESFP.GT.NUCDEP) THEN
            WRITE(LVPRI,9050) NESFP,NUCDEP
            WRITE(*,*) '3'
            STOP
         END IF
C
         CALL DZERO_(DERCEN,MXSP*MXCENT*3*3)
         CALL DZERO_(DERRAD,MXSP*MXCENT*3)
         DO NSFE = 1,NUCDEP
            NATSPH=0
            NSFER=NSFE
            IF(ICESPH.EQ.2)THEN
               NATSPH=1
               DO JJ=1,NESFP
                  IF(INA(JJ).EQ.NSFE)THEN
                     NATSPH=0
                     NSFER=JJ
                  ENDIF
               ENDDO
            ENDIF
            IF(NATSPH.EQ.0) then
               DO ICOORD = 1, 3
                  CALL CAVDER_(NSFE,NSFER,ICOORD,INTSPH,NEWSPH)
               ENDDO
            ENDIF
         ENDDO
      END IF
      IF(IPRSOL .GT. 5) THEN
         WRITE(LVPRI,9070)
         WRITE(LVPRI,9090) (I,ISPHE(I),AS(I),XTSCOR(I),YTSCOR(I),
     *                   ZTSCOR(I),XTSCOR(NTS+I),YTSCOR(NTS+I),
     *                   ZTSCOR(NTS+I), I=1,NTS)
      END IF
C
C     Call the routine which checks if the cavity is single or divided.
C
      CALL CAVSPL_(ICAV1,ICAV2,NCAV1,NCAV2,NPCMN,SOME)
C
C     The dispersion calculation is allowed only in the case of
C     single cavity.
C
      IF(NCAV2.NE.0) IDISP=0
C
      IDISREP = 0
      IRETCAV = 0
      IF (IDISP.EQ.2) IDISP = 1
      RETURN
C
 9000 FORMAT(10X,'-- CENTER OF CHARGE --'/
     *        ' X =',F8.4,' A  Y =',F8.4,' A  Z =',F8.4,' A')
 9010 FORMAT(/' WARNING: The distance between center of tessera',I4,
     *         'and',I4,' is ',F8.6,', less than ',F8.6,' AU'/)
 9020 FORMAT(/' Total number of spheres =',I5/
     *        ' Sphere             Center  (X,Y,Z) (AU)            ',
     *           '   Radius (AU)      Area (AU^2)')
 9030 FORMAT(I4,4F15.9,F15.9)
 9040 FORMAT(/' Total number of tesserae =',I8
     *       /' Surface area =',F14.8,' (AU^2)    Cavity volume =',
     *            F14.8,' (A^3)')
 9050 FORMAT( ' PEDRA: CONFUSION ABOUT SPHERE COUNTS. NESFP, NATM=',2I6)
 9060 FORMAT(/' ADDITIONAL MEMORY NEEDED TO SETUP GRADIENT RUN=',I10)
 9061 FORMAT(/' ADDITIONAL MEMORY NEEDED TO SETUP IEF RUN=',I10)
 9070 FORMAT(/' ***  PARTITION OF THE SURFACE  ***'
     &      //' TESSERA_  SPHERE   AREA   X Y Z TESSERA CENTER  ',
     &        'X Y Z NORMAL VECTOR')
 9090 FORMAT(2I4,7F12.7)
      END
C/* Deck SphPer */
      SUBROUTINE SPHPER_(NESF,NESF0,NUMSPH,NPERM,NEWSPH,XE,YE,ZE,RE)
C
#include <pcm_implicit.h>
#include <pcm_maxorb.h>
#include <pcm_maxaqn.h>
#include <pcm_mxcent.h> 
#include <pcm_priunit.h>
#include <pcm_symmet.h>
C
Clf#include <pcm_trans.h>
C
      DIMENSION XE(*),YE(*),ZE(*),RE(*),V1(3),V2(3)
      DIMENSION NPERM(NUMSPH,*),NEWSPH(NUMSPH,2)
C
#include <pcm_ibtfun.h>
C
c  crea la tabella di permutazione delle sfere
c  e completa per simmetria le sfere aggiunte

      
      NESF1=NESF
C
      DO I=1,NESF
        V1(1)=XE(I)
        V1(2)=YE(I)
        V1(3)=ZE(I)
        R1=RE(I)
C
c  ciclo sulle operazioni di simmetria
C
        
        DO J=1,MAXREP
           DO K=1,NESF
              DO L=1,3
                 V2(L)=PT(IBTAND(ISYMAX(L,1),J))*V1(L)
              ENDDO
Clf
C            DIFF1=
C     &        DSQRT((XE(K)-V2(1))**2+(YE(K)-V2(2))**2+(ZE(K)-V2(3))**2)
C            DIFF2=DABS(R1-RE(K))
C            IF ((DIFF1.LT.1.0d-3).AND.(DIFF2.LT.1.0d-3)) THEN
Clf
              DIFF1 = DSQRT((XE(K)-V2(1))**2 + (YE(K)-V2(2))**2 +
     $             (ZE(K)-V2(3))**2 + (RE(K)-R1)**2)
              IF ((DIFF1.LT.1.0d-3)) THEN
                 NPERM(I,J+1)=K
                 GOTO 10
              ENDIF 
           ENDDO 
C
C  caso in cui occorre completare la creazione di sfere aggiunte
C  oppure le sfere iniziali sono messe male
C
          If (I.LT.NESF0) THEN
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
            IF (NewSph(i,1).lt.0) THEN
              N1=-N1
              N2=-N2
            ENDIF 
            NewSph(NESF1,1)=N1
            NewSph(NESF1,2)=N2
            WRITE(LVPRI,*) 'new!',NESF1,N1,N2
          ENDIF
 10       CONTINUE
        ENDDO
      ENDDO
c
c  eventuale completamento della tabella delle permutazioni per le
c  nuove sfere
c  
      DO I=NESF+1,NESF1
        V1(1)=XE(I)
        V1(2)=YE(I)
        V1(3)=ZE(I)
        R1=RE(I)
        
c  ciclo sulle operazioni di simmetria
        
        DO J=1,MAXREP
           DO K=1,NESF1
              DO L=1,3
                 V2(L)=PT(IBTAND(ISYMAX(L,1),J))*V1(L)
              ENDDO
              DIFF1 = DSQRT((XE(K)-V2(1))**2 + (YE(K)-V2(2))**2 +
     $             (ZE(K)-V2(3))**2 + (RE(K)-R1)**2)
              IF ((DIFF1.LT.1.0d-3)) THEN
                 NPERM(I,J+1)=K
                 GOTO 20
              ENDIF
           ENDDO
           WRITE(*,*) '4'
           STOP

 20        CONTINUE
        ENDDO
      ENDDO
      
      NESF=NESF1
c
c  vengono segnate le sfere da tassellare: la prima in ordine di ogni
c  classe di sfere equivalenti per simmetria
c
      DO I=1,NESF
         NPERM(I,1)=1
         DO K=1,MAXREP
            IF (NPERM(I,K+1).LT.I) THEN
               NPERM(I,1)=0
               GOTO 30
            ENDIF
         ENDDO
 30      CONTINUE
      ENDDO
 1000 format(8I3)
 2000 format (I4,4F18.8)
      RETURN
      END
C/* Deck Polygen */
      SUBROUTINE POLYGEN_(IPFLAG,TSARE,ITSNUM,XEN,YEN,ZEN,REN,
     $  ITSEFF,CV,JTR,NPERM,NSFE,NUMTS,NUMSPH,NUMVER,ROTCAV)
C
#include <pcm_implicit.h>
      integer ITRVO(60,3),ITREO(60,3),IEDO(90,2)
      integer oldtr(100,100),ednew(90,100),trnew(60,100,100)
#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_pi.h>
#include <pcm_mxcent.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_symmet.h>
C
C Polygen: a program to generate spherical polyhedra with triangular
C faces. An equilater division algorithm is used.
C
      DIMENSION NPERM(NUMSPH,*),CV(NUMVER,*),JTR(NUMTS,*)
      DIMENSION V1(3),V2(3),V3(3),ROTCAV(3,3)
      DATA ITRVO /180*0/
      DATA ITREO /180*0/
      DATA IEDO  /180*0/
      DATA D0/0.0d0/,D1/1.0d0/
C
c  Se la sfera non deve essere tassellata ritorna a pedra con 0 tessere
c
#include <pcm_ibtfun.h>
Clf      IF (NPerm(NSFE,1).EQ.0.or.nsfe.ne.2) THEN      
      IF (NPerm(NSFE,1).EQ.0) THEN
         ITsEff=0
Clf         print *, 'return!', ITsEff,NPerm(NSFE,1)
         RETURN
      ENDIF
c
c  A seconda delle due opzioni di funzionamento (Area ottimale o
c  numero di tessere ottimale ) vengono stabilite il tipo di poliedro
c  e la  frequenza di divisione (NF)
c
c  IPflag = 0 Numero di tessere: si genera il poliedro disponibile
c           con il numero di tessere piu' prossimo a ITSNUM
c  IPflag = 1  Area: si genera il poliedro con il numero di tessere tali
c           da avere un area il piu' simile possibile
C
      IPFLAG = 1
      IF (IPflag.EQ.1) THEN
         ITSNUM = INT (4.0d0 * pi * REN**2 / TsAre + 0.5D0 )
      ENDIF
C
      IF (ITSNUM.gt.MxTs) THEN
         WRITE(*,*) '5', ITSNUM, MxTs
         STOP
         ENDIF

C  costruisce la tassellazione iniziale che dipende dal gruppo di
c  simmetria e (per C1) dal numero di tessere richieste
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
C     
      CV(1,1)=D1
      CV(1,2)=D0
      CV(1,3)=D0
C     
      CV(2,1)=D0
      CV(2,2)=D1
      CV(2,3)=D0
C
      CV(3,1)=D0
      CV(3,2)=D0
      CV(3,3)=D1
c
c  determination of NF
c
      
      NF = INT(DSQRT(D1 * ITSNUM / ( D1 * NT0  * 8 )) + 0.5d0)
      IF (NF.le.1) NF=2
      IF (NF.LE.2) WRITE(LVPRI,1000)
      IF (NF.GE.8) WRITE(LVPRI,1100)
      DNF = DFLOAT(NF)
check WRITE(LVPRI,*) 'dopo polydata',NT0,NE0,NV,NF
C
      ITsEff=NT0*NF**2
C
c  -nuovi vertici posti lungo i vecchi spigoli
c  -le regole di calcolo derivano da  calcoli di algebra
c   lineare e trigonometria
c  -i vertici vengono memorizzati in ordine progressivo ed riferiti
c  tramite ednew al vertice di appartenenza
C
      NVPT=NV+1
      DO J = 1,NE0
        DO K = 1,3
          v1(k)=CV(IEDO(j,1),k)
          v2(k)=CV(IEDO(j,2),k)
        ENDDO
        costheta=
     &    ( v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)) /
     &    (sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2) *
     &     sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2))
        theta=acos(costheta)
        sintheta=sin(theta)
        DO l=1,NF-1
          DL = DFLOAT(L)
CLF          m=NF-l
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
      
C
c  -nuovi vertici non posti lungo i vecchi spigoli
c  -a partire dai vertici in ednew secondo regole analoghe alle
c  precedenti, vengono memorizzati a seconda del triangolo in trnew
c  trnew(triangolo,fila,n ordine)=nvertice
C
Clf      DO J=1,NT0
      j=1
      ii=1
      jj=3
      DO L =3,NF
         DL = DFLOAT(L)
         DO N=1,l-2
            DM = DFLOAT(NF - L + 1)
            DL = DFLOAT(L - 1 - N)
            DN = DFLOAT(N)
C            CV(NVPT,1) = DSQRT(DM/DNF)
C            CV(NVPT,2) = DSQRT(DL/DNF)
C            CV(NVPT,3) = DSQRT(DN/DNF)
            CV(NVPT,1) = DM
            CV(NVPT,2) = DL
            CV(NVPT,3) = DN
            DNORM = DSQRT(CV(NVPT,1)**2 + CV(NVPT,2)**2 + CV(NVPT,3)**2)
            CV(NVPT,1) = CV(NVPT,1) / DNORM 
            CV(NVPT,2) = CV(NVPT,2) / DNORM 
            CV(NVPT,3) = CV(NVPT,3) / DNORM 
            trnew(j,l,n+1)=NVPT
            NVPT=NVPT+1
         ENDDO
       ENDDO
Clf       ENDDO
      NV=NVPT-1
C
c  -ora per ogni triangolo originario vengono posti nella matrice oldtr
c  i numeri di ordine dei vertici originali,creati lungo gli spigoli e
c  no dei vecchi triangoli secondo lo schema:
c
c          11                    La disposizione e' secondo
c          | \                   la posizione geometrica del
c          |  \                  vertice.
c          21--22
c          | \  |\               i vertici (1,1)-(NF+1,1)-(NF+1,NF+1) sono
c          |  \ | \              quelli originari.
c          31--32--33
c          | \  |\ | \           i nuovi triangoli sono:
c          |  \ | \|  \          (i,j)-(i+1,j)-(i+1,j+1)
c          41--42--43--44        i=1,...,NF j=1,...,i
c
c                                (i,j)--(i,j+1),(i+1,j+1)
c                                i=2,...,NF j=1,...,i-1
c
      NTPT=1
      do 2310 n=1,NT0
C
c  -1 vecchi spigoli
C
        oldtr(1,1)=ITRVO(n,1)
        oldtr(NF+1,1)=ITRVO(n,2)
        oldtr(NF+1,NF+1)=ITRVO(n,3)
C
c  -2 nuovi vertici lungo i vecchi spigoli
C
        DO l=2,NF
          oldtr(l,1)=ednew(ITREO(n,1),l)
          oldtr(NF+1,l)=ednew(ITREO(n,2),l)
          oldtr(l,l)=ednew(ITREO(n,3),l)
        ENDDO
C
c  -3 nuovi vertici non lungo i vecchi spigoli
C
        DO l=3,NF
          DO m=2,l-1
            oldtr(l,m)=trnew(n,l,m)
          ENDDO
        ENDDO
C
c  -ora si creano i nuovi triangoli
C
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
 2310 CONTINUE
      NV=NVPT-1
      NT=NTPT-1
C
Clf Replication generate the right irreducible part of the cavity
C   for the selected group
C
      CALL PREREP_(NV,NT,ITSEFF,CV,JTR,NUMVER,NUMTS)
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
Clf         WRITE(LVPRI,*)'CV',I
Clf         WRITE(LVPRI,'(3F10.4)')  (CV(i,II),II=1,3)
      ENDDO

c
c  replicazione con operatori locali  
c        nei casi in cui l'o.l. trasforma la sfera in un altra sfera  

      NOpPt=0
      DO ISYMOP = 1, MAXREP
         NTRA=NPerm(NSFE,ISYMOP + 1)
C
C     the If statment is entered only if the sphere which is equivalent
C     to the one we are tesselating now under the operation ISYMOP, is
C     not itself.
C
         IF (NTRA .NE. NSFE) THEN
C
C     If the replication of this part of the cavity has already been
C     done we start checking a new simmetry operation
C
            DO JSYMOP = 1 , ISYMOP-1
               IF(Nperm(NSFE, JSYMOP+1) .EQ. NTRA) Goto 100
            ENDDO
            NOpPt=NOpPt+1
c     
c     riproduzione vertici 
c     
            Do I=1,NV
               II=I+NOpPt*NV
               Do K=1,3
                  CV(II,K) = PT(IBTAND(ISYMAX(K,1),ISYMOP)) * CV(I,K)
               ENDDO
            ENDDO
c
c  riproduzione topologia 
c
            Do i=1,NT
               ii=i+NOpPt*NT
               jj=NOpPt*NV
               Do k=1,3
                  JTR(ii,k)=JTR(i,k)+jj
               ENDDO
            ENDDO
 100        CONTINUE 
         ENDIF
      ENDDO
c
c  aggiornamento indici
c
      NV=NV*(NOpPt+1)
      ITsEff=ITsEff*(NOpPt+1)
C
c scrittura della CV
C
CcheckWRITE(LVPRI,*) 'OFF'
CcheckWRITE(LVPRI,*) NV,ITsEff,NV
      DO i=1,NV
         CV(I,1) = CV(I,1) * REN + XEN
         CV(I,2) = CV(I,2) * REN + YEN
         CV(I,3) = CV(I,3) * REN + ZEN
Clf         WRITE(LVPRI,*)'CV',I
Clf         WRITE(LVPRI,'(3F10.4)')  (CV(i,II),II=1,3)
      ENDDO
C
      RETURN
 1000 FORMAT(/' ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN',
     $       /' IT IS VALUABLE ALMOST ONLY FOR TESTING')
 1100 FORMAT(/' ** WARNING ** A VERY EXPENSIVE TESSELATION ',
     $   'HAS BEEN CHOSEN',
     $       /' IT WILL PROBABLY PRODUCE A HUGE AMOUNT OF TESSERA_E')
      END
c/* Deck Repcav */
      SUBROUTINE REPCAV_(VERT,CENTR,NPERM,NUMTS,NUMSPH)
C
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
C
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_nuclei.h>
#include <pcm_symmet.h>
Clf#include <pcm_trans.h>
c
c  reproduce the irreducibile part of the cavity
c
      DIMENSION VERT(NUMTS,10,3),CENTR(NUMTS,10,3)
      DIMENSION NPERM(NUMSPH,*)
      DATA D0/0.0D0/
C
#include <pcm_ibtfun.h>
C
c  replication of the cavity
c
CcheckWRITE(LVPRI,*)'IN REPCAV_',NTS
      NTSIRR = NTS
      NTS = NTS * (MAXREP + 1)

      IF (NTS.GT.MXTS) THEN 
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
               VERT(II,K,L)  = 
     $              PT(IBTAND(ISYMAX(L,1),ISYMOP)) * VERT(I,K,L)
               CENTR(II,K,L) = 
     $              PT(IBTAND(ISYMAX(L,1),ISYMOP)) * CENTR(I,K,L)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
Clf Moving the normal points in a safe location where they will not be overwritten 
C
      DO I=1,NTSIRR
        II=I+NTSIRR
        II2=I+(MAXREP+1)*NTSIRR
        XTSCOR(II2)=XTSCOR(II)
        YTSCOR(II2)=YTSCOR(II)
        ZTSCOR(II2)=ZTSCOR(II)
      ENDDO
C
Clf      write(lvpri,*) 'see if the pt are correct',nts
      DO ISYMOP=1,MAXREP
         DO I=1,NTSIRR
            II=I+NTSIRR*ISYMOP
            XTSCOR(II)=PT(IBTAND(ISYMAX(1,1),ISYMOP)) * XTSCOR(I)
            YTSCOR(II)=PT(IBTAND(ISYMAX(2,1),ISYMOP)) * YTSCOR(I)
            ZTSCOR(II)=PT(IBTAND(ISYMAX(3,1),ISYMOP)) * ZTSCOR(I)
C            II=I+NTS*(ISYMOP+MAXREP+1)
C            II2=I+NTS
C            XTSCOR(II)=PT(IBTAND(ISYMAX(1,1),ISYMOP)) * XTSCOR(II2)
C            YTSCOR(II)=PT(IBTAND(ISYMAX(2,1),ISYMOP)) * YTSCOR(II2)
C            ZTSCOR(II)=PT(IBTAND(ISYMAX(3,1),ISYMOP)) * ZTSCOR(II2)
         ENDDO
      ENDDO
Clf
C      do i=1,nts
C         do j=1,nts
C            distance=sqrt((xtscor(j)-xtscor(i))**2
C     $                   +(ytscor(j)-ytscor(i))**2
C     $                   +(ztscor(j)-ztscor(i))**2)
C            if(i.ne.j .and. distance.lt.1.0d-8) then
C               write(lvpri,*) 'tesserae too close:',i,j
C            endif
C         enddo
C      enddo
Clf
      RETURN
      END
C/* Deck tessera */
      SUBROUTINE TESSERA_(NS,NV,PTS,CCC,PP,PP1,AREA,INTSPH,NUMTS)
C
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_infpri.h>
C
      DIMENSION PTS(3,10),CCC(3,10),PP(3),PP1(3),INTSPH(NUMTS,10)
      DIMENSION P1(3),P2(3),P3(3),P4(3),P4BIS(3),POINT(3),P5(3),P6(3),
     *     PSCR(3,10),CCCP(3,10),POINTL(3,10),
     *     IND(10),LTYP(10),INTSCR(10),NTRHSO(10)
      LOGICAL LAN
      DIMENSION DISTCK(10,10)
      CHARACTER*16 TYPLAB_
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
C
C     Coord. del centro che sottende l`arco tra i vertici
C     n e n+1 (per i primi tre vertici e' sicuramente il centro della
C     sfera) e sfera alla cui intersezione con NS appartiene l'arco (se
C     appartiene alla sfera originaria INTSPH(numts,N)=NS)
C
C
      IPRCAV = 0
      LAN = .FALSE.
      AREA = 0.0D+00
      DO J=1, 3
         CCC(1,J) = XE(NS)
         CCC(2,J) = YE(NS)
         CCC(3,J) = ZE(NS)
      ENDDO
      IF(IPRCAV.GE.10) THEN
         CALL AROUND_('INPUT DATA IN TESSERA_')
         WRITE(LVPRI,1000) NS,NV,NUMTS
         WRITE(LVPRI,*) '=======PTS=======' 
         CALL OUTPUT_(PTS,1,3,1,3,3,3,1,LVPRI)
         WRITE(LVPRI,*) '=======CCC=======' 
         CALL OUTPUT_(CCC,1,3,1,3,3,3,1,LVPRI)
      END IF
C     
C     INTSPH viene riferito alla tessera -numts-, e in seguito riceve il
C     numero corretto.
C
      DO N = 1, 3
         INTSPH(NUMTS,N) = NS
      ENDDO
C
C     Loop sulle altre sfere
C
      DO 150 NSFE1=1,NESF
         IF(NSFE1.EQ.NS) GO TO 150
         IF(IPRCAV.GE.10) THEN
            WRITE(LVPRI,1005) NSFE1
            WRITE(LVPRI,1010) XE(NSFE1),YE(NSFE1),ZE(NSFE1),RE(NSFE1)
         END IF
C     
C     Memorizza i vertici e i centri che sottendono gli archi
C
         DO J =1, NV
            INTSCR(J) = INTSPH(NUMTS,J)
            DO I = 1,3
               PSCR(I,J) = PTS(I,J)
               CCCP(I,J) = CCC(I,J)
            ENDDO
         ENDDO
         IF(IPRCAV.GE.10) THEN
            CALL AROUND_('ACTUAL TESSERA_ STATUS')
            WRITE(LVPRI,1000) NS,NV,NUMTS
            WRITE(LVPRI,*) '=======PSCR=======' 
            CALL OUTPUT_(PSCR,1,3,1,NV,3,NV,1,LVPRI)
            WRITE(LVPRI,*) '=======CCCP=======' 
            CALL OUTPUT_(CCCP,1,3,1,NV,3,NV,1,LVPRI)
         END IF
         
C     
         ICOP = 0
         DO J =1, 10
            IND(J) = 0
            LTYP(J) = 0
         ENDDO
C     
C     Loop sui vertici della tessera considerata
C
         DO 100 I=1,NV
            DELR2=(PTS(1,I)-XE(NSFE1))**2+(PTS(2,I)-YE(NSFE1))**2+
     *           (PTS(3,I)-ZE(NSFE1))**2
            DELR=SQRT(DELR2)
            IF (IPRPCM.GE.10) THEN
               WRITE(LVPRI,1015) DELR,RE(NSFE1)
            ENDIF
            DIFFDR = DELR - RE(NSFE1)
            IF(DIFFDR.LT.0.0D0) THEN
               IND(I) = 1
               ICOP = ICOP+1
            END IF
 100     CONTINUE
C     Se la tessera e' completamente coperta, la trascura
         IF(ICOP.EQ.NV) THEN
            IF(IPRCAV.GE.10) WRITE(LVPRI,1020) NSFE1
            RETURN
C           ******
         END IF

C
C     Controlla e classifica i lati della tessera: LTYP = 0 (coperto),
C     1 (tagliato con il II vertice coperto), 2 (tagliato con il I
C     vertice coperto), 3 (bitagliato), 4 (libero)
C     Loop sui lati
C
         DO L = 1, NV
            IV1 = L
            IV2 = L+1
            IF(L.EQ.NV) IV2 = 1
            IF(IND(IV1).EQ.1.AND.IND(IV2).EQ.1) THEN
               LTYP(L) = 0
            ELSE IF(IND(IV1).EQ.0.AND.IND(IV2).EQ.1) THEN
               LTYP(L) = 1
            ELSE IF(IND(IV1).EQ.1.AND.IND(IV2).EQ.0) THEN
               LTYP(L) = 2
            ELSE IF(IND(IV1).EQ.0.AND.IND(IV2).EQ.0) THEN
               LTYP(L) = 4
               RC2 = (CCC(1,L)-PTS(1,L))**2 + (CCC(2,L)-PTS(2,L))**2 +
     *              (CCC(3,L)-PTS(3,L))**2
               RC = SQRT(RC2)
C
C     Su ogni lato si definiscono 11 punti equispaziati, che vengono
C     controllati
C
               TOL = - 1.0D-10
               DO II = 1, 11
                  POINT(1) = PTS(1,IV1) + 
     $                 II * (PTS(1,IV2)-PTS(1,IV1)) / 11
                  POINT(2) = PTS(2,IV1) + 
     $                 II * (PTS(2,IV2)-PTS(2,IV1)) / 11
                  POINT(3) = PTS(3,IV1) + 
     $                 II * (PTS(3,IV2)-PTS(3,IV1)) / 11
                  POINT(1) = POINT(1) - CCC(1,L)
                  POINT(2) = POINT(2) - CCC(2,L)
                  POINT(3) = POINT(3) - CCC(3,L)
                  DNORM = SQRT(POINT(1)**2 + POINT(2)**2 + POINT(3)**2)
                  POINT(1) = POINT(1) * RC / DNORM + CCC(1,L)
                  POINT(2) = POINT(2) * RC / DNORM + CCC(2,L)
                  POINT(3) = POINT(3) * RC / DNORM + CCC(3,L)
                  DIST = SQRT( (POINT(1)-XE(NSFE1))**2 +
     $                         (POINT(2)-YE(NSFE1))**2 + 
     $                         (POINT(3)-ZE(NSFE1))**2 )
                  IF((DIST - RE(NSFE1)) .LT. TOL) THEN
C     IF(DIST.LT.RE(NSFE1)) then
                     LTYP(L) = 3
                     DO JJ = 1, 3
                        POINTL(JJ,L) = POINT(JJ)
                     ENDDO
                     GO TO 160
                  END IF
               ENDDO
            END IF
 160        CONTINUE
            IF(IPRCAV.GE.10) THEN
               WRITE(LVPRI,1025) L,TYPLAB_(LTYP(L))
            END IF
         ENDDO
C
C     Se la tessera e' spezzata in due o piu' tronconi, la trascura
C
         ICUT = 0
         DO L = 1, NV
            IF(LTYP(L).EQ.1.OR.LTYP(L).EQ.2) ICUT = ICUT + 1
            IF(LTYP(L).EQ.3) ICUT = ICUT + 2
         ENDDO
         ICUT = ICUT / 2
         IF(ICUT.GT.1) THEN
            WRITE(LVPRI,*) 'Tessera cut in pieces and removed.'
            RETURN
         END IF
C     
C     Creazione dei nuovi vertici e lati della tessera
C     Loop sui lati
C     
         N = 1
         IF(IPRCAV.GE.10) WRITE(LVPRI,*) 
     $        'NOW CREATING NEW VERTICES AND EDGES IF NEED BE....'
         DO 300 L = 1, NV
C     Se il lato L e' coperto:
            IF(LTYP(L).EQ.0) GO TO 300
            IV1 = L
            IV2 = L+1
            IF(L.EQ.NV) IV2 = 1
            IF (IPRCAV.GE.10) THEN
               WRITE(LVPRI,1030) L,TYPLAB_(LTYP(L)),LTYP(L)
               WRITE(LVPRI,1035) IV1,IV2
            ENDIF
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Se il lato L e' tagliato (con il I vertice scoperto):
            IF(LTYP(L).EQ.1) THEN
Clf store first point in the final set
               DO JJ = 1, 3
                  PTS(JJ,N) = PSCR(JJ,IV1)
                  CCC(JJ,N) = CCCP(JJ,IV1)
               ENDDO
               INTSPH(NUMTS,N) = INTSCR(IV1)
               N = N+1
C     
C     Trova l'intersezione tra i due vertici del lato L
C     
C     P1 = coord. del primo vertice
C     P2 = coord. del secondo vertice
C     P3 = coord. del centro dell`arco sotteso
C     P4 = coord. dell'intersezione
C
               DO JJ = 1, 3
                  P1(JJ) = PSCR(JJ,IV1)
                  P2(JJ) = PSCR(JJ,IV2)
                  P3(JJ) = CCCP(JJ,IV1)
                  P4(JJ) = 0.0D0
               ENDDO
               INTCAS=0
               IF(IPRCAV.GE.10) THEN
                  WRITE(LVPRI,1040) INTCAS,NSFE1
                  WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                  WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                  WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                  WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
               END IF
               CALL INTER_(P1,P2,P3,P4,NSFE1,INTCAS)
               IF(IPRCAV.GE.10) THEN
                  WRITE(LVPRI,1050) 
                  WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                  WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                  WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                  WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
               END IF
               DIST1 = DSQRT((P4(1)-P1(1))**2 + (P4(2)-P1(2))**2 +
     $                       (P4(3)-P1(3))**2)
               DIST2 = DSQRT((P4(1)-P2(1))**2 + (P4(2)-P2(2))**2 +
     $                       (P4(3)-P2(3))**2)
C     Aggiorna i vertici della tessera e il centro dell'arco
               DO JJ = 1,3
                  PTS(JJ,N) = P4(JJ)
               ENDDO
C     
C     Il nuovo arco sara' sotteso tra questo e il prossimo punto
C     di intersezione: il centro che lo sottende
C     sara' il centro del cerchio di intersezione tra la sfera NS
C     e la sfera NSFE1.
C
               DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+
     *              (ZE(NSFE1)-ZE(NS))**2
               CCC(1,N)=XE(NS)+(XE(NSFE1)-XE(NS))*
     *              (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
               CCC(2,N)=YE(NS)+(YE(NSFE1)-YE(NS))*
     *              (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
               CCC(3,N)=ZE(NS)+(ZE(NSFE1)-ZE(NS))*
     *              (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
               INTSPH(NUMTS,N) = NSFE1
               N = N+1
            END IF
C     
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Se il lato L e' tagliato (con il II vertice scoperto):
            IF(LTYP(L).EQ.2) THEN
C     Trova l'intersezione tra i due vertici del lato L
C
C     P1 = coord. del primo vertice
C     P2 = coord. del secondo vertice
C     P3 = coord. del centro dell`arco sotteso
C     P4 = coord. dell'intersezione
C
               DO JJ = 1, 3
                  P1(JJ) = PSCR(JJ,IV1)
                  P2(JJ) = PSCR(JJ,IV2)
                  P3(JJ) = CCCP(JJ,IV1)
                  P4(JJ) = 0.0D0
               ENDDO
               INTCAS=1
               IF(IPRCAV.GE.10) THEN
                  WRITE(LVPRI,1040) INTCAS,NSFE1
                  WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                  WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                  WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                  WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
               END IF
               CALL INTER_(P1,P2,P3,P4,NSFE1,INTCAS)
               IF(IPRCAV.GE.10) THEN
                  WRITE(LVPRI,1050) 
                  WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                  WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                  WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                  WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
               END IF
C     Aggiorna i vertici della tessera e il centro dell'arco
               DO JJ = 1,3
                  PTS(JJ,N) = P4(JJ)
                  CCC(JJ,N) = CCCP(JJ,IV1)
               ENDDO
               INTSPH(NUMTS,N) = INTSCR(IV1)
               N = N+1
            END IF
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Se il lato e' intersecato due volte:
            IF(LTYP(L).EQ.3) THEN
               DO JJ = 1, 3
                  PTS(JJ,N) = PSCR(JJ,IV1)
                  CCC(JJ,N) = CCCP(JJ,IV1)
               ENDDO
               INTSPH(NUMTS,N) = INTSCR(IV1)
               N = N+1
C
C     Trova l'intersezione tra il primo vertice e un punto intermedio
C     coperto
C
C     P1 = coord. del primo vertice
C     P2 = coord. del secondo vertice
C     P3 = coord. del centro dell`arco sotteso
C     P4 = coord. dell'intersezione
C
               DO JJ = 1, 3
                  P1(JJ) = PSCR(JJ,IV1)
                  P2(JJ) = POINTL(JJ,L)
                  P3(JJ) = CCCP(JJ,IV1)
                  P4(JJ) = 0.0D0
               ENDDO
               INTCAS=0
               IF(IPRCAV.GE.10) THEN
                  WRITE(LVPRI,1040) INTCAS,NSFE1
                  WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                  WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                  WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                  WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
               END IF
               CALL INTER_(P1,P2,P3,P4,NSFE1,INTCAS)
               IF(IPRCAV.GE.10) THEN
                  WRITE(LVPRI,1050) 
                  WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                  WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                  WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                  WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
               END IF
C     Aggiorna i vertici della tessera e il centro dell'arco
               DO JJ = 1,3
                  PTS(JJ,N) = P4(JJ)
               ENDDO
C     
C     Il nuovo arco sara' sotteso tra questo e il prossimo punto
C     di intersezione: il centro che lo sottende
C     sara' il centro del cerchio di intersezione tra la sfera NS
C     e la sfera NSFE1.
C     
               DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+
     *              (ZE(NSFE1)-ZE(NS))**2
               CCC(1,N)=XE(NS)+(XE(NSFE1)-XE(NS))*
     *              (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
               CCC(2,N)=YE(NS)+(YE(NSFE1)-YE(NS))*
     *              (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
               CCC(3,N)=ZE(NS)+(ZE(NSFE1)-ZE(NS))*
     *              (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.0D+00*DE2)
               INTSPH(NUMTS,N) = NSFE1
               N = N+1
C     
C     Trova l'intersezione tra un punto intermedio coperto e il
C     secondo vertice
C     
C     P1 = coord. del primo vertice
C     P2 = coord. del secondo vertice
C     P3 = coord. del centro dell`arco sotteso
C     P4 = coord. dell'intersezione
C     
               DO JJ = 1, 3
                  P1(JJ) = POINTL(JJ,L)
                  P2(JJ) = PSCR(JJ,IV2)
                  P3(JJ) = CCCP(JJ,IV1)
                  P4(JJ) = 0.0D0
               ENDDO
               INTCAS=1
               IF(IPRCAV.GE.10) THEN
                  WRITE(LVPRI,1040) INTCAS,NSFE1
                  WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                  WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                  WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                  WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
               END IF
               CALL INTER_(P1,P2,P3,P4,NSFE1,INTCAS)
               IF(IPRCAV.GE.10) THEN
                  WRITE(LVPRI,1050) 
                  WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                  WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                  WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                  WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
               END IF
               DIST1 = DSQRT((P4(1)-P1(1))**2 + (P4(2)-P1(2))**2 +
     $              (P4(3)-P1(3))**2)
               DIST2 = DSQRT((P4(1)-P2(1))**2 + (P4(2)-P2(2))**2 +
     $              (P4(3)-P2(3))**2)
C     Aggiorna il vertice e il centro dell'arco
               DO JJ = 1,3
                  PTS(JJ,N) = P4(JJ)
                  CCC(JJ,N) = CCCP(JJ,IV1)
               ENDDO
               INTSPH(NUMTS,N) = INTSCR(IV1)
               N = N + 1
            END IF
C     
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Se il lato e' scoperto:
            IF(LTYP(L).EQ.4) THEN
               DO JJ = 1, 3
                  PTS(JJ,N) = PSCR(JJ,IV1)
                  CCC(JJ,N) = CCCP(JJ,IV1)
               ENDDO
               INTSPH(NUMTS,N) = INTSCR(IV1)
               N = N+1
            END IF
 300     CONTINUE
         NV = N - 1
         IF(IPRCAV.GE.10) THEN
            CALL AROUND_('AFTER INTER_SECTION TESSERA_ STATUS')
            WRITE(LVPRI,1000) NS,NV,NUMTS
            WRITE(LVPRI,*) '=======PTS=======' 
            CALL OUTPUT_(PTS,1,3,1,NV,3,NV,1,LVPRI)
            WRITE(LVPRI,*) '=======CCC=======' 
            CALL OUTPUT_(CCC,1,3,1,NV,3,NV,1,LVPRI)
            DO IDX=1,NV
               IDX2=IDX+1
               IF (IDX2.GT.NV) IDX2=1
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
C     Controlla che il numero di vertici creati non sia eccessivo
         IF(NV.GT.10) THEN
            WRITE(LVPRI,*)'TOO MANY VERTICES IN TESSERA_: BYE BYE...'
            STOP
         END IF
 150  CONTINUE
      IF(IPRCAV.GE.10) THEN
         CALL AROUND_('FINAL TESSERA_ STATUS')
         WRITE(LVPRI,1000) NS,NV,NUMTS
         WRITE(LVPRI,*) '=======PSCR=======' 
         CALL OUTPUT_(PSCR,1,3,1,NV,3,NV,1,LVPRI)
         WRITE(LVPRI,*) '=======CCCP=======' 
         CALL OUTPUT_(CCCP,1,3,1,NV,3,NV,1,LVPRI)
      END IF
C     
C     Se la tessera non e' stata scartata, a questo punto ne troviamo
C     l'area e il punto rappresentativo
C
C
Clf We label the edges whose lenght is lower than a certain threshold
C
      NVNEGL = 0
      DO I=1,NV
         NTRHSO(I) = 0
         II=I+1
         IF (I .EQ. NV) II = 1
         X1 = PTS(1,I)
         Y1 = PTS(2,I)
         Z1 = PTS(3,I)
         X2 = PTS(1,II)
         Y2 = PTS(2,II)
         Z2 = PTS(3,II)
         DIST = DSQRT((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)
         IF (DIST .LT. 1.0D-6) THEN 
            NTRHSO(I) = 1
            NVNEGL = NVNEGL + 1
clf            NTRHSO(I) = 0
         END IF
      ENDDO
      NVLEFT = NV - NVNEGL
C
C If we end up with a tessera with less than three long enough edges
C we discard it.
C
      IF (NVLEFT .LT. 3) THEN
         WRITE(LVPRI,*) '*******WARNING*******\\',
     $        ' NEGLECTED TESSERA_: TOO FEW ACCEPTABLE EDGES'
         RETURN
C
C Removal of too small edges to avoid numerical problems
C in area calculation (done only if one ore more edges are 
C to be neglected).
C
      ELSE IF (NVNEGL .GT. 0) THEN
         IVNEW = 0
         DO IVOLD = 1,NV
            IF(NTRHSO(IVOLD).EQ.0) THEN
               IVNEW = IVNEW + 1
               INTSPH(NUMTS,IVNEW) = INTSPH(NUMTS,IVOLD)
               DO K =1,3
                  CCC(K,IVNEW) = CCC(K,IVOLD)
                  PTS(K,IVNEW) = PTS(K,IVOLD)
               ENDDO
            END IF
         END DO
         IF (IVNEW .NE. NVLEFT) THEN
            WRITE(LVPRI,*) 'SIRCAV: BADLY MADE ALGORITHM!'
            STOP
         ENDIF
         NV = NVLEFT
      END IF
C
C Finally we calculate area and center of the tessera!
C
      CALL GAUBON_(NV,NS,PTS,CCC,PP,PP1,AREA,INTSPH,NUMTS)
      RETURN
 1000 FORMAT('NS=',I4,' NV=',I4,' NUMTS=', I4)
 1005 FORMAT('CHECKING INTER_SECTIONS WITH SPHERE N.',I4)
 1010 FORMAT('Coordinates: X=',F12.8,'  Y=',F12.8,'  Z=',F12.8,
     $     ' R=',F12.8)
 1015 FORMAT('VERT-OTHER CENTER DIST:',F10.8,' SPHERE RADIUS',F10.8)
 1020 FORMAT('THIS TESSERA_ IS COMPLETELY COVERED BY SPHERE N.',I4)
 1025 FORMAT('EDGE N.',I4,' STATUS: ',A18)
 1030 FORMAT('EDGE N.',I2,' IS IN STATUS ',A16,' CASE N.',I2)
 1035 FORMAT('VERTICES NUMBERS:',I3,' AND',I3)
 1040 FORMAT('CALLING INTER_, CASE N.',I2,' INTERSECTING SPHERE:',I4)
 1045 FORMAT('P',I1,' X = ',F12.9,' Y = ',F12.9,' Z = ',F12.9)
 1050 FORMAT('AFTER INTER_')
 1055 FORMAT('DIST CHECK EDGE N.',I2,' DIFFERENCE:',F12.10)
      END
C/* Deck inter */
      SUBROUTINE INTER_(P1,P2,P3,P4,NS,I)
C
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_mxcent.h>
#include <pcm_pcmdef.h>
#include <pcm_infpri.h>
C
      DIMENSION P1(3),P2(3),P3(3),P4(3)
      LOGICAL LALOW,LBLOW
C
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
      IPRCAV = 0
C
C     Trova il punto P4, sull`arco P1-P2 sotteso dal centro P3, che
C     si trova sulla superficie della sfera NS
C     P4 e' definito come combinazione lineare di P1 e P2, con
C     il parametro ALPHA ottimizzato per tentativi.
C
      R2 = (P1(1)-P3(1))**2+(P1(2)-P3(2))**2+(P1(3)-P3(3))**2
      R = SQRT(R2)
Clf Don't change this threshold unless you know exactly what you are doing.
      TOL = 1.0D-12
      ALPHAT = 0.5D+00
      DELTA = 0.0D+00
Clf Check distances of P1 and P2 from P3
      P1P3 = SQRT((P1(1)-P3(1))**2+(P1(2)-P3(2))**2+(P1(3)-P3(3))**2)
      P2P3 = SQRT((P2(1)-P3(1))**2+(P2(2)-P3(2))**2+(P2(3)-P3(3))**2)
      IF(IPRCAV.GE.10) THEN
         WRITE(LVPRI,1000) P1P3-P2P3
      END IF
Clf Check if one or both starting points are below the threshold
      LALOW=.FALSE.
      LBLOW=.FALSE.
      DIFF2A= (P1(1)-XE(NS))**2 + (P1(2)-YE(NS))**2 + (P1(3)-ZE(NS))**2
      DIFFA = DABS(DSQRT(DIFF2A) - RE(NS))
      DIFF2B= (P2(1)-XE(NS))**2 + (P2(2)-YE(NS))**2 + (P2(3)-ZE(NS))**2
      DIFFB = DABS(DSQRT(DIFF2B) - RE(NS))
      LALOW=(DIFFA.LT.TOL)
      LBLOW=(DIFFB.LT.TOL)
      IF(LALOW.AND..NOT.LBLOW) THEN
         IF(IPRCAV.GE.10) WRITE (LVPRI,*) 
     $        'INTER_: TAKEN FIRST POINT'
         DO J=1,3
            P4(J)=P1(J)
         END DO
         GOTO 100
      END IF
      IF(LBLOW.AND..NOT.LALOW) THEN
         IF(IPRCAV.GE.10) WRITE (LVPRI,*) 
     $        'INTER_: TAKEN SECOND POINT'
         DO J=1,3
            P4(J)=P2(J)
         END DO
         GOTO 100
      END IF
      IF(LALOW.AND.LBLOW) THEN
         IF(DIFFA.LE.DIFFB) THEN
            IF(IPRCAV.GE.10) WRITE (LVPRI,*) 
     $           'INTER_: TAKEN FIRST POINT'
            DO J=1,3
               P4(J)=P1(J)
            END DO
         ELSE
            IF(IPRCAV.GE.10) WRITE (LVPRI,*) 
     $           'INTER_: TAKEN SECOND POINT'
            DO J=1,3
               P4(J)=P2(J)
            END DO
         END IF
         GOTO 100
      END IF
C
C     Start iterations
C
      M = 1
 10   CONTINUE
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
      IF(ABS(DIFF).LT.TOL) GOTO 100
        IF(I.EQ.0) THEN
         IF(DIFF.GT.0.0D+00) DELTA = 1.0D+00/(2.0D+00**DBLE(M+1))
         IF(DIFF.LT.0.0D+00) DELTA = - 1.0D+00/(2.0D+00**DBLE(M+1))
         M = M + 1
        ELSE IF(I.EQ.1) THEN
         IF(DIFF.GT.0.0D+00) DELTA = - 1.0D+00/(2.0D+00**DBLE(M+1))
         IF(DIFF.LT.0.0D+00) DELTA = 1.0D+00/(2.0D+00**DBLE(M+1))
         M = M + 1
        ELSE
         WRITE(*,*) '7'
         STOP
        END IF
        IF (M.GT.300)THEN
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
C
C Final printing and return
C
 100  CONTINUE
      IF(IPRPCM.GE.10) THEN
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
 1015 FORMAT('SPHERE',' X = ',F12.9,' Y = ',F12.9,' Z = ',F12.9,
     $     ' R = ',F12.9)
      END
C/* Deck gaubon */
      SUBROUTINE GAUBON_(NV,NS,PTS,CCC,PP,PP1,AREA,INTSPH,NUMTS)
C
#include <pcm_implicit.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0)
#include <pcm_priunit.h>
#include <pcm_pcmdef.h>      
#include <pcm_mxcent.h>
#include <pcm_pi.h>
#include <pcm_infpri.h>
C
      DIMENSION PTS(3,10),CCC(3,10),PP(3),PP1(3),INTSPH(NUMTS,10),
     $     NTRHSO(10),BETA(10)
      DIMENSION P1(3),P2(3),P3(3),U1(3),U2(3),PHIN(10),WEIGHT(0:10)
C
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
C
C     Sfrutta il teorema di Gauss-Bonnet per calcolare l'area
C     della tessera con vertici PTS(3,NV). Consideriamo sempre
C     che il lato N della tessera e' quello compreso tra i vertici
C     N e N+1 (oppure NV e 1). In CCC(3,NV) sono le posizioni dei
C     centri degli archi che sottendono i vari lati della tessera.
C     La formula di Gauss-Bonet per le sfere e':
C            Area=R^2[2pi+S(Phi(N)cosT(N))-S(Beta(N))]
C     dove Phi(N) e' la lunghezza d'arco (in radianti) del lato N,
C     T(N) e' l'angolo polare del lato N, Beta(N) l'angolo esterno
C     relativo al vertice N.
C
      TPI=2*PI
C
C     Calcola la prima sommatoria
      SUM1 = D0
      SUMPHI = D0
      DO 100 N = 1, NV
         PHIN(N) = 0.0D0
C
C         IF (NTRHSO(N) .EQ. 1) GOTO 100
C     
         X1 = PTS(1,N) - CCC(1,N)
         Y1 = PTS(2,N) - CCC(2,N)
         Z1 = PTS(3,N) - CCC(3,N)
         IF(N.LT.NV) THEN
            X2 = PTS(1,N+1) - CCC(1,N)
            Y2 = PTS(2,N+1) - CCC(2,N)
            Z2 = PTS(3,N+1) - CCC(3,N)
         ELSE
            X2 = PTS(1,1) - CCC(1,N)
            Y2 = PTS(2,1) - CCC(2,N)
            Z2 = PTS(3,1) - CCC(3,N)
         END IF
         DNORM1 = X1*X1 + Y1*Y1 + Z1*Z1
         DNORM2__ = X2*X2 + Y2*Y2 + Z2*Z2
         SCAL = X1*X2 + Y1*Y2 + Z1*Z2
         COSPHIN = SCAL / (SQRT(DNORM1*DNORM2__))
         IF(COSPHIN.GT.1.0D+00) COSPHIN = 1.0D+00
         PHIN(N) = ACOS(COSPHIN)
         SUMPHI = SUMPHI + PHIN(N)
C     
C     NSFE1 e' la sfera con cui la sfera NS si interseca (eventualmente)
         NSFE1 = INTSPH(NUMTS,N)
         X1 = XE(NSFE1) - XE(NS)
         Y1 = YE(NSFE1) - YE(NS)
         Z1 = ZE(NSFE1) - ZE(NS)
         DNORM1 = SQRT(X1*X1 + Y1*Y1 + Z1*Z1)
         IF(DNORM1.EQ.D0) DNORM1 = 1.0D+00
         X2 = PTS(1,N) - XE(NS)
         Y2 = PTS(2,N) - YE(NS)
         Z2 = PTS(3,N) - ZE(NS)
         DNORM2__ = SQRT(X2*X2 + Y2*Y2 + Z2*Z2)
         COSTN = (X1*X2+Y1*Y2+Z1*Z2)/(DNORM1*DNORM2__)
         SUM1 = SUM1 + PHIN(N) * COSTN
 100  CONTINUE
C     
C WEIGHTS GENERATION FOR CENTER DEFINITION
C 
C the use of old weights is only for backward compatibility
C the new weights give a smooth variation of the 
C tessera change during geometry optimizations especially when
C the number of vertices of a tessera changes.
C
C the difference in the sum has no importance since at the end 
C the point is renormalized in order to be on the sphere
C
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
C
C     Calcola la seconda sommatoria: l'angolo esterno Beta(N) e'
C     definito usando i versori (u(N-1),u(N)) tangenti alla sfera
C     nel vertice N lungo le direzioni dei lati N-1 e N:
C                cos( Pi-Beta(N) )=u(N-1)*u(N)
C            u(N-1) = [V(N) x (V(N) x V(N-1))]/NORM
C            u(N) = [V(N) x (V(N) x V(N+1))]/NORM
C     dove V(I) e' il vettore posizione del vertice I RISPETTO AL
C     CENTRO DELL'ARCO CHE SI STA CONSIDERANDO.
C
      SUM2 = 0.0D+00
C     Loop sui vertici
      DO 200 N = 1, NV
         DO JJ = 1, 3
            P1(JJ) = 0.0D+00
            P2(JJ) = 0.0D+00
            P3(JJ) = 0.0D+00
         ENDDO
         N1 = N
Clf
C         IF(N.GT.1) N0 = N - 1
C         IF(N.EQ.1) N0 = NV
C         IF(N.LT.NV) N2 = N + 1
C         IF(N.EQ.NV) N2 = 1
Clf
C
C If one or both the edges are labelled "small we jump to the preavious and/or next and so on
C until we reach a long enough edge.
C
         N0 = N
Clf 220     CONTINUE
         N0 = MOD(NV+N0-1,NV)
         IF(N0 .EQ. 0) N0 = NV
Clf         write(lvpri,*) "N0 is", N0, NTRHSO(N0)
Clf         IF(NTRHSO(N0).EQ.1) GOTO 220
         N2 = N
Clf         NCONT = 0
Clf 230     CONTINUE
         N2 = MOD(N2+1,NV)
         IF(N2 .EQ. 0) N2 = NV
Clf         write(lvpri,*) "N2 is", N2, NTRHSO(N1)
Clf         IF(NTRHSO(N1 + NCONT).EQ.1) THEN
Clf            NCONT = NCONT + 1
Clf            GOTO 230
CLF         END IF
C
C     Trova i vettori posizione rispetto ai centri corrispondenti
C     e i versori tangenti
C     
C     Lato N0-N1:
         DO JJ = 1, 3
            P1(JJ) = PTS(JJ,N1) - CCC(JJ,N0)
            P2(JJ) = PTS(JJ,N0) - CCC(JJ,N0)
         ENDDO
C     
         CALL VECP_(P1,P2,P3,DNORM3)
         DO JJ = 1, 3
            P2(JJ) = P3(JJ)
         ENDDO
         CALL VECP_(P1,P2,P3,DNORM3)
         DO JJ = 1, 3
            U1(JJ) = P3(JJ)/DNORM3
         ENDDO
C     
C     Lato N1-N2:
         DO JJ = 1, 3
            P1(JJ) = PTS(JJ,N1) - CCC(JJ,N1)
            P2(JJ) = PTS(JJ,N2) - CCC(JJ,N1)
         ENDDO
C     
         CALL VECP_(P1,P2,P3,DNORM3)
         DO JJ = 1, 3
            P2(JJ) = P3(JJ)
         ENDDO
         CALL VECP_(P1,P2,P3,DNORM3)
         DO JJ = 1, 3
            U2(JJ) = P3(JJ)/DNORM3
         ENDDO
C     
         BETA(N) = ACOS(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
CLF         SUM2 = SUM2 + (PI - BETAN)
 200  CONTINUE
      do I=1,NV
CLF         II = I - 1
Clf         IF (II .EQ. 0) II = NV
Clf         SUM2 = SUM2 - BETA(I) * (D1 - DP5 * DBLE(NTRHSO(I)+NTRHSO(II)))
Clf     $               + PI * (D1 - DBLE(NTRHSO(I)))
         SUM2 = SUM2 - BETA(I)
      enddo
C     Calcola l'area della tessera
      AREA = RE(NS)*RE(NS)*(DBLE(2-NV) * PI + SUM1 - SUM2)
C     Trova il punto rappresentativo (come media dei vertici)
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
C     Trova il punto sulla normale (interna!) distante DR dal punto
C     rappresentativo
      PP1(1) = XE(NS) + (PP(1) - XE(NS)) * (RE(NS) - DR) / RE(NS)
      PP1(2) = YE(NS) + (PP(2) - YE(NS)) * (RE(NS) - DR) / RE(NS)
      PP1(3) = ZE(NS) + (PP(3) - ZE(NS)) * (RE(NS) - DR) / RE(NS)
C
C     A causa delle approssimazioni numeriche, l'area di alcune piccole
C     tessere puo' risultare negativa, e viene in questo caso trascurata
      IF(AREA.LT.0.0D+00) THEN
         WRITE(LVPRI,1000) NS,area
         AREA = 0.0D+00
 1000    FORMAT(/,'WARNING: THE AEREA OF ONE TESSERA_ ON SPHERE ',I3,
     *  ' IS IGNORED',F16.10)
      END IF
      RETURN
      END
C/* Deck vecp */
      SUBROUTINE VECP_(P1,P2,P3,DNORM3)
C
#include <pcm_implicit.h>
C
      DIMENSION P1(3),P2(3),P3(3)
C
C     Esegue il prodotto vettoriale P3 = P1 x P2
C
      P3(1) = P1(2)*P2(3) - P1(3)*P2(2)
      P3(2) = P1(3)*P2(1) - P1(1)*P2(3)
      P3(3) = P1(1)*P2(2) - P1(2)*P2(1)
      DNORM3 = SQRT(P3(1)*P3(1) + P3(2)*P3(2) + P3(3)*P3(3))
      RETURN
      END
C/* Deck cavspl */
      SUBROUTINE CAVSPL_(ICAV1,ICAV2,NCAV1,NCAV2,NATM,SOME)
C
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_mxcent.h>
#include <pcm_pcmdef.h>
#include <pcm_infpri.h>
C
      LOGICAL SOME
C
      DIMENSION ICAV1(MXCENT),ICAV2(MXCENT)
C
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
C
C  QUESTA ROUTINE CONTROLLA SE IL SOLUTO E' CONTENUTO IN UNA
C  UNICA CAVITA' O IN CAVITA' DISTINTE: IN TAL CASO I NUMERI
C  D'ORDINE DELLE SFERE CHE COSTITUISCONO LA PRIMA E LA SECONDA
C  CAVITA' VENGONO MEMORIZZATI RISPETTIVAMENTE IN ICAV1 E ICAV2
C
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
  50  DO I=1,MXCENT
        ICAV2(I)=0
      ENDDO
      N2=1
      N=ICAV1(NN)
      DO 100 ICEN=1,NESF
        IF(ICEN.EQ.N) GO TO 100
        DO I=1,NESF
          IF(ICEN.EQ.ICAV1(I)) GO TO 100
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
        IF(SUM.GT.R) THEN
          ICAV1(N1)=ICEN
          N1=N1+1
          ELSE
          ICAV2(N2)=ICEN
          N2=N2+1
        END IF
 100  CONTINUE
      NN=NN+1
      IF(ICAV1(NN).NE.0) GO TO 50
      NCAV1=NN-1
      NCAV2=NESF-NCAV1
      IF(SOME) THEN
         IF(NCAV2.EQ.0) THEN
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
C
 200  FORMAT(/10X,'THE SOLUTE IS ENCLOSED IN ONE CAVITY')
 300  FORMAT(/10X,'THE SOLUTE IS ENCLOSED IN TWO DISTINCT CAVITIES'/
     1        10X,'OF',I3,3X,'E',I3,3X,'SPHERE(S),  RESPECTIVELY')
 400  FORMAT(/10X,'THE FIRST CAVITY IS FORMED BY SPHERE(S) :'/)
 500  FORMAT(/10X,'THE SECOND CAVITY IS FORMED BY SPHERE(S) :'/)
      END
C/* Deck plotcav */
      SUBROUTINE PLOTCAV_(Vert,NUMTS)
C
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_iratdef.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_infpri.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
      DIMENSION VERT(NUMTS,10,*),IVTS(MXTS,10)
C
C     Prepare the input file for GeomView (coloured polyhedra)
C
      IF (ICESPH .EQ. 0 .OR. ICESPH .EQ. 2 .OR. ICESPH .EQ. 3) INDEX = 1
      LUCAV = 12121201
      OPEN(UNIT   = LUCAV,
     $     FILE   = 'cavity.off',
     $     ACCESS = 'SEQUENTIAL',
     $     STATUS = 'UNKNOWN')
      Numv = 0
      DO 10 I = 1, NTs
  10    NUMV = NUMV + NVert(i)
      WRITE(LUCAV,500)'COFF'
      WRITE(LUCAV,1000)NUMV, NTs, NumV
      K = 0
      LAST = 0
 100  CONTINUE
      DO 2010 I = 1, NTs
        N=ISphe(I)
        IF(N.ne.Last) WRITE(LUCAV,1500)N
        Last = N
        IF(Index.EQ.1) THEN
           Call COLOUR_(N,C1,C2,C3)
        ELSE
           C1=1.0D0
           C2=1.0D0
           C3=1.0D0
        END IF
        Do 2020 j = 1, NVert(i)
          IVTS(I,J) = K
          K = K + 1
          WRITE(LUCAV,2001)(VERT(I,J,JCORD),JCORD=1,3),
     $         C1,C2,C3,0.75,I
 2020 CONTINUE
 2010 CONTINUE
      Do 20 I = 1, NTs
 20      WRITE(LUCAV,3000)NVert(I),(IVTS(I,J),J=1,NVERT(I))
 500  FORMAT(1x,a)
1000  FORMAT(3i10)
1500  FORMAT('# Sphere number ',i4)
 2001 FORMAT('  ',3f16.9,4f5.2,' # Tess. ',i4)
3000  FORMAT('  ',14i10)
      RETURN
      END
C/* Deck Colour */
      SUBROUTINE COLOUR_(N,C1,C2,C3)
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_pcmdef.h>
#include <pcm_mxcent.h>
#include <pcm_nuclei.h>
#include <pcm_pcmnuclei.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_codata.h>
      Character*20 Col
      Data Col /'                    '/
C
C     Assign tesserae colours for GeomView:
C     Carbon: green, nitrogen: blue, oxygen: red, hydrogen: light blue,
C     others: violet, added spheres: gray.
C
      Delta=1.d-03
      IF(N.gt.NESFP)THEN
         Col = 'Gray'
         Call COLTSS_(Col,C1,C2,C3)
         Return
      END IF
      NUM = 0
      DO I = 1, NUCDEP
         DIFF=( PCMCORD(1,I) - XE(N) )**2
     $       +( PCMCORD(2,I) - YE(N) )**2
     *       +( PCMCORD(3,I) - ZE(N) )**2
         DDIFF=Sqrt(DIFF)
         IF(DDIFF.LE.Delta) THEN
            NUM=INT(PCMCHG(I))
            GOTO 11
         END IF
      END DO
 11   CONTINUE
      IF(NUM.EQ.6)THEN
         COL = 'Green'
      ELSEIF(NUM.EQ.7)THEN
         COL = 'Blue'
      ELSEIF(NUM.EQ.8)THEN
         COL = 'Red'
      ELSEIF(NUM.EQ.1)THEN
         COL = 'Light Blue'
      ELSE
         COL = 'Fuchsia'
      ENDIF
      CALL COLTSS_(Col,C1,C2,C3)
      RETURN
      END
C/* ComdDeck ColTss */
      SUBROUTINE COLTSS_(COLOUR_,C1,C2,C3)
      Implicit Real*8(A-H,O-Z)
C
C     Assigne tesserae colours for GeomView:
C
      Character*20 COLOUR_
C
      IF(COLOUR_.EQ.'White') THEN
        C1 = 1.0
        C2 = 1.0
        C3 = 1.0
      ELSE IF(COLOUR_.EQ.'Gray') THEN
        C1 = 0.66
        C2 = 0.66
        C3 = 0.66
      ELSE IF(COLOUR_.EQ.'Blue'.or.COLOUR_.EQ.'Dark Blue') THEN
        C1 = 0.0 
        C2 = 0.0 
        C3 = 1.0 
      ELSE IF(COLOUR_.EQ.'Light Blue') THEN
        C1 = 0.0 
        C2 = 1.0 
        C3 = 1.0 
      ELSE IF(COLOUR_.EQ.'Green') THEN
        C1 = 0.0 
        C2 = 1.0 
        C3 = 0.0 
      ELSE IF(COLOUR_.EQ.'Yellow') THEN
        C1 = 1.0 
        C2 = 1.0 
        C3 = 0.0 
      ELSE IF(COLOUR_.EQ.'Orange') THEN
        C1 = 1.0 
        C2 = 0.5 
        C3 = 0.0 
      ELSE IF(COLOUR_.EQ.'Violet') THEN
        C1 = 0.6 
        C2 = 0.0 
        C3 = 1.0 
      ELSE IF(COLOUR_.EQ.'Pink'.or.COLOUR_.EQ.'Light Red') THEN
        C1 = 1.0 
        C2 = 0.5 
        C3 = 1.0 
      ELSE IF(COLOUR_.EQ.'Fuchsia') THEN
        C1 = 1.0 
        C2 = 0.0 
        C3 = 1.0 
      ELSE IF(COLOUR_.EQ.'Red'.or.COLOUR_.EQ.'Dark Red') THEN
        C1 = 1.0 
        C2 = 0.0 
        C3 = 0.0  
      ELSE IF(COLOUR_.EQ.'Black') THEN
        C1 = 0.0 
        C2 = 0.0 
        C3 = 0.0 
      ELSE
        WRITE(*,*) '8'
        STOP
        ENDIF
      RETURN
      END
C/*DECK PREREP_*/
      SUBROUTINE PREREP_(NV,NT,ITS,CV,JTR,NVERT,NUMTS)
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_mxcent.h>
C
#include <pcm_symmet.h>
#include <pcm_pgroup.h>
      LOGICAL LSYMOP(0:7)
      DIMENSION JTR(NUMTS,*),CV(NVERT,*)
#include <pcm_ibtfun.h>
C
C      SYMOP(0) = ' E '
C      SYMOP(1) = 'Oyz'
C      SYMOP(2) = 'Oxz'
C      SYMOP(3) = 'C2z'
C      SYMOP(4) = 'Oxy'
C      SYMOP(5) = 'C2y'
C      SYMOP(6) = 'C2x'
C      SYMOP(7) = ' i '
      DO I = 0, 7
         LSYMOP(I) = .FALSE.
      END DO
      LSYMOP(0) = .TRUE.
      IF(GROUP.EQ.'C1 ') THEN
         LSYMOP(1) = .TRUE.
         LSYMOP(2) = .TRUE.
         LSYMOP(4) = .TRUE.
      ELSEIF(GROUP.EQ.'C2 '.OR.GROUP.EQ.'Ci ') THEN
         LSYMOP(1) = .TRUE.
         LSYMOP(2) = .TRUE.
      ELSEIF(GROUP.EQ.'Cs ') THEN
         LSYMOP(2) = .TRUE.
         LSYMOP(4) = .TRUE.
      ELSEIF(GROUP.EQ.'D2 '.OR.GROUP.EQ.'C2v') THEN
         LSYMOP(1) = .TRUE.
      ELSEIF(GROUP.EQ.'C2h') THEN
         LSYMOP(2) = .TRUE.
      ELSEIF(GROUP.NE.'D2h') THEN
         WRITE(*,*) 'Check symmetry group.'
         STOP
      ENDIF
      DO ISYMOP = 1,7
         IF(LSYMOP(ISYMOP)) THEN
c     
c     riproduzione vertici
c     
            DO I = 1,NV
               II = I + NV
               DO K = 1,3
                  CV(II,K) =
     $              PT(IBTAND(ISYMOP,2**(K-1))) * CV(I,K)
               ENDDO
            ENDDO
Clf            WRITE(LVPRI,*) 'I,II,JJ,K,JTR(II,K),JTR(I,K)'
c     
c     riproduzione topologia 
c     
            DO I = 1,NT
               II = I + NT
               JJ = NV
               Do K = 1,3
                  JTR(II,K) = JTR(I,K) + JJ
Clf                  WRITE(LVPRI,*) I,II,JJ,K,JTR(II,K),JTR(I,K)
               ENDDO
            ENDDO
c
c  aggiornamento indici
c
            NT = NT * 2
            NV = NV * 2
            ITS = ITS * 2
Clf            write(lvpri,*) 'nt,nv,its',nt,nv,its
         ENDIF
      END DO
      DO I = 1,ITS
         DO J = 1,3
Clf            WRITE(LVPRI,*) I,JTR(I,J),(CV(JTR(I,J),K),K=1,3)
         ENDDO
      ENDDO
      RETURN
      END
C  /* Deck pcmtns */ 
      SUBROUTINE PCMTNS_(VMAT,GEOM,AMASS,KATOM)
#include <pcm_implicit.h>
#include <pcm_priunit.h>
#include <pcm_mxcent.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_iratdef.h>
#include <pcm_infpri.h>
#include <pcm_pcmdef.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_symmet.h>
#include <pcm_nuclei.h>
#include <pcm_orgcom.h>
      DIMENSION EIGVAL(3),EIGVEC(3,3),GEOM(KATOM,*),AMASS(*),TINERT(3,3)
      DIMENSION ANGMOM(3),OMEGAD(3),EIGINV(3,3),VMAT(3,3),SCAL(3)
      DIMENSION IAX(6),NORDER_(3)
      LOGICAL PLANAR, LINEAR
      DATA IAX/1,2,3,3,2,1/
#include <pcm_ibtfun.h>      
      JATOM = 1
C
C As correctly pointed out by Kenneth, we want two molecules
C like HCCD and HCCH to have the same cavity so NUMIS is always 1
C
C      NUMIS = 1
C      DO 100 IATOM = 1, NUCIND
C         NATTYP = NINT(CHARGE(IATOM))
Clf         NUMIS  = ISOTOP(IATOM)
C         DO 110 ISYMOP = 0, MAXOPR
C            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
C               AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')
C               DO 120 ICOOR = 1, 3
C                  GEOM(JATOM,ICOOR) = PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
C     &                               *CORD(ICOOR,IATOM) - CMXYZ(ICOOR)
C 120           CONTINUE
Clf               WRITE(LVPRI,*) 'NATTYP, NUMIS,JATOM,AMASS',
Clf     $              NATTYP, NUMIS,JATOM,AMASS(JATOM)
C               JATOM = JATOM + 1
C            END IF
C 110     CONTINUE
C 100  CONTINUE
Clf
C      DO I=1,JATOM-1
C         WRITE(LVPRI,*) 'ATOM ',I,' COORD ',(GEOM(I,J),J=1,3),AMASS(I)
C      ENDDO 
Clf
      DO I = 1, NESFP
         AMASS(I) = RIN(I)
         GEOM(I,1)  = XE(I)
         GEOM(I,2)  = YE(I)
         GEOM(I,3)  = ZE(I)
      ENDDO

      ANGMOM(1) = 1.0D0
      ANGMOM(2) = 1.0D0
      ANGMOM(3) = 1.0D0
      CALL PCM_WLKDIN(GEOM,AMASS,NESFP,ANGMOM,TINERT,OMEGAD,EIGVAL,
     &            EIGVEC,.TRUE.,PLANAR,LINEAR)
Clf      WRITE(LVPRI,*) 'INERTIA TENSOR EIGENVECTORS AND EIGENVALUES'
      DO I = 1, 3
         DO J = 1, 3
            EIGINV(I,J) = EIGVEC(J,I)
         ENDDO
Clf         WRITE(LVPRI,*) (EIGINV(I,J),J=1,3), EIGVAL(I), OMEGA(I), QUAD
      ENDDO
C      do i=1,3
C         write(lvpri,*) 'tinert',(tinert(i,j),j=1,3)
C      enddo
      DO I = 1,3
         DO J= 1,3
            VMAT(J,I) = 0.0D0
         ENDDO
      ENDDO
      NOPAX = NROTS + NREFL
      NORDER_(1) = 1
      NORDER_(2) = 2
      NORDER_(3) = 3
      IF (NOPAX.GE.3) THEN
         DO I = 1,3
C            write(lvpri,*) 'iax',i,isop(i),iax(isop(i))
            DO J = 1,3
               DIJ = 0.0D0
               IF (I.EQ.J) DIJ = 1.0D0
               VMAT(J,IAX(ISOP(I))) = DIJ
            ENDDO
         ENDDO
      ELSEIF (NOPAX.GE.1) THEN
         JAX = IAX(ISOP(1))
         SCAL(1) = DABS(EIGINV(1,JAX))
         SCAL(2) = DABS(EIGINV(2,JAX))
         SCAL(3) = DABS(EIGINV(3,JAX))
         NMAX = 1
         DO J = 2,3
            IF (SCAL(J).GT.SCAL(NMAX)) NMAX = J
         ENDDO
         NSHIFT = MOD(NMAX-1,3)
         DO I = 0,2
            K = MOD(I + NSHIFT,3) + 1
            DO J =1,3
               VMAT(I+1,J) = EIGINV(K,J)
            ENDDO
         ENDDO
      ELSEIF (NOPAX.EQ.0) THEN
         DO I = 1,3
            DO J = 1,3
               VMAT(I,J) = EIGINV(I,J)
            ENDDO
         ENDDO
      ELSE
         WRITE(*,*) '9'
         STOP
      ENDIF
      RETURN
      END
ClfC/* Deck ANLINT */
Clf      SUBROUTINE ANLINT(H0,K0,A0,B0,P0,Q0,O0,RP,RQ,L,NEXIT)
ClfC
ClfC     Luca Frediani February 2 2002. This is a replacemement for the old
ClfC     inter subroutine with an analytic algorithm instead of the old
ClfC     numeric algorithm. It is necessary in order to get a good
ClfC     symmetric cavity.
ClfC      
ClfC
Clf#include <pcm_implicit.h>
Clf#include <pcm_priunit.h>
ClfC
Clf      DOUBLE PRECISION MAT(3,3),W(3),V(3)
Clf      DOUBLE PRECISION H(3), A(3), B(3), P(3), Q(3)
Clf      DOUBLE PRECISION H0(3),A0(3),B0(3),P0(3),Q0(3),O0(3)
Clf      LOGICAL INT1, INT2, RIGHT
ClfC
ClfC A: First point on the initial arc
ClfC B: Second point on the initial arc
ClfC P: Center of the sphere where the arc is located (radii = RP)
ClfC Q: Center of the intersecting sphere (radii = RQ)
ClfC O: Center of the arc
ClfC W: Scratch vector
ClfC
ClfC
ClfClf We copy the original values to the working variables
ClfC   and we translate them at the same time
ClfC
Clf      DO I =1,3
Clf         A(I) = A0(I) - O0(I)
Clf         B(I) = B0(I) - O0(I)
Clf         P(I) = P0(I) - O0(I)
Clf         Q(I) = Q0(I) - O0(I)
Clf      ENDDO
ClfC
ClfC Determination of the normal vector to the arc
ClfC
Clf      W(1) = A(2) * B(3) - A(3) * B(2)
Clf      W(2) = A(3) * B(1) - A(1) * B(3)
Clf      W(3) = A(1) * B(2) - A(2) * B(1)
Clf      write(lvpri,*) 'normal vector', (w(i),i=1,3)
Clf      WM2  = W(1) ** 2 + W(2) ** 2 + W(3) ** 2
Clf      WM   = DSQRT ( WM2 )
Clf      WMXY = DSQRT ( W(1) ** 2 + W(2) ** 2 )
ClfC
ClfC Magnitudes and square magnitudes
ClfC
Clf      PM2 = P(1) ** 2 + P(2) ** 2 + P(3) ** 2
Clf      QM2 = Q(1) ** 2 + Q(2) ** 2 + Q(3) ** 2
Clf      AM2 = A(1) ** 2 + A(2) ** 2 + A(3) ** 2
Clf      BM2 = B(1) ** 2 + B(2) ** 2 + B(3) ** 2
Clf      PM  = DSQRT ( PM2 )
Clf      QM  = DSQRT ( QM2 )
Clf      AM  = DSQRT ( AM2 )
Clf      BM  = DSQRT ( BM2 )
ClfC      write(lvpri,*) 'moduli'
ClfC      write(lvpri,*) pm2,qm2,am2,bm2
ClfC      write(lvpri,*) pm,qm,am,bm
ClfC
ClfClf determination of sines and cosines of the euler angles
ClfC
Clf      IF(WMXY.GT.1.0D-12) THEN
Clf         write(lvpri,*) 'wmxy grande'
Clf         COSA = W(1) / DSQRT( W(1) ** 2 + W(2) ** 2 )
Clf         SINA = W(2) / DSQRT( W(1) ** 2 + W(2) ** 2 )
Clf         COSB = W(3) / WM
Clf         SINB = DSQRT( 1 - COSB ** 2 )
Clf      ELSEIF(WM.GT.1.0D-8) THEN
Clf         write(lvpri,*) 'wm grande'
Clf         COSA = 1.0D0
Clf         SINA = 0.0D0
Clf         COSB = W(3)/DABS(W(3))
Clf         SINB = 0.0D0
Clf      ELSE
Clf         CALL QUIT('ANLINT: NORMAL VECTOR TOO SMALL')
Clf      ENDIF
ClfC
Clf      V(1) = COSB * COSA * W(1) + COSB * SINA * W(2) - SINB * W(3)
Clf      V(2) =      - SINA * W(1) +        COSA * W(2)              
Clf      V(3) = SINB * COSA * W(1) + SINB * SINA * W(2) + COSB * W(3)
Clf      IF(DABS(V(1)).GT.1.0D-8.OR.DABS(V(2)).GT.1.0D-8) 
Clf     $     CALL QUIT('ANLINT: ROTATION ERROR 1')
Clf      IF(V(3).LE.0.0D0) THEN 
Clf         WRITE(LVPRI,*) 'ROT ERR 2', (V(I),I=1,3)
Clf         CALL QUIT('ANLINT: ROTATION ERROR 2')
Clf      ENDIF
ClfC
Clf      W(1) = COSB * COSA * A(1) + COSB * SINA * A(2) - SINB * A(3)
Clf      W(2) =      - SINA * A(1) +        COSA * A(2)              
Clf      W(3) = SINB * COSA * A(1) + SINB * SINA * A(2) + COSB * A(3)
Clf      COSG = W(1) / AM
Clf      SING = W(2) / AM
Clf      IF(DABS(W(3)).GT.1.0D-8) CALL QUIT('ANLINT: ROTATION ERROR 3')
Clf      write(lvpri,*) 'DOPO a e b rot',(w(i),i=1,3)
Clf      write(lvpri,*) 'rp e rq', rp,rq
Clf      write(lvpri,*) 'sin e cos'
Clf      write(lvpri,*) sina,sinb,sing
Clf      write(lvpri,*) cosa,cosb,cosg
ClfC
ClfClf rotation matrix
ClfC
Clf      MAT(1,1) =   COSG * COSB * COSA - SING * SINA
Clf      MAT(1,2) =   COSG * COSB * SINA + SING * COSA
Clf      MAT(1,3) = - COSG * SINB
Clf      MAT(2,1) = - SING * COSB * COSA - COSG * SINA
Clf      MAT(2,2) = - SING * COSB * SINA + COSG * COSA
Clf      MAT(2,3) =   SING * SINB
Clf      MAT(3,1) =   SINB * COSA
Clf      MAT(3,2) =   SINB * SINA
Clf      MAT(3,3) =   COSB
Clf      write(lvpri,*)'mat'
Clf      do i=1,3
Clf         write(lvpri,*) (mat(I,J),J=1,3)
Clf      enddo
ClfC
ClfClf rotation of vectors and intersection calculations
ClfC
Clf      write(lvpri,*)'AAAAAA',AM2
Clf      CALL TRIROT(MAT,A)
Clf      write(lvpri,*)'BBBBBB',BM2
Clf      CALL TRIROT(MAT,B)
Clf      write(lvpri,*)'PPPPPP',PM2
Clf      CALL TRIROT(MAT,P)
Clf      write(lvpri,*)'QQQQQQ',QM2
Clf      CALL TRIROT(MAT,Q)
Clf      R1R1 = RP**2 - PM2
Clf      R1 = DSQRT(R1R1)
Clf      COST = (RP**2 - RQ**2 + QM2 - PM2)/2.0D0
Clf      IF(Q(1)**2 + Q(2)**2.LT.1.0D-20) THEN
Clf         CALL QUIT('NO HOPE!')
Clf      END IF
Clf      IF(DABS(Q(1)) .GE. DABS(Q(2))) THEN
Clf         CALL SECGRA_(COST,R1,Q(1),Q(2),X1,Y1,X2,Y2,NSOL)
Clf      ELSE
Clf         CALL SECGRA_(COST,R1,Q(2),Q(1),Y1,X1,Y2,X2,NSOL)
Clf      END IF
Clf      WRITE(LVPRI,*) 'SOL1',X1,Y1,X1**2+Y1**2
Clf      WRITE(LVPRI,*) 'SOL2',X2,Y2,X2**2+Y2**2
Clf      IF(NSOL.EQ.-1) THEN
Clf         CALL QUIT('ANLINT: SECGRA_ EXITED ABNORMALLY')
Clf      ELSE IF (NSOL.EQ.0) THEN
Clf         WRITE(LVPRI,*) 
Clf     $        'NO INTER_SECTION BETWEEN THIS EDGE AND THIS SPHERE'
ClfC
ClfC     PUT SOME CODE HERE!!!
ClfC
Clf         RETURN
Clf      ELSE IF (NSOL.EQ.1) THEN
Clf         WRITE(LVPRI,*) 'THE TWO SPHERES ARE (ALMOST) TANGENT'
ClfC
ClfC     PUT SOME CODE HERE!!!
ClfC
Clf      ELSE IF (NSOL.EQ.2) 
Clf         NCASE = 0
Clf         THRESH = 1.0D-12
Clf         DA1 = DSQRT((X1 - A(1))**2 + (Y1-A(2))**2)
Clf         DB1 = DSQRT((X1 - B(1))**2 + (Y1-B(2))**2)
Clf         DA2 = DSQRT((X2 - A(1))**2 + (Y2-A(2))**2)
Clf         DB2 = DSQRT((X2 - B(1))**2 + (Y2-B(2))**2)
ClfC
ClfC Check some inconsistencies
ClfC
Clf         IF ((DA1.LT.THRESH.AND.DB1.LT.THRESH).OR.
Clf     $       (DA2.LT.THRESH.AND.DB2.LT.THRESH) THEN
Clf            CALL QUIT('DISTANCE INCONSISTENCY IN ANLINT')
Clf         ENDIF
Clf         INT1 = (X1 .LT. A(1)) .AND. (X1 .GT. B(1)).AND.
Clf     $          (Y1 .GT. A(2)) .AND. (Y1 .LT. B(2))
Clf         INT2 = (X2 .LT. A(1)) .AND. (X2 .GT. B(1)).AND.
Clf     $          (Y2 .GT. A(2)) .AND. (Y2 .LT. B(2))
Clf         IF(DB1.LT.THRESH) THEN
Clf            NCASE = NCASE + 3
Clf         ELSE IF(DA1.LT.THRESH) THEN
Clf            NCASE = NCASE + 2
Clf         ELSE IF(.NOT.INT1)
Clf            NCASE = NCASE + 1
Clf         END IF
Clf         IF(DB2.LT.THRESH) THEN
Clf            NCASE = NCASE + 12
Clf         ELSE IF(DA2.LT.THRESH) THEN
Clf            NCASE = NCASE + 8
Clf         ELSE IF(.NOT.INT2)
Clf            NCASE = NCASE + 4
Clf         END IF
Clf         WRITE(LVPRI,*) 'NCASE =' NCASE
Clf         GO TO (100,101,102,103,104,105,106,107,108,109,
Clf     $          110,111,112,113,114,115), NCASE
Clf 100     
Clf      RIGHT = INT1.XOR.INT2
Clf      IF(RIGHT) THEN
Clf         H(3) = 0.0D0
Clf         IF(INT1) THEN
Clf            H(1) = X1
Clf            H(2) = Y1
Clf         ELSE
Clf            H(1) = X2
Clf            H(2) = Y2
Clf         ENDIF
Clf      ELSE
Clf         CALL QUIT('ANLINT: INTER_SECTION ERROR 2')
Clf      ENDIF         
ClfC
ClfC     transpose matrix elements
ClfC
Clf      MAT(2,1) =   COSG * COSB * SINA + SING * COSA
Clf      MAT(3,1) = - COSG * SINB
Clf      MAT(1,2) = - SING * COSB * COSA - COSG * SINA
Clf      MAT(3,2) =   SING * SINB
Clf      MAT(1,3) =   SINB * COSA
Clf      MAT(2,3) =   SINB * SINA
Clf      write(lvpri,*)'matinv'
Clf      do i=1,3
Clf         write(lvpri,*) (mat(I,J),J=1,3)
Clf      enddo
Clf      CALL TRIROT(MAT,H)
Clf      DO I = 1,3
Clf         H0(I) = H(I) + O0(I) 
Clf      ENDDO
Clf      RETURN
Clf      END
ClfC
ClfC
ClfC/* Deck TRIROT*/
Clf      SUBROUTINE TRIROT(M,V)
Clf#include <pcm_implicit.h>
Clf#include <pcm_priunit.h>
Clf      DOUBLE PRECISION M(3,*),V(*),W(3)
Clf      write(lvpri,*) 'PRIMA', (V(I),I=1,3)
Clf      DO I = 1,3
Clf         W(I) = M(I,1) * V(1) + M(I,2) * V(2) + M(I,3) * V(3)
Clf      ENDDO
Clf      DO I = 1,3
Clf         V(I) = W(I)
Clf      ENDDO
Clf      write(lvpri,*) 'DOPO', (V(I),I=1,3)
Clf      RETURN
Clf      END

C/* Deck UPDCAV_*/
      SUBROUTINE UPDCAV_(COORD)
#include <pcm_implicit.h>
#include <pcm_mxcent.h>
#include <pcm_maxaqn.h>
#include <pcm_maxorb.h>
#include <pcm_pcmdef.h>
#include <pcm_pcm.h>
#include <pcm_pcmlog.h>
#include <pcm_nuclei.h>
#include <pcm_pcmnuclei.h>
#include <pcm_symmet.h>
      DIMENSION COORD(3,*)
#include <pcm_ibtfun.h>
      IF (ICESPH .EQ. 0 .OR. ICESPH .EQ. 2) THEN
         JATOM = 0
         DO I = 1, NUCIND
            MULCNT = ISTBNU(I)
            DO ISYM = 0, MAXREP
               IF(IBTAND(ISYM,MULCNT) .EQ. 0) THEN
                  JATOM = JATOM + 1
                  DO JCORD = 1,3
                     PCMCORD(JCORD,JATOM) = 
     $                    PT(IBTAND(ISYMAX(JCORD,1),ISYM)) 
     $                    *COORD(JCORD,I)
                  END DO
               END IF
            END DO
         END DO
      END IF
      RETURN
      END
C/* Deck ordpcm*/
      SUBROUTINE ORDPCM_(lvpri,nts,xtscor,ytscor,ztscor,as,privec,
     $     idxpri,work,lwork)
#include <pcm_implicit.h>
      DIMENSION PRIVEC(4,*), XTSCOR(*), YTSCOR(*), ZTSCOR(*),
     &          WORK(*), AS(*)
      DIMENSION IDXPRI(*)
      LOGICAL CHKTSS_, LCHK, LSWTCH
      LSWTCH = .FALSE.
      DO I = 1,NTS
         PRIVEC(1,I) = XTSCOR(I)
         PRIVEC(2,I) = YTSCOR(I)
         PRIVEC(3,I) = ZTSCOR(I)
         PRIVEC(4,I) = AS(I)
         IDXPRI(I) = I 
      ENDDO
 222  CONTINUE
      DO I = 1,NTS - 1
         II = I + 1
         LCHK = CHKTSS_(PRIVEC(1,I),PRIVEC(2,I),PRIVEC(3,I),PRIVEC(4,I),
     $        PRIVEC(1,II),PRIVEC(2,II),PRIVEC(3,II),PRIVEC(4,II))
         LSWTCH = LSWTCH .OR. LCHK
         IF (LCHK) THEN
            XBAK         = PRIVEC(1,I) 
            YBAK         = PRIVEC(2,I) 
            ZBAK         = PRIVEC(3,I)  
            ABAK         = PRIVEC(4,I) 
            IBAK         = IDXPRI(I)
            PRIVEC(1,I)  = PRIVEC(1,II) 
            PRIVEC(2,I)  = PRIVEC(2,II) 
            PRIVEC(3,I)  = PRIVEC(3,II) 
            PRIVEC(4,I)  = PRIVEC(4,II) 
            IDXPRI(I)    = IDXPRI(II)  
            PRIVEC(1,II) = XBAK             
            PRIVEC(2,II) = YBAK 
            PRIVEC(3,II) = ZBAK 
            PRIVEC(4,II) = ABAK 
            IDXPRI(II)   = IBAK  
         END IF
      ENDDO
      IF(LSWTCH) THEN
         LSWTCH = .FALSE.
         GOTO 222
      ENDIF
      DO I = 1, NTS
         WRITE(LVPRI,1250) (dabs(PRIVEC(J,I)),J=1,4)
      ENDDO
      RETURN
 1240 format(i3,4f15.9)
 1250 format(4f15.9)
      END

C/* Deck chktss*/
      LOGICAL FUNCTION CHKTSS_(X1,Y1,Z1,A1,X2,Y2,Z2,A2)
#include <pcm_implicit.h>
      IF(DABS(X1).LT.DABS(X2)) THEN
         CHKTSS_ = .TRUE.
         RETURN
      ELSE IF(DABS(X1).GT.DABS(X2)) THEN
         CHKTSS_ = .FALSE.
         RETURN
      ELSE IF(DABS(Y1).LT.DABS(Y2)) THEN
         CHKTSS_ = .TRUE.
         RETURN
      ELSE IF(DABS(Y1).GT.DABS(Y2)) THEN
         CHKTSS_ = .FALSE.
         RETURN
      ELSE IF(DABS(Z1).LT.DABS(Z2)) THEN
         CHKTSS_ = .TRUE.
         RETURN
      ELSE IF(DABS(Z1).GT.DABS(Z2)) THEN
         CHKTSS_ = .FALSE.
         RETURN
      ELSE IF(DABS(A1).LT.DABS(A2)) THEN
         CHKTSS_ = .TRUE.
         RETURN
      ELSE
         CHKTSS_ = .FALSE.
         RETURN
      ENDIF
      END
      
C/* Deck SECGRA_*/
      SUBROUTINE SECGRA_(COST,RP,Q1,Q2,X1,Y1,X2,Y2,NSOL)
#include <pcm_implicit.h>
      PARAMETER (D0 =0.0D0)
      THRESH = 1.0D-20
      A = Q1 ** 2 + Q2 ** 2
      B = - 2.0D0 * COST * Q2 
      C = COST ** 2 - RP ** 2 * Q1 ** 2
      DELTA = B ** 2 - 4.0D0 * A * C
      IF(DELTA.GT.THRESH) THEN
         NSOL = 2
         Y1 = (- B + DSQRT(DELTA)) / (2.0D0 * A)
         Y2 = (- B - DSQRT(DELTA)) / (2.0D0 * A)
         X1 = (COST - Q2 * Y1) / Q1
         X2 = (COST - Q2 * Y2) / Q1
         RETURN
      ELSE IF(DELTA.LT.-THRESH) THEN
         NSOL = 0
         X1 = D0
         X2 = D0
         Y1 = D0
         Y2 = D0
         RETURN
      ELSE
         NSOL = 1
         Y1 = - B / (2.0D0 * A)
         Y2 = Y1
         X1 = (COST - Q2 * Y1) / Q1
         X2 = X1
         RETURN
      END IF
      NSOL = - 1
      RETURN
      END
C
C/* Deck TYPLAB_*/
C
      CHARACTER*16 FUNCTION TYPLAB_(I)
      INTEGER I
      IF(I.EQ.4) THEN
         TYPLAB_='ALL EDGE IS FREE'
      ELSE IF(I.EQ.1) THEN
         TYPLAB_='2ND VERT COVERED'
      ELSE IF(I.EQ.2) THEN
         TYPLAB_='1ST VERT COVERED'
      ELSE IF(I.EQ.3) THEN
         TYPLAB_='EDGE  CUT  TWICE'
      ELSE IF(I.EQ.0) THEN
         TYPLAB_='ALL EDGE COVERED'
      ELSE
         TYPLAB_='UNDEFINED CASE!!'
      END IF
      RETURN
      END
