c/* Deck Deriva */
      SUBROUTINE CAVDER_(NSJ,NSJR,ICOORD,INTSPH,NEWSPH)
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
C
      INTEGER ALGE(63),CASCA(10)
      DIMENSION INTSPH(MXTS,10), NEWSPH(MXSP,2)
      PARAMETER (D0=0.0D0)
C
C     Le derivate contengono termini dovuti direttamente allo
C     spostamento del centro della sfera NSJ, e termini "mediati" dagli
C     spostamenti del centro e dal cambiamento del raggio delle sfere
C     "aggiunte" (create da PEDRA_, oltre a quelle originarie).
C
C     Memorizza in DERRAD(NS,NSJ,ICOORD) la derivata del raggio di
C     NS e in DERCEN(NS,NSJ,ICOORD,3) le derivate delle
C     coordinate del centro di NS rispetto alla coord. ICOORD della
C     sfera NSJ.
C
C     Se NS e' una sfera originaria queste derivate sono 0, tranne
C     DERCEN(NSJR,NSJ,ICOORD,ICOORD)=1:
C
C
C
      DERCEN(NSJR,NSJ,ICOORD,ICOORD) = 1.0D+00
C
C     2) Effetti indiretti.
C     Loop sulle sfere aggiunte
C
      DO 500 NSA = NESFP+1, NESF
         DO II = 1, 63
            ALGE(II) = 0
         ENDDO
C
C     Costruiamo l'"albero genealogico" della sfera NSA
C
         ALGE(1) = NSA
         ALGE(2) = ABS(NEWSPH(NSA,1))
         ALGE(3) = ABS(NEWSPH(NSA,2))
         LIVEL = 3
         NUMBER = 2
 510     NSUB = 1
         DO II = LIVEL-NUMBER+1, LIVEL
            IF(ALGE(II).GT.NESFP) THEN
               ALGE(LIVEL+NSUB)   = ABS(NEWSPH(ALGE(II),1))
               ALGE(LIVEL+NSUB+1) = ABS(NEWSPH(ALGE(II),2))
            END IF
            NSUB = NSUB + 2
         ENDDO
         NUMBER = NUMBER * 2
         LIVEL = LIVEL + NUMBER
         IF(NUMBER.LT.32) GO TO 510
C
C     Quando un elemento di ALGE e' = NSJR, costruisce la corrispondente
C     "cascata" di sfere aggiunte che collega NSJR a NSA
C
         DO 600 LIVEL = 2, 6
            MIN = 2**(LIVEL-1)
            MAX = (2**LIVEL) - 1
            DO 700 II = MIN, MAX
               IF(ALGE(II).NE.NSJR) GO TO 700
               DO K = 1, 10
                  CASCA(K) = 0
               ENDDO
               CASCA(1) = NSJR
               INDEX = II
               K = 2
               DO LL = LIVEL, 2, -1
                  FACT = (INDEX - 2**(LL-1)) / 2.0D+00
                  INDEX = INT(2**(LL-2) + FACT)
                  CASCA(K) = ALGE(INDEX)
                  K = K + 1
               ENDDO
C     Contiamo gli elementi diversi da 0 in CASCA
               ICONT = 0
               DO K = 1, 10
                  IF(CASCA(K).NE.0) ICONT = ICONT + 1
               ENDDO
C
C     Costruiamo le derivate composte del raggio e delle coordinate di
C     NSA (ultimo elemento di CASCA)
C     rispetto alla coordinata ICOORD di NSJ (primo elemento di CASCA)
C
               NS1 = CASCA(1)
               NS2 = CASCA(2)
               CALL DRRDCN_(NS2,ICOORD,NS1,DR1,NEWSPH)
               CALL DRCNCN_(1,NS2,ICOORD,NS1,DX,NEWSPH)
               CALL DRCNCN_(2,NS2,ICOORD,NS1,DY,NEWSPH)
               CALL DRCNCN_(3,NS2,ICOORD,NS1,DZ,NEWSPH)
               DO 800 I = 3, ICONT
                  DDR = D0
                  DDX = D0
                  DDY = D0
                  DDZ = D0
                  NS1 = CASCA(I-1)
                  NS2 = CASCA(I)
C     
C     Derivata del raggio dell'elemento I di CASCA rispetto
C     alla coord. ICOORD dell'elemento 1 di CASCA
C
                  CALL DRRDRD_(NS2,NS1,DER,NEWSPH)
                  DDR = DER * DR1
                  CALL DRRDCN_(NS2,1,NS1,DER,NEWSPH)
                  DDR = DDR + DER * DX
                  CALL DRRDCN_(NS2,2,NS1,DER,NEWSPH)
                  DDR = DDR + DER * DY
                  CALL DRRDCN_(NS2,3,NS1,DER,NEWSPH)
                  DDR = DDR + DER * DZ
C     
C     Derivata della coord. X dell'elemento I di CASCA rispetto
C     alla coord. ICOORD dell'elemento 1 di CASCA
C
                  CALL DRCNRD_(1,NS2,NS1,DER,NEWSPH)
                  DDX = DER * DR1
                  CALL DRCNCN_(1,NS2,1,NS1,DER,NEWSPH)
                  DDX = DDX + DER * DX
                  CALL DRCNCN_(1,NS2,2,NS1,DER,NEWSPH)
                  DDX = DDX + DER * DY
                  CALL DRCNCN_(1,NS2,3,NS1,DER,NEWSPH)
                  DDX = DDX + DER * DZ
C     
C     Derivata della coord. Y dell'elemento I di CASCA rispetto
C     alla coord. ICOORD dell'elemento 1 di CASCA
C     
                  CALL DRCNRD_(2,NS2,NS1,DER,NEWSPH)
                  DDY = DER * DR1
                  CALL DRCNCN_(2,NS2,1,NS1,DER,NEWSPH)
                  DDY = DDY + DER * DX
                  CALL DRCNCN_(2,NS2,2,NS1,DER,NEWSPH)
                  DDY = DDY + DER * DY
                  CALL DRCNCN_(2,NS2,3,NS1,DER,NEWSPH)
                  DDY = DDY + DER * DZ
C     
C     Derivata della coord. Z dell'elemento I di CASCA rispetto
C     alla coord. ICOORD dell'elemento 1 di CASCA
C     
                  CALL DRCNRD_(3,NS2,NS1,DER,NEWSPH)
                  DDZ = DER * DR1
                  CALL DRCNCN_(3,NS2,1,NS1,DER,NEWSPH)
                  DDZ = DDZ + DER * DX
                  CALL DRCNCN_(3,NS2,2,NS1,DER,NEWSPH)
                  DDZ = DDZ + DER * DY
                  CALL DRCNCN_(3,NS2,3,NS1,DER,NEWSPH)
                  DDZ = DDZ + DER * DZ
C     
                  DR1 = DDR
                  DX = DDX
                  DY = DDY
                  DZ = DDZ
 800           CONTINUE
C     
C     Se NS e' una sfera aggiunta, memorizza le derivate del raggio
C     e delle coordinate del centro:
C     
               DERRAD(NSA,NSJ,ICOORD) = DR1
               DERCEN(NSA,NSJ,ICOORD,1) = DX
               DERCEN(NSA,NSJ,ICOORD,2) = DY
               DERCEN(NSA,NSJ,ICOORD,3) = DZ
 700        CONTINUE
 600     CONTINUE
 500  CONTINUE
C     
      RETURN
      END
C***********************************************************************
c/* Deck Dcdr */
      SUBROUTINE DRCNRD_(JJ,NSI,NSJ,DC,NEWSPH)
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
      DIMENSION COORDJ(3), COORDK(3)
      DIMENSION INTSPH(MXTS,10), NEWSPH(MXSP,2)
      PARAMETER (D0=0.D0)
C
C     Trova la derivata della coordinata JJ del centro della sfera
C     NSI rispetto al raggio dellla sfera NSJ.
C
C     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA_)
C     dipende dalle due sfere "precedenti" NSJ e NSK
C
C     Se NSJ o NSK sono negativi, la sfera aggiunta e' di tipo C
C     e la generatrice "principale" corrisponde al label negativo
C     (cfr. JCC 11, 1047 (1990))
C
      IF(NEWSPH(NSI,1).LT.0.OR.NEWSPH(NSI,2).LT.0) GO TO 100
      NSK = NEWSPH(NSI,1)
      IF(NSK.EQ.NSJ) NSK = NEWSPH(NSI,2)
      COORDJ(1) = XE(NSJ)
      COORDJ(2) = YE(NSJ)
      COORDJ(3) = ZE(NSJ)
      COORDK(1) = XE(NSK)
      COORDK(2) = YE(NSK)
      COORDK(3) = ZE(NSK)
      D2 = (XE(NSJ)-XE(NSK))**2 + (YE(NSJ)-YE(NSK))**2 +
     *     (ZE(NSJ)-ZE(NSK))**2
      D = SQRT(D2)
      DC = - (COORDJ(JJ) - COORDK(JJ)) / (2.0D+00 * D)
      GO TO 200
C
 100  CONTINUE
      NSK = NEWSPH(NSI,1)
      IF(ABS(NSK).EQ.NSJ) NSK = NEWSPH(NSI,2)
      DC = D0
      IF(NSK.LT.D0) GO TO 200
      COORDJ(1) = XE(NSJ)
      COORDJ(2) = YE(NSJ)
      COORDJ(3) = ZE(NSJ)
      COORDK(1) = XE(NSK)
      COORDK(2) = YE(NSK)
      COORDK(3) = ZE(NSK)
      D2 = (XE(NSJ)-XE(NSK))**2 + (YE(NSJ)-YE(NSK))**2 +
     *     (ZE(NSJ)-ZE(NSK))**2
      D = SQRT(D2)
      DC = - ( COORDJ(JJ) - COORDK(JJ) ) / D
C
 200  CONTINUE
      RETURN
      END
c/* Deck Dcdc */
      SUBROUTINE DRCNCN_(JJ,NSI,ICOORD,NSJ,DC,NEWSPH)
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
      DIMENSION COORDJ(3), COORDK(3)
      DIMENSION NEWSPH(MXSP,2)
      PARAMETER (D0=0.D0)
C
C     Trova la derivata della coordinata JJ del centro della sfera
C     NSI rispetto alla coordinata ICOORD di NSJ, che interseca NSI.
C
C     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA_)
C     dipende dalle due sfere "precedenti" NSJ e NSK
C
C     Se NSJ o NSK sono negativi, la sfera aggiunta e' di tipo C
C     e la generatrice "principale" corrisponde al label negativo
C     (cfr. JCC 11, 1047 (1990))
C
      IF(NEWSPH(NSI,1).LT.0.OR.NEWSPH(NSI,2).LT.0) GO TO 100
      K = NEWSPH(NSI,1)
      IF(K.EQ.NSJ) K = NEWSPH(NSI,2)
      COORDJ(1) = XE(NSJ)
      COORDJ(2) = YE(NSJ)
      COORDJ(3) = ZE(NSJ)
      COORDK(1) = XE(K)
      COORDK(2) = YE(K)
      COORDK(3) = ZE(K)
      D2 = (XE(NSJ)-XE(K))**2 + (YE(NSJ)-YE(K))**2 + (ZE(NSJ)-ZE(K))**2
      D = SQRT(D2)
      DC = (RE(NSJ)-RE(K)) * (COORDJ(ICOORD)-COORDK(ICOORD)) *
     *        (COORDJ(JJ) - COORDK(JJ)) / (2.0D+00 * D**3)
      IF(JJ.EQ.ICOORD)DC = DC + 0.5D+00 - (RE(NSJ)-RE(K)) / (2.0D+00*D)
      GO TO 200
C
 100  CONTINUE
      NSK = NEWSPH(NSI,1)
      IF(ABS(NSK).EQ.NSJ) NSK = NEWSPH(NSI,2)
      COORDJ(1) = XE(NSJ)
      COORDJ(2) = YE(NSJ)
      COORDJ(3) = ZE(NSJ)
      COORDK(1) = XE(ABS(NSK))
      COORDK(2) = YE(ABS(NSK))
      COORDK(3) = ZE(ABS(NSK))
      D2 = (COORDJ(1)-COORDK(1))**2 + (COORDJ(2)-COORDK(2))**2 +
     *     (COORDJ(3)-COORDK(3))**2
      D = SQRT(D2)
      IF(NSK.GT.0) THEN
        DC = RE(NSJ) * (COORDJ(JJ)-COORDK(JJ)) * (COORDJ(ICOORD)-
     *  COORDK(ICOORD)) / D**3
        IF(ICOORD.EQ.JJ) DC = DC + 1.0D+00 - RE(NSJ) / D
      ELSE
        DC = - RE(ABS(NSK)) * (COORDK(JJ)-COORDJ(JJ)) * (COORDK(ICOORD)-
     *  COORDJ(ICOORD)) / D**3
        IF(ICOORD.EQ.JJ) DC = DC + RE(ABS(NSK)) / D
      END IF
C
 200  CONTINUE
      RETURN
      END
c/* Deck Drdr */
      SUBROUTINE DRRDRD_(NSI,NSJ,DR1,NEWSPH)
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
      DIMENSION NEWSPH(MXSP,2)
      PARAMETER (D0=0.D0)
C
C     Trova la derivata del raggio della sfera NSI rispetto al raggio
C     della sfera NSJ.
C
C     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA_)
C     dipende dalle due sfere "precedenti" NSJ e NSK
C     Se NSJ o NSK sono negativi, la sfera aggiunta e' di tipo C
C     e la generatrice "principale" corrisponde al label negativo
C     (cfr. JCC 11, 1047 (1990))
C
      IF(NEWSPH(NSI,1).LT.0 .OR. NEWSPH(NSI,2).LT.0) GO TO 100
      NSK = NEWSPH(NSI,1)
      IF(NSK.EQ.NSJ) NSK = NEWSPH(NSI,2)
      RS = RSOLV
      RJ = RE(NSJ) + RS
      RK = RE(NSK) + RS
      RI = RE(NSI) + RS
      D2 = (XE(NSJ)-XE(NSK))**2 + (YE(NSJ)-YE(NSK))**2 +
     *   (ZE(NSJ)-ZE(NSK))**2
      D = SQRT(D2)
      DR1 = (-3.0D+00*RJ*RJ + RK*RK + 2.0D+00*RJ*RK
     *      + 3.0D+00*D*RJ - D*RK) / (4.0D+00*D*RI)
      GO TO 200
C
 100  CONTINUE
      NSK = NEWSPH(NSI,1)
      IF(ABS(NSK).EQ.NSJ) NSK = NEWSPH(NSI,2)
C
      IF(NSK.GT.0) THEN
        RS = RSOLV
        RJ = RE(NSJ) + RS
        RK = RE(NSK) + RS
        RI = RE(NSI) + RS
        D2 = (XE(NSJ)-XE(NSK))**2 + (YE(NSJ)-YE(NSK))**2 +
     *     (ZE(NSJ)-ZE(NSK))**2
        D = SQRT(D2)
        DR1 = ( 2.0D+00*D*RJ + 2.0D+00*D*RE(NSJ) - 2.0D+00*RJ*RE(NSJ) +
     *       D*D - RJ*RJ - RK*RK ) / (2.0D+00*D*RI)
      ELSE
        RS = RSOLV
        RJ = RE(NSJ) + RS
        RI = RE(NSI) + RS
        D2 = (XE(NSJ)-XE(ABS(NSK)))**2 + (YE(NSJ)-YE(ABS(NSK)))**2 +
     *     (ZE(NSJ)-ZE(ABS(NSK)))**2
        D = SQRT(D2)
        DR1 = ( RE(ABS(NSK)) * RJ ) / ( D*RI)
      END IF
 200  CONTINUE
      RETURN
      END
c/* Deck Drdc */
      SUBROUTINE DRRDCN_(NSI,ICOORD,NSJ,DR1,NEWSPH)
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
      DIMENSION NEWSPH(MXSP,2)
      DIMENSION COORDJ(3),COORDK(3)
      PARAMETER (D0=0.D0)
C
C     Trova la derivata del raggio della sfera NSI rispetto alla
C     coordinata ICOORD (1=X, 2=Y, 3=Z) della sfera NSJ, che interseca
C     NSI.
C
C     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA_)
C     dipende dalle due sfere "precedenti" NSJ e K
C
C     Se NSJ o NSK sono negativi, la sfera aggiunta e' di tipo C
C     e la generatrice "principale" corrisponde al label negativo
C     (cfr. JCC 11, 1047 (1990))
C
      IF(NEWSPH(NSI,1).LT.0.OR.NEWSPH(NSI,2).LT.0) GO TO 100
      K = NEWSPH(NSI,1)
      IF(K.EQ.NSJ) K = NEWSPH(NSI,2)
      COORDJ(1) = XE(NSJ)
      COORDJ(2) = YE(NSJ)
      COORDJ(3) = ZE(NSJ)
      COORDK(1) = XE(K)
      COORDK(2) = YE(K)
      COORDK(3) = ZE(K)
      D2 = (XE(NSJ)-XE(K))**2 + (YE(NSJ)-YE(K))**2 + (ZE(NSJ)-ZE(K))**2
      D = SQRT(D2)
      B = 0.5D+00 * (D + RE(NSJ) - RE(K))
      RS = RSOLV
      A = ((RE(NSJ)+RS)**2 + D2 - (RE(K)+RS)**2) / D
      DR1 = (2.0D+00*A*B - 2.0D+00*B*D - A*D) *
     *     (COORDJ(ICOORD)-COORDK(ICOORD)) / (4.0D+00*D2*(RE(NSI)+RS))
      GO TO 200
C
 100  CONTINUE
      NSK = NEWSPH(NSI,1)
      IF(ABS(NSK).EQ.NSJ) NSK = NEWSPH(NSI,2)
      COORDJ(1) = XE(NSJ)
      COORDJ(2) = YE(NSJ)
      COORDJ(3) = ZE(NSJ)
      COORDK(1) = XE(ABS(NSK))
      COORDK(2) = YE(ABS(NSK))
      COORDK(3) = ZE(ABS(NSK))
      RI = RE(NSI) + RSOLV
      RJ = RE(NSJ) + RSOLV
      RK = RE(ABS(NSK)) + RSOLV
      DIFF = COORDJ(ICOORD) - COORDK(ICOORD)
      D2 = (COORDJ(1)-COORDK(1))**2 + (COORDJ(2)-COORDK(2))**2 +
     *     (COORDJ(3)-COORDK(3))**2
      D = SQRT(D2)
      FAC = RE(NSJ) * ( RJ*RJ - D*D - RK*RK )
      IF(NSK.LT.0) FAC = RE(ABS(NSK)) * (RK*RK - D*D - RJ*RJ )
      DR1 = DIFF * FAC / ( 2.0D+00 * D**3 * RI )
C
 200  CONTINUE
      RETURN
      END
