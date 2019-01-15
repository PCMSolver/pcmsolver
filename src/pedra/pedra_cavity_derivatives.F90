!
! PCMSolver, an API for the Polarizable Continuum Model
! Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
!
! This file is part of PCMSolver.
!
! PCMSolver is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PCMSolver is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
!
! For information on the complete list of contributors to the
! PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
!

    module pedra_cavity_derivatives

    use pedra_precision

    implicit none

    public cavder

    private

    contains

    subroutine cavder(nsj, nsjr, icoord, intsph, newsph)

#include "pcm_mxcent.inc"
#include "pcm_nuclei.inc"
#include "pcm_pcmdef.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: nsj, nsjr, icoord
    integer(kind=regint_k) :: intsph(mxts, 10), newsph(mxsp, 2)

    real(kind=dp), parameter :: d0 = 0.0d0
    integer(kind=regint_k) :: alge(63), casca(10)
    real(kind=dp) :: ddr, ddx, ddy, ddz, der, dr1, dx, dy, dz
    real(kind=dp) :: fact
    integer(kind=regint_k) :: i, ii, index, k, livel, ll, max, icont, min
    integer(kind=regint_k) :: ns1, ns2, nsa, nsub, number


!     Le derivate contengono termini dovuti direttamente allo
!     spostamento del centro della sfera NSJ, e termini "mediati" dagli
!     spostamenti del centro e dal cambiamento del raggio delle sfere
!     "aggiunte" (create da PEDRA, oltre a quelle originarie).

!     Memorizza in DERRAD(NS,NSJ,ICOORD) la derivata del raggio di
!     NS e in DERCEN(NS,NSJ,ICOORD,3) le derivate delle
!     coordinate del centro di NS rispetto alla coord. ICOORD della
!     sfera NSJ.

!     Se NS e' una sfera originaria queste derivate sono 0, tranne
!     DERCEN(NSJR,NSJ,ICOORD,ICOORD)=1:



    DERCEN(NSJR,NSJ,ICOORD,ICOORD) = 1.0D+00

!     2) Effetti indiretti.
!     Loop sulle sfere aggiunte

    DO 500 NSA = NESFP+1, NESF
        DO II = 1, 63
            ALGE(II) = 0
        end do

    !     Costruiamo l'"albero genealogico" della sfera NSA

        ALGE(1) = NSA
        ALGE(2) = ABS(NEWSPH(NSA,1))
        ALGE(3) = ABS(NEWSPH(NSA,2))
        LIVEL = 3
        NUMBER = 2
        510 NSUB = 1
        DO II = LIVEL-NUMBER+1, LIVEL
            IF(ALGE(II) > NESFP) THEN
                ALGE(LIVEL+NSUB)   = ABS(NEWSPH(ALGE(II),1))
                ALGE(LIVEL+NSUB+1) = ABS(NEWSPH(ALGE(II),2))
            end if
            NSUB = NSUB + 2
        end do
        NUMBER = NUMBER * 2
        LIVEL = LIVEL + NUMBER
        IF(NUMBER < 32) go to 510

    !     Quando un elemento di ALGE e' = NSJR, costruisce la corrispondente
    !     "cascata" di sfere aggiunte che collega NSJR a NSA

        DO 600 LIVEL = 2, 6
            MIN = 2**(LIVEL-1)
            MAX = (2**LIVEL) - 1
            DO 700 II = MIN, MAX
                IF(ALGE(II) /= NSJR) go to 700
                DO K = 1, 10
                    CASCA(K) = 0
                end do
                CASCA(1) = NSJR
                INDEX = II
                K = 2
                DO LL = LIVEL, 2, -1
                    FACT = (INDEX - 2**(LL-1)) / 2.0D+00
                    INDEX = INT(2**(LL-2) + FACT)
                    CASCA(K) = ALGE(INDEX)
                    K = K + 1
                end do
            !     Contiamo gli elementi diversi da 0 in CASCA
                ICONT = 0
                DO K = 1, 10
                    IF(CASCA(K) /= 0) ICONT = ICONT + 1
                end do

            !     Costruiamo le derivate composte del raggio e delle coordinate di
            !     NSA (ultimo elemento di CASCA)
            !     rispetto alla coordinata ICOORD di NSJ (primo elemento di CASCA)

                NS1 = CASCA(1)
                NS2 = CASCA(2)
                CALL DRRDCN(NS2,ICOORD,NS1,DR1,NEWSPH)
                CALL DRCNCN(1,NS2,ICOORD,NS1,DX,NEWSPH)
                CALL DRCNCN(2,NS2,ICOORD,NS1,DY,NEWSPH)
                CALL DRCNCN(3,NS2,ICOORD,NS1,DZ,NEWSPH)
                DO 800 I = 3, ICONT
                    DDR = D0
                    DDX = D0
                    DDY = D0
                    DDZ = D0
                    NS1 = CASCA(I-1)
                    NS2 = CASCA(I)

                !     Derivata del raggio dell'elemento I di CASCA rispetto
                !     alla coord. ICOORD dell'elemento 1 di CASCA

                    CALL DRRDRD(NS2,NS1,DER,NEWSPH)
                    DDR = DER * DR1
                    CALL DRRDCN(NS2,1,NS1,DER,NEWSPH)
                    DDR = DDR + DER * DX
                    CALL DRRDCN(NS2,2,NS1,DER,NEWSPH)
                    DDR = DDR + DER * DY
                    CALL DRRDCN(NS2,3,NS1,DER,NEWSPH)
                    DDR = DDR + DER * DZ

                !     Derivata della coord. X dell'elemento I di CASCA rispetto
                !     alla coord. ICOORD dell'elemento 1 di CASCA

                    CALL DRCNRD(1,NS2,NS1,DER,NEWSPH)
                    DDX = DER * DR1
                    CALL DRCNCN(1,NS2,1,NS1,DER,NEWSPH)
                    DDX = DDX + DER * DX
                    CALL DRCNCN(1,NS2,2,NS1,DER,NEWSPH)
                    DDX = DDX + DER * DY
                    CALL DRCNCN(1,NS2,3,NS1,DER,NEWSPH)
                    DDX = DDX + DER * DZ

                !     Derivata della coord. Y dell'elemento I di CASCA rispetto
                !     alla coord. ICOORD dell'elemento 1 di CASCA

                    CALL DRCNRD(2,NS2,NS1,DER,NEWSPH)
                    DDY = DER * DR1
                    CALL DRCNCN(2,NS2,1,NS1,DER,NEWSPH)
                    DDY = DDY + DER * DX
                    CALL DRCNCN(2,NS2,2,NS1,DER,NEWSPH)
                    DDY = DDY + DER * DY
                    CALL DRCNCN(2,NS2,3,NS1,DER,NEWSPH)
                    DDY = DDY + DER * DZ

                !     Derivata della coord. Z dell'elemento I di CASCA rispetto
                !     alla coord. ICOORD dell'elemento 1 di CASCA

                    CALL DRCNRD(3,NS2,NS1,DER,NEWSPH)
                    DDZ = DER * DR1
                    CALL DRCNCN(3,NS2,1,NS1,DER,NEWSPH)
                    DDZ = DDZ + DER * DX
                    CALL DRCNCN(3,NS2,2,NS1,DER,NEWSPH)
                    DDZ = DDZ + DER * DY
                    CALL DRCNCN(3,NS2,3,NS1,DER,NEWSPH)
                    DDZ = DDZ + DER * DZ

                    DR1 = DDR
                    DX = DDX
                    DY = DDY
                    DZ = DDZ
                800 end do

            !     Se NS e' una sfera aggiunta, memorizza le derivate del raggio
            !     e delle coordinate del centro:

                DERRAD(NSA,NSJ,ICOORD) = DR1
                DERCEN(NSA,NSJ,ICOORD,1) = DX
                DERCEN(NSA,NSJ,ICOORD,2) = DY
                DERCEN(NSA,NSJ,ICOORD,3) = DZ
            700 end do
        600 end do
    500 end do

    end subroutine cavder

    subroutine drcnrd(jj,nsi,nsj,dc,newsph)

#include "pcm_mxcent.inc"
#include "pcm_nuclei.inc"
#include "pcm_pcmdef.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: jj, nsi, nsj
    real(kind=dp) :: coordj(3), coordk(3)
    integer(kind=regint_k) :: intsph(mxts, 10), newsph(mxsp, 2)

    real(kind=dp), parameter :: d0 = 0.0d0
    real(kind=dp) :: dc, d, d2
    integer(kind=regint_k) :: nsk

!     Trova la derivata della coordinata JJ del centro della sfera
!     NSI rispetto al raggio dellla sfera NSJ.

!     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA_)
!     dipende dalle due sfere "precedenti" NSJ e NSK

!     Se NSJ o NSK sono negativi, la sfera aggiunta e' di tipo C
!     e la generatrice "principale" corrisponde al label negativo
!     (cfr. JCC 11, 1047 (1990))

    IF(NEWSPH(NSI,1) < 0 .OR. NEWSPH(NSI,2) < 0) go to 100
    NSK = NEWSPH(NSI,1)
    IF(NSK == NSJ) NSK = NEWSPH(NSI,2)
    COORDJ(1) = XE(NSJ)
    COORDJ(2) = YE(NSJ)
    COORDJ(3) = ZE(NSJ)
    COORDK(1) = XE(NSK)
    COORDK(2) = YE(NSK)
    COORDK(3) = ZE(NSK)
    D2 = (XE(NSJ)-XE(NSK))**2 + (YE(NSJ)-YE(NSK))**2 + &
    (ZE(NSJ)-ZE(NSK))**2
    D = SQRT(D2)
    DC = - (COORDJ(JJ) - COORDK(JJ)) / (2.0D+00 * D)
    go to 200

    100 CONTINUE
    NSK = NEWSPH(NSI,1)
    IF(ABS(NSK) == NSJ) NSK = NEWSPH(NSI,2)
    DC = D0
    IF(NSK < D0) go to 200
    COORDJ(1) = XE(NSJ)
    COORDJ(2) = YE(NSJ)
    COORDJ(3) = ZE(NSJ)
    COORDK(1) = XE(NSK)
    COORDK(2) = YE(NSK)
    COORDK(3) = ZE(NSK)
    D2 = (XE(NSJ)-XE(NSK))**2 + (YE(NSJ)-YE(NSK))**2 + &
    (ZE(NSJ)-ZE(NSK))**2
    D = SQRT(D2)
    DC = - ( COORDJ(JJ) - COORDK(JJ) ) / D

    200 CONTINUE
    end subroutine drcnrd

    subroutine drcncn(jj,nsi,icoord,nsj,dc,newsph)

#include "pcm_mxcent.inc"
#include "pcm_nuclei.inc"
#include "pcm_pcmdef.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: jj, nsi, icoord, nsj
    real(kind=dp) :: dc
    integer(kind=regint_k) :: newsph(mxsp,2)

    real(kind=dp) :: coordj(3), coordk(3)
    real(kind=dp), parameter :: d0 = 0.0d0
    real(kind=dp) :: d, d2
    integer(kind=regint_k) :: k, nsk

!     Trova la derivata della coordinata JJ del centro della sfera
!     NSI rispetto alla coordinata ICOORD di NSJ, che interseca NSI.

!     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA_)
!     dipende dalle due sfere "precedenti" NSJ e NSK

!     Se NSJ o NSK sono negativi, la sfera aggiunta e' di tipo C
!     e la generatrice "principale" corrisponde al label negativo
!     (cfr. JCC 11, 1047 (1990))

    IF(NEWSPH(NSI,1) < 0 .OR. NEWSPH(NSI,2) < 0) go to 100
    K = NEWSPH(NSI,1)
    IF(K == NSJ) K = NEWSPH(NSI,2)
    COORDJ(1) = XE(NSJ)
    COORDJ(2) = YE(NSJ)
    COORDJ(3) = ZE(NSJ)
    COORDK(1) = XE(K)
    COORDK(2) = YE(K)
    COORDK(3) = ZE(K)
    D2 = (XE(NSJ)-XE(K))**2 + (YE(NSJ)-YE(K))**2 + (ZE(NSJ)-ZE(K))**2
    D = SQRT(D2)
    DC = (RE(NSJ)-RE(K)) * (COORDJ(ICOORD)-COORDK(ICOORD)) * &
    (COORDJ(JJ) - COORDK(JJ)) / (2.0D+00 * D**3)
    IF(JJ == ICOORD)DC = DC + 0.5D+00 - (RE(NSJ)-RE(K)) / (2.0D+00*D)
    go to 200

    100 CONTINUE
    NSK = NEWSPH(NSI,1)
    IF(ABS(NSK) == NSJ) NSK = NEWSPH(NSI,2)
    COORDJ(1) = XE(NSJ)
    COORDJ(2) = YE(NSJ)
    COORDJ(3) = ZE(NSJ)
    COORDK(1) = XE(ABS(NSK))
    COORDK(2) = YE(ABS(NSK))
    COORDK(3) = ZE(ABS(NSK))
    D2 = (COORDJ(1)-COORDK(1))**2 + (COORDJ(2)-COORDK(2))**2 + &
    (COORDJ(3)-COORDK(3))**2
    D = SQRT(D2)
    IF(NSK > 0) THEN
        DC = RE(NSJ) * (COORDJ(JJ)-COORDK(JJ)) * (COORDJ(ICOORD)- &
        COORDK(ICOORD)) / D**3
        IF(ICOORD == JJ) DC = DC + 1.0D+00 - RE(NSJ) / D
    ELSE
        DC = - RE(ABS(NSK)) * (COORDK(JJ)-COORDJ(JJ)) * (COORDK(ICOORD)- &
        COORDJ(ICOORD)) / D**3
        IF(ICOORD == JJ) DC = DC + RE(ABS(NSK)) / D
    end if

    200 CONTINUE
    end subroutine drcncn

    subroutine drrdrd(nsi,nsj,dr1,newsph)

#include "pcm_mxcent.inc"
#include "pcm_nuclei.inc"
#include "pcm_pcmdef.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: nsi, nsj, newsph(mxsp, 2)
    real(kind=dp) :: dr1

    real(kind=dp), parameter :: d0 = 0.0d0
    real(kind=dp) :: d, d2, ri, rj, rk, rs
    integer(kind=regint_k) :: nsk

!     Trova la derivata del raggio della sfera NSI rispetto al raggio
!     della sfera NSJ.

!     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA_)
!     dipende dalle due sfere "precedenti" NSJ e NSK
!     Se NSJ o NSK sono negativi, la sfera aggiunta e' di tipo C
!     e la generatrice "principale" corrisponde al label negativo
!     (cfr. JCC 11, 1047 (1990))

    IF(NEWSPH(NSI,1) < 0 .OR. NEWSPH(NSI,2) < 0) go to 100
    NSK = NEWSPH(NSI,1)
    IF(NSK == NSJ) NSK = NEWSPH(NSI,2)
    RS = RSOLV
    RJ = RE(NSJ) + RS
    RK = RE(NSK) + RS
    RI = RE(NSI) + RS
    D2 = (XE(NSJ)-XE(NSK))**2 + (YE(NSJ)-YE(NSK))**2 + &
    (ZE(NSJ)-ZE(NSK))**2
    D = SQRT(D2)
    DR1 = (-3.0D+00*RJ*RJ + RK*RK + 2.0D+00*RJ*RK &
    + 3.0D+00*D*RJ - D*RK) / (4.0D+00*D*RI)
    go to 200

    100 CONTINUE
    NSK = NEWSPH(NSI,1)
    IF(ABS(NSK) == NSJ) NSK = NEWSPH(NSI,2)

    IF(NSK > 0) THEN
        RS = RSOLV
        RJ = RE(NSJ) + RS
        RK = RE(NSK) + RS
        RI = RE(NSI) + RS
        D2 = (XE(NSJ)-XE(NSK))**2 + (YE(NSJ)-YE(NSK))**2 + &
        (ZE(NSJ)-ZE(NSK))**2
        D = SQRT(D2)
        DR1 = ( 2.0D+00*D*RJ + 2.0D+00*D*RE(NSJ) - 2.0D+00*RJ*RE(NSJ) + &
        D*D - RJ*RJ - RK*RK ) / (2.0D+00*D*RI)
    ELSE
        RS = RSOLV
        RJ = RE(NSJ) + RS
        RI = RE(NSI) + RS
        D2 = (XE(NSJ)-XE(ABS(NSK)))**2 + (YE(NSJ)-YE(ABS(NSK)))**2 + &
        (ZE(NSJ)-ZE(ABS(NSK)))**2
        D = SQRT(D2)
        DR1 = ( RE(ABS(NSK)) * RJ ) / ( D*RI)
    end if
    200 CONTINUE
    end subroutine drrdrd

    subroutine drrdcn(nsi, icoord, nsj, dr1, newsph)

#include "pcm_mxcent.inc"
#include "pcm_nuclei.inc"
#include "pcm_pcmdef.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: nsi, icoord, nsj, newsph(mxsp, 2)
    real(kind=dp) :: dr1

    real(kind=dp) :: coordj(3), coordk(3)
    real(kind=dp), parameter :: d0 = 0.0d0
    real(kind=dp) :: a, b, d, d2, diff, fac, ri, rj, rk, rs
    integer(kind=regint_k) :: k, nsk

!     Trova la derivata del raggio della sfera NSI rispetto alla
!     coordinata ICOORD (1=X, 2=Y, 3=Z) della sfera NSJ, che interseca
!     NSI.

!     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA_)
!     dipende dalle due sfere "precedenti" NSJ e K

!     Se NSJ o NSK sono negativi, la sfera aggiunta e' di tipo C
!     e la generatrice "principale" corrisponde al label negativo
!     (cfr. JCC 11, 1047 (1990))

    IF(NEWSPH(NSI,1) < 0 .OR. NEWSPH(NSI,2) < 0) go to 100
    K = NEWSPH(NSI,1)
    IF(K == NSJ) K = NEWSPH(NSI,2)
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
    DR1 = (2.0D+00*A*B - 2.0D+00*B*D - A*D) * &
    (COORDJ(ICOORD)-COORDK(ICOORD)) / (4.0D+00*D2*(RE(NSI)+RS))
    go to 200

    100 CONTINUE
    NSK = NEWSPH(NSI,1)
    IF(ABS(NSK) == NSJ) NSK = NEWSPH(NSI,2)
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
    D2 = (COORDJ(1)-COORDK(1))**2 + (COORDJ(2)-COORDK(2))**2 + &
    (COORDJ(3)-COORDK(3))**2
    D = SQRT(D2)
    FAC = RE(NSJ) * ( RJ*RJ - D*D - RK*RK )
    IF(NSK < 0) FAC = RE(ABS(NSK)) * (RK*RK - D*D - RJ*RJ )
    DR1 = DIFF * FAC / ( 2.0D+00 * D**3 * RI )

    200 CONTINUE
    end subroutine drrdcn

    end module pedra_cavity_derivatives
