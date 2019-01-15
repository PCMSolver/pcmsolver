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

! Originally written for DALTON by Luca Frediani (ca. 2003-2004)
! Extracted from DALTON by Krzysztof Mozgawa and Ville Weijo (ca. 2010-2012)
! RDR 0114 Wrap up in a F90 module, with a major clean-up of all
!          the DALTON stuff that still lingered but was unused.
!
! NOTE: the ibtfun module is gone, as we can safely use Fortran
!       standard intrisic functions.
!       Mapping between intrinsics and ibtfun:
!           ibtand(i, j) <-> iand(i, j)
!           ibtor(i, j)  <-> ior(i, j)
!           ibtshl(i, j) <-> ishft(i, j)
!           ibtshr(i, j) <-> ishft(i, -j) !WARNING!
!           ibtxor(i, j) <-> ieor(i, j)
    module pedra_cavity

    use pedra_precision
    use pedra_symmetry, only: point_group

    implicit none

    public polyhedra_driver

    private

    ! The point group
    type(point_group) :: group
    ! Some print levels
    integer(kind=regint_k) :: iprsol
    ! The global print unit
    integer(kind=regint_k) :: lvpri
    ! Error code
    !    0: everything went OK
    !    1: exceeded maximum number of spheres
    !       [subroutine polyhedra]
    !    2: exceeded maximum number of tesserae
    !       [subroutine polyhedra]
    !    3: number of spheres bigger than number of nuclei [only for gradients]
    !       [subroutine polyhedra]
    !    4: cavity not consistent with symmetry
    !       [subroutine sphper]
    !    5: additional spheres not consistent with symmetry
    !       [subroutine sphper]
    !    6: itsnum > mxts ?
    !       [subroutine polygen]
    !    7: nts > mxts ?
    !       [subroutine repcav]
    !    8: too many vertices in tessera
    !       [subroutine tessera]
    !    9: ?
    !       [subroutine tessera]
    !   10: ?
    !       [subroutine inter]
    !   11: symmetry string not recognized
    !       [subroutine prerep]
    integer(kind=regint_k) :: pedra_error_code

    contains

    subroutine polyhedra_driver(pgroup, vert, centr, masses, global_print_unit, error_code)

#include "pcm_pcmdef.inc"
#include "pcm_mxcent.inc"
#include "pcm_pcm.inc"

    type(point_group) :: pgroup
    real(kind=dp)     :: masses(mxts)
    integer(kind=regint_k)           :: global_print_unit
    integer(kind=regint_k)           :: error_code

    logical :: some
    integer(kind=regint_k) :: numts, numsph, natm, numver

    integer(kind=regint_k), allocatable :: intsph(:, :), newsph(:, :)
    integer(kind=regint_k), allocatable :: icav1(:), icav2(:)
    integer(kind=regint_k), allocatable :: jtr(:, :)
    real(kind=dp), allocatable :: vert(:, :, :), centr(:, :, :)
    real(kind=dp), allocatable :: xval(:), yval(:), zval(:)
    real(kind=dp), allocatable :: cv(:, :)

    SOME = IPRPCM.NE.-5
!     ----- set the print level
    iprsol = 0
!     ----- set lvpri
    lvpri = global_print_unit
!     ----- set the point group
    group = pgroup

!     maximum number of tessalations is not known, take worst case
!     maximum number of spheres is not known, take worst case

    numts  = mxts
    numsph = mxsp
    natm   = mxcent
    numver = mxver

    write(lvpri, '(a)') "Memory management through standard Fortran 90 allocate/deallocate."

    allocate(intsph(numts, 10))
    intsph = 0
    allocate(newsph(numsph, 2))
    newsph = 0
    allocate(icav1(natm))
    icav1 = 0
    allocate(icav2(natm))
    icav2 = 0
    allocate(jtr(numts, 3))
    jtr = 0
    allocate(xval(numts))
    xval = 0.0d0
    allocate(yval(numts))
    yval = 0.0d0
    allocate(zval(numts))
    zval = 0.0d0
    allocate(cv(numver, 3))
    cv = 0.0d0

    call polyhedra(intsph, vert, centr, newsph, icav1, icav2, xval, yval, zval, &
    jtr, cv, numts, numsph, numver, natm, some, masses)

! Bring the error code back home
    error_code = pedra_error_code

    deallocate(intsph)
    deallocate(newsph)
    deallocate(icav1)
    deallocate(icav2)
    deallocate(jtr)
    deallocate(xval)
    deallocate(yval)
    deallocate(zval)
    deallocate(cv)

    end subroutine polyhedra_driver

    subroutine polyhedra(intsph,vert,centr,newsph,icav1,icav2,xval,yval, &
    zval,jtr,cv,numts,numsph,numver,natm,some,masses)
!
! We list variables of interest, that might be declared inside a common block.
! nesf: total number of spheres
! nesfp: number of original spheres
!

    use pedra_dblas, only: dzero
    use pedra_print, only: output
!    use pedra_cavity_derivatives, only: cavder

#include "pcm_pcmdef.inc"
#include "pcm_mxcent.inc"
#include "pcm_pcm.inc"
! Import nucdep from pcm_nuclei.inc: to set up gradient calculation...
#include "pcm_nuclei.inc"

    integer(kind=regint_k) :: numts, natm, numsph, numver
    integer(kind=regint_k) :: intsph(numts, 10), newsph(numsph, 2), icav1(natm), icav2(natm)
    real(kind=dp) :: masses(mxts)
    real(kind=dp) :: vert(numts, 10, 3), centr(numts, 10, 3), cv(numver, 3)
    real(kind=dp) :: xval(numts), yval(numts), zval(numts)
    real(kind=dp) :: pp(3), pp1(3), pts(3, 10), ccc(3, 10)
    logical :: some
    integer(kind=regint_k) :: jtr(numts, 3)

    real(kind=dp), parameter :: d0 = 0.0d0
    real(kind=dp) :: area, cosom2, fc, fc1, hh, omg, prod, r2gn
    real(kind=dp) :: reg, reg2, regd2, ren, rend2, reo, reo2
    real(kind=dp) :: rep, rep2, repd2, rgn, rij, rij2, rik, rik2
    real(kind=dp) :: rjd, rjk, rjk2, rtdd, rtdd2, senom, sp
    real(kind=dp) :: test, test1, test2, test3, test7, test8
    real(kind=dp) :: xen, xi, xj, xn, yen, yi, yj, yn, zen, zi, zj, zn
    integer(kind=regint_k) :: i, idisp, ii, ipflag, iprcav, iptype
    integer(kind=regint_k) :: its, itseff, itsnum, itypc, iv, iver, j, jj, jcor
    integer(kind=regint_k) :: k, kg, idisrep
    integer(kind=regint_k) :: kp, n, n1, n2, n3
    integer(kind=regint_k) :: ncav1, ncav2, ne, nes, net, nev, nn
    integer(kind=regint_k) :: nsfe, nv
    real(kind=dp) :: rotcav(3, 3)
    real(kind=dp), allocatable :: mass(:), geom(:, :), ssfe(:)
    integer(kind=regint_k), allocatable :: permutation_table(:, :)
    integer(kind=regint_k) :: untesselated

!     Se stiamo costruendo una nuova cavita' per il calcolo del
!     contributo dispersivo e repulsivo:

    idisrep = 0
    iprcav = 0
    rotcav = reshape([1.0d0, 0.0d0, 0.0d0,  &
                      0.0d0, 1.0d0, 0.0d0,  &
                      0.0d0, 0.0d0, 1.0d0], [3, 3])


! se icesph=0
!     legge i raggi dall'input e fa coincidere i centri con gli atomi
! se icesph=1
!     legge centri e raggi dall'input
! se icesph=2
!     legge i raggi dall'input e fa coincidere i centri con
!     alcuni atomi definiti dall'indice ina(k) con k=1,NESFP
!     es: xe(k)=c(1,ina(k))
! icesph is always 1, centers and radii are passed explicitly from
!  the module. PEDRA is completely ignorant about the existence of molecules.
! All data is already in bohr. No unit conversions are performed by PEDRA

    write(lvpri, '(a)') "Input spheres"
    write(lvpri, '(a, i5)') "Initial number of spheres = ", nesfp
    write(lvpri, '(a)') "Sphere                      Center  (X,Y,Z) (AU)                       Radius (AU)"
    do i = 1, nesfp
        re(i) = rin(i) ! The scaling has been done C++-side
        write(lvpri, '(i4, 4f20.14)') i, xe(i), ye(i), ze(i), re(i)
    end do

! Creation of new spheres
    do n = 1, nesf
        newsph(n,1) = 0
        newsph(n,2) = 0
    end do

    itypc = 0
    ! Degrees-to-Radians conversion: pi/180
    omg = omega * (acos(-1.0d0) / 180.0d0)
    senom = sin(omg)
    cosom2 = (cos(omg))**2
    rtdd = ret + rsolv
    rtdd2 = rtdd * rtdd
    net = nesf
    nn = 2
    ne = nesf
    nev = nesf
    go to 100
    110 NN = NE + 1
    NE = NET
    100 CONTINUE

!     check on the number of spheres

    write(lvpri, *) 'Number of extra spheres = ', ne
    if (ne > mxsp) then
        write(lvpri, *) ne
        write(lvpri, *) 'Exceeded limit on the maximum number of spheres:', mxsp
        pedra_error_code = 1
        stop
    end if

! Addition of spheres happens in the following
! Criteria for addition of spheres are listed and explained in:
! Paper 1:
! J. L. Pascual-Ahuir, E. Silla, J. Tomasi, R. Bonaccorsi, J. Comput. Chem. 8, 778 (1987)
! Paper 2:
! J. L. Pascual-Ahuir, E. Silla, J. Comput. Chem. 11, 1047 (1990)
    do 120 i = nn, ne
        nes = i - 1
        do 130 j=1,nes
            rij2 = (xe(i) - xe(j))**2 + (ye(i) - ye(j))**2 + (ze(i) - ze(j))**2
            rij = sqrt(rij2)
            rjd = re(j) + rsolv
            test1 = re(i) + rjd + rsolv
            if (rij >= test1) cycle
            reg = max(re(i), re(j))
            rep = min(re(i), re(j))
            reg2 = reg * reg
            rep2 = rep * rep
            test2 = rep * senom + sqrt(reg2 - rep2 * cosom2)
            if (rij <= test2) cycle
            regd2 = (reg + rsolv) * (reg + rsolv)
            test3 = (regd2 + reg2 - rtdd2) / reg
            if (rij >= test3) cycle
            do k = 1, nev
                if (k == j .or. k == i) cycle
                rjk2 = (xe(j) - xe(k))**2 + (ye(j) - ye(k))**2 + (ze(j) - ze(k))**2
                if (rjk2 >= rij2) cycle
                rik2 = (xe(i) - xe(k))**2 + (ye(i) - ye(k))**2 + (ze(i) - ze(k))**2
                if (rik2 >= rij2) cycle
                rjk = sqrt(rjk2)
                rik = sqrt(rik2)
                sp = (rij + rjk + rik) / 2.0d0
                hh = 4 * (sp * (sp - rij) * (sp - rik) * (sp - rjk)) / rij2
                reo = re(k) * fro
                if (k >= ne) reo = 0.0002d0
                reo2 = reo * reo
                if (hh < reo2) go to 130
            end do
            repd2 = (rep + rsolv)**2
            test8 = sqrt(repd2 - rtdd2) + sqrt(regd2 - rtdd2)
            if (rij <= test8) go to 150
            rend2 = regd2 + reg2 - (reg / rij) * (regd2 + rij2 - repd2)
            if (rend2 <= rtdd2) go to 130
            ren = sqrt(rend2) - rsolv
            fc = reg / (rij - reg)
            test7 = reg - re(i)
            kg = i
            kp = j
            if (test7 <= 1.0d-09) go to 160
            kg = j
            kp = i
            160 fc1 = fc + 1.0
            xen = (xe(kg) + fc * xe(kp)) / fc1
            yen = (ye(kg) + fc * ye(kp)) / fc1
            zen = (ze(kg) + fc * ze(kp)) / fc1
            itypc = 1
            go to 170
            150 r2gn = rij - rep + reg
            rgn = r2gn / 2.0d0
            fc = r2gn / (rij + rep - reg)
            fc1 = fc + 1.0d0
            test7 = reg - re(i)
            kg = i
            kp = j
            if (test7 <= 1.0d-09) go to 180
            kg = j
            kp = i
            180 xen = (xe(kg) + fc * xe(kp)) / fc1
            yen = (ye(kg) + fc * ye(kp)) / fc1
            zen = (ze(kg) + fc * ze(kp)) / fc1
            ren = sqrt(regd2 + rgn * (rgn - (regd2 + rij2 - repd2) / rij)) - rsolv
            170 net = net + 1
            xe(net) = xen
            ye(net) = yen
            ze(net) = zen
            re(net) = ren

            ! The matrix newsph(nesf, 2) stores the numbers of spheres
            ! "generating" the new sphere with index net: if the new sphere is
            ! of type A or B both numbers are positive. If the new sphere is of
            ! type C the number of the "principal" sphere is negative.
            ! Definition of the new sphere types are given in Paper 2.
            if(itypc == 0) then
                newsph(net, 1) = kg
                newsph(net, 2) = kp
            elseif (itypc == 1) then
                newsph(net, 1) = - kg
                newsph(net, 2) = kp
            end if

        130 end do
        nev = net
    120 end do
    if (net /= ne) go to 110
    nesf = net

! Build the spheres permutation table
! Allocate the space for the permutation table and clean it up
    allocate(permutation_table(nesf, (group%maxrep+1)))
    permutation_table = 0
    call sphper(nesf, nesfp, numsph, permutation_table, newsph, xe, ye, ze, re)

! Determination of eigenvalues and eigenvectors of the inertia tensor.
! We have to construct the tensor of inertia and diagonalize it.
! Our tensor of inertia uses the center of the spheres as coordinates
! and the radii as masses.

    allocate(mass(nesfp))
    allocate(geom(nesfp, 3))

    mass = 0.0d0
    geom = 0.0d0

    do i = 1, nesfp
        mass(i)   = masses(i)
        geom(i,1) = xe(i)
        geom(i,2) = ye(i)
        geom(i,3) = ze(i)
    end do

    if (nesfp > 1) then
        call pcmtns(rotcav, geom, mass, nesfp)
    end if
    deallocate(geom)
    deallocate(mass)

! Division of the surface into tesserae
! For every tessera, check if there's other tesserae covering it
! or not. If yes, cut it.
    NN = 0
    do nsfe = 1, nesf
        ! First of all check if the sphere is Unique For Symmetry (UFS)
        ! and hence needs to be tesselated.
        ! If the sphere is not UFS, then go to the next sphere.
        if (permutation_table(nsfe, 1) == 0) then
                write(lvpri, '(a, i4, a)') "Sphere ", nsfe, " is not Unique-For-Symmetry: skipping tesselation."
!                write(lvpri, '(a, i4, a, i2)') "permutation_table(", nsfe, ", 1) = ", permutation_table(nsfe, 1)
                cycle
        end if

        ! Transfer center and radius to temporaries
        xen = xe(nsfe)
        yen = ye(nsfe)
        zen = ze(nsfe)
        ren = re(nsfe)
        if(iprcav >= 10) then
            write(lvpri, '(a, i4)') "Now tesselating sphere no. ", nsfe
            write(lvpri, '(a, f12.8, a, f12.8, a, f12.8, a, f12.8)') "Coordinates: X = ", xen, &
                            " Y = ", yen, " Z = ",  zen, " R = ", ren
        end if
        ! Default options (same as traditional GEPOL)
        iptype = 2
        ipflag = 0
        itsnum = 60
    ! f Which type of tessellation?
    !         IF(IPOLYG.GT.0) THEN
    !            IPFlag = 0
    !            ITSNUM = IPolyg
    !         ELSEIF(IPolyg.lt.0) THEN
    !            IPFlag = 1
    !            TsAre = -1.0D-03*IPolyg
    ! f         end if
        ! polygen: where the magic happens. Generate the spherical polyhedra.
        call polygen(ipflag, areats, itsnum, xen, yen, zen, ren, itseff, cv,  &
        jtr, permutation_table, nesf, nsfe, numts, numver, rotcav)
        if (iprcav >= 10) then
              write(lvpri, '(a, i8)') "After polygen, itseff=", itseff
        end if
        do 310 its = 1, itseff
            n1 = jtr(its, 1)
            n2 = jtr(its, 2)
            n3 = jtr(its, 3)
            pts(1, 1) = cv(n1, 1)
            pts(2, 1) = cv(n1, 2)
            pts(3, 1) = cv(n1, 3)
            pts(1, 2) = cv(n2, 1)
            pts(2, 2) = cv(n2, 2)
            pts(3, 2) = cv(n2, 3)
            pts(1, 3) = cv(n3, 1)
            pts(2, 3) = cv(n3, 2)
            pts(3, 3) = cv(n3, 3)
            nv = 3
            if(iprcav >= 10) then
                write(lvpri, '(a, 3i6)') "Vertices indices: ", n1, n2, n3
                write(lvpri, *) 'Vertices'' Coordinates'
                call output(pts, 1_regint_k, 3_regint_k, 1_regint_k, 3_regint_k, 3_regint_k, 3_regint_k, 1_regint_k, lvpri)
            end if
            do jj = 1, 3
                pp(jj) = d0
                pp1(jj) = d0
            end do
            ! For every tessera, find the free portion (not covered by any other
            ! tesserae) and calculate its surface with the Gauss-Bonnet Theorem.
            ! The representative point is defined as the average of vertices
            ! belonging to the free portion of the tessera. The representative
            ! point is stored in the pp variable, while pp1 contains the
            ! coordinates of the point on the inward pointing normal vector.
            ! The vertices of each tessera are stored in vert(mxts, 10, 3), the
            ! number of vertices of each tessera is in nvert(mxts), the centers
            ! of the circles on each edge are in centr(mxts, 10, 3).
            ! The spheres to which teh edges of the tesserae belong to are
            ! stored in intsph(mxts, 10)
            if (iprcav > 15) then
                write(lvpri, '(a, i4)') "Before Tessera n.", nn + 1
                do iver = 1, 3
                    write(lvpri, '(a, i4, 3f15.9)') "XYZ vert n.", iver, (pts(jcor, iver), jcor=1, 3)
                end do
            end if
            call tessera(nsfe, nv, pts, ccc, pp, pp1, area, intsph, numts)
            if (iprcav >= 10) then
                write(lvpri, *) 'After tessera subroutine'
                write(lvpri, '(a, i6)') "Number of vertices: ", nv
                write(lvpri,*) 'Vertices'' Coordinates'
                call output(pts, 1_regint_k, 3_regint_k, 1_regint_k, nv, 3_regint_k, nv, 1_regint_k, lvpri)
            end if
            if(abs(area) <= 1.0d-14) then
                if (iprpcm >= 10) write(lvpri, '(a, i8)') "Zero area in tessera ", nn + 1
                go to 310
            end if
            if (nn >= mxts) then
                write(lvpri, '(a)') 'Too many tesserae in polyhedra!'
                write(lvpri, '(a, i4, a, i4)') "nn = ", nn, "  mxts = ", mxts
                pedra_error_code = 2
                stop
            end if
            nn = nn + 1
            xtscor(nn) = pp(1)
            ytscor(nn) = pp(2)
            ztscor(nn) = pp(3)
            xval(nn) = pp1(1)
            yval(nn) = pp1(2)
            zval(nn) = pp1(3)
            as(nn) = area
            isphe(nn) = nsfe
            nvert(nn) = nv
            do iv = 1, nv
                do jj = 1, 3
                    vert(nn, iv, jj) = pts(jj, iv)
                    centr(nn, iv, jj) = ccc(jj, iv)
                end do
            end do
            do iv = 1, nv
                intsph(nn, iv) = intsph(numts, iv)
            end do
            if(iprcav >= 10) then
                write(lvpri, '(a, i6)') "Adding tessera data to the list. nn = ", nn
                write(lvpri, '(a, 3f15.11)') "Center: ", xtscor(nn), ytscor(nn), ztscor(nn)
                write(lvpri, '(a, 3f15.11)') "Boh??? ", xval(nn), yval(nn), zval(nn)
                write(lvpri, '(a, f10.8, a, i4, a, i3)') "Area = ", as(nn), " Isphe = ", isphe(nn), " nvert = ", nvert(nn)
            end if
        310 end do
    end do
    nts = nn
    ! Check if two tesserae are too close
    test = 0.02d0
    test2 = test * test
    do i = 1, nts-1
        if(abs(as(i)) <= 1.0d-14) cycle
        xi = xtscor(i)
        yi = ytscor(i)
        zi = ztscor(i)
        ii = i + 1
        do j = ii , nts
            if(isphe(i) == isphe(j)) cycle
            if(abs(as(j)) <=  1.0d-14) cycle
            xj = xtscor(j)
            yj = ytscor(j)
            zj = ztscor(j)
            rij = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
            if(rij > test2) cycle
            write(lvpri, '(a, i8, a, i8, a, f8.6, a, f8.6, a)') " WARNING: The distance between center of tessera",  &
                            i, " and ", j, " is ", rij, ", less than ", test2, " AU"
        ! Check  if(as(i).lt.as(j)) as(i) = d0
        ! Check  if(as(i).ge.as(j)) as(j) = d0
        end do
    end do

! Replication of geometrical parameters according to simmetry
    call repcav(vert, centr, permutation_table, numts)

! Prepare data for geomview
    call ordpcm(nts, xtscor, ytscor, ztscor, as)

! Calculate cavity volume using Gauss Theorem:
!       V = sum_i {area(i) * center(i) * normal(i)} / 3
! the sum runs on the tesserae.
    vol  = 0.0d0
    do its = 1, nts
        nsfe = isphe(its)
        ! Find the unit normal vector
        xn = (xtscor(its) - xe(nsfe)) / re(nsfe)
        yn = (ytscor(its) - ye(nsfe)) / re(nsfe)
        zn = (ztscor(its) - ze(nsfe)) / re(nsfe)
        ! Calculate scalar product
        prod = xtscor(its)*xn + ytscor(its)*yn + ztscor(its)*zn
        ! And finally the volume
        vol = vol + as(its) * prod / 3.d0
    end do
! Print out some cavity data
    allocate(ssfe(mxsp))
    ssfe = 0.0_dp
    do i = 1, nts
        k = isphe(i)
        ssfe(k) = ssfe(k) + as(i)
    end do

     untesselated = 0
    if (some) then
        write(lvpri, '(/a)') "Final list of spheres"
        write(lvpri, '(a)') "Sphere             Center  (X,Y,Z) (AU)              Radius (AU)      Exposed surface (AU^2)"
        stot = 0.0_dp
        do i = 1, nesf
          if (abs(ssfe(i)) <= 1.0d-14) then
            untesselated = untesselated + 1
          else
            write(lvpri, '(i4, 4f15.9, f15.9)') i, xe(i), ye(i), ze(i), re(i), ssfe(i)
          end if
          stot = stot + ssfe(i)
        end do
        write(lvpri, '(/a, i8)') "Initial number of spheres = ", nesfp
        write(lvpri, '(a, i8)')  "Total number of spheres   = ", nesf
        write(lvpri, '(a, i8, a)') 'GePol added ', nesf-nesfp, ' spheres'
        write(lvpri, '(a, i8)')  ' Untesselated (internally buried) spheres = ', untesselated
        write(lvpri, '(a, i8)') " Total number of tesserae = ", nts
        write(lvpri, '(a, i8)') " Number of irreducible tesserae = ", ntsirr
        write(lvpri, '(a, f20.14, a)') " Cavity surface = ", stot, " AU^2"
        write(lvpri, '(a, f20.14, a/)') " Cavity volume  = ", vol, " AU^3"
    end if
    deallocate(ssfe)

!     ----- set up for possible gradient calculation -----
!           DONE ONLY IF WE HAVE SPHERES ON ATOMS
! RDR 191014 Deactivate call to derivative code
!   if(icesph /= 1) then
!
!       if(nesfp > nucdep) then
!           write(lvpri, '(a)') "PEDRA: confusion about the sphere count."
!           write(lvpri, '(a, i6, a, i6)') "nesfp = ", nesfp, "natm = ", nucdep
!           pedra_error_code = 3
!           stop
!       end if
!
!       call dzero(dercen, mxsp*mxcent*3*3)
!       call dzero(derrad, mxsp*mxcent*3)
!       do nsfe = 1, nucdep
!           natsph = 0
!           nsfer = nsfe
!           if (icesph == 2)then
!               natsph = 1
!               do jj = 1, nesfp
!                   if (ina(jj) == nsfe)then
!                       natsph = 0
!                       nsfer = jj
!                   end if
!               end do
!           end if
!           if(natsph == 0) THEN
!               do icoord = 1, 3
!                   call cavder(nsfe, nsfer, icoord, intsph, newsph)
!               end do
!           end if
!       end do
!   end if

    if (iprsol > 5) then
        write(lvpri, '(a)') " ***  Partition of the surface  ***"
        write(lvpri, '(a)') " Tessera  Sphere   Area   x y z tessera center  x y z normal vector"
        write(lvpri, '(2i4, 7f12.7)') (i,isphe(i),as(i),xtscor(i),ytscor(i), &
        ztscor(i),xtscor(nts+i),ytscor(nts+i), &
        ztscor(nts+i), i=1,nts)
    end if

! Check if the cavity is made up of disjoint spheres
    call cavspl(icav1, icav2, ncav1, ncav2, some)
! The dispersion calculation is allowed only in the case of single cavity.
    idisp = 0
    if(ncav2 /= 0) idisp=0

    idisrep = 0
    iretcav = 0
    if (idisp == 2) idisp = 1

    deallocate(permutation_table)

    end subroutine polyhedra

    subroutine sphper(nesf, nesf0, numsph, permutation_table, newsph, xe, ye, ze, re)
!
! Create permutation table for the spheres and complete the added
! spheres by symmetry.
! The permutation table permutation_table has dimension numsph*8, since 8 is the number of
! operations in D2h the biggest Abelian group.
! It reports which sphere is transformed into which under the operation
! considered.
! The ordering of operations is the DALTON ordering: E, Oyz, Oxz, C2z, Oxy, C2y,
! C2x, i
! The first column does not report how the sphere transforms under the identity,
! as it is trivially known. The first column does instead report if the sphere
! is Unique-For-Symmetry or not and hence needs to be tesselated or not.
! A class is composed by the spheres that are equivalent (i.e. overlapping
! exactly) under a given operation of symmetry.
! See also C. S. Pomelli, J. Tomasi and R. Cammi, J. Comput. Chem. 22, 1262 (2001)
!       Operations  1  2  3  4  5  6  7  8
!                  ------------------------
!       Spheres
!
!
    use pedra_symmetry, only: get_pt

#include "pcm_mxcent.inc"

    integer(kind=regint_k) :: nesf, nesf0, numsph
    real(kind=dp) :: xe(*), ye(*), ze(*), re(*), v1(3), v2(3)
    integer(kind=regint_k) :: permutation_table(nesf, *), newsph(numsph, 2)

    real(kind=dp) :: diff1, r1
    integer(kind=regint_k) :: i, j, k, l, n1, n2, nesf1

    nesf1 = nesf

    ! Loop over spheres
    do i = 1, nesf
        v1(1) = xe(i)
        v1(2) = ye(i)
        v1(3) = ze(i)
        r1    = re(i)
        ! Loop over symmetry operations
        do j = 1, group%maxrep
            ! Loop over spheres
            do k = 1, nesf
                do l = 1, 3
                    v2(l) = get_pt(iand(group%isymax(l, 1), j)) * v1(l)
                end do
                diff1 = sqrt((xe(k) - v2(1))**2 + (ye(k) - v2(2))**2 + (ze(k) - v2(3))**2 + (re(k) - r1)**2)
                if ((diff1 < 1.0D-3)) then
                    permutation_table(i, j+1) = k
                    go to 10
                end if
            end do
            ! Check if we need to complete the creation of additional sphere
            ! or the initial spheres are ill-defined.
            if (i < nesf0) then
                write(lvpri, '(a)') "Cavity is not consistent with symmetry"
                pedra_error_code = 4
                stop
            else
                nesf1 = nesf1 + 1
                permutation_table(i, j+1) = nesf1
                xe(nesf1) = v2(1)
                ye(nesf1) = v2(2)
                ze(nesf1) = v2(3)
                re(nesf1) = r1
                n1 = permutation_table(abs(newsph(i, 1)), j+1)
                n2 = permutation_table(abs(newsph(i, 2)), j+1)
                if (newsph(i, 1) < 0) then
                    n1 = -n1
                    n2 = -n2
                end if
                newsph(nesf1, 1) = n1
                newsph(nesf1, 2) = n2
                write(lvpri,*) "New sphere added in sphper ", nesf1, n1, n2
            end if
            10 continue
        end do
    end do

! If needed complete the permutation table for the new spheres
    do i = nesf + 1, nesf1
        v1(1) = xe(i)
        v1(2) = ye(i)
        v1(3) = ze(i)
        r1    = re(i)

        ! Loop over symmetry operations
        do j = 1, group%maxrep
            do k = 1, nesf1
                do l = 1, 3
                    v2(l) = get_pt(iand(group%isymax(l, 1), j)) * v1(l)
                end do
                diff1 = sqrt((xe(k) - v2(1))**2 + (ye(k) - v2(2))**2 + (ze(k) - v2(3))**2 + (re(k) - r1)**2)
                if ((diff1 < 1.0D-3)) then
                    permutation_table(i, j+1) = k
                    go to 20
                end if
            end do
            write(lvpri, '(a)') "Additional spheres not consistent with symmetry"
            pedra_error_code = 5
            stop
            20 continue
        end do
    end do

    nesf = nesf1
! Now define which spheres are to be tesselated: it will be the first
! in each class of spheres equivalent by symmetry (first column of permutation_table array)
    do i = 1, nesf
        permutation_table(i, 1) = 1
        do k = 1, group%maxrep
            if (permutation_table(i, k+1) < i) then
                permutation_table(i, 1) = 0
                cycle
            end if
        end do
    end do

    end subroutine sphper

    subroutine polygen(ipflag, tsare, itsnum, xen, yen, zen, ren, itseff, cv,  &
                       jtr, permutation_table, nesf, nsfe, numts, numver, rotcav)
!
! Polygen: a program to generate spherical polyhedra with triangular faces.
! An equilateral division algorithm is used.
! Polygen can generate an optimal tesselation based on used input.
! Using the ipflag integer(kind=regint_k), the user can request:
!    - a spherical polyhedron with an optimal number of tesserae (ipflag = 0);
!    - a spherical polyhedron with an optimal average tesserae area (ipflag = 1)
! The polyhedron whose vertices are to be projected on the sphere and the
! subdivision frequency (nf) for the subsequent equilateral division procedure are
! chosen based on ipflag.
!
    use pedra_symmetry, only: get_pt

#include "pcm_pcmdef.inc"
#include "pcm_mxcent.inc"

    integer(kind=regint_k) :: ipflag, itsnum, itseff, nsfe, numts, numver
    integer(kind=regint_k), intent(in) :: nesf
    integer(kind=regint_k), intent(in) :: permutation_table(nesf, *)
    integer(kind=regint_k) :: jtr(numts, *)
    real(kind=dp) :: cv(numver, 3)
    real(kind=dp) :: tsare, xen, yen, zen, ren

    real(kind=dp) :: v1(3), v2(3), v3(3), rotcav(3, 3)
    integer(kind=regint_k) :: itrvo(60,3), itreo(60,3), iedo(90,2)
    integer(kind=regint_k) :: oldtr(100,100), ednew(90,100), trnew(60,100,100)

    real(kind=dp), parameter :: d0 = 0.0d0
    real(kind=dp), parameter :: d1 = 1.0d0
    real(kind=dp), parameter :: pi = acos(-1.0d0)
    real(kind=dp) :: alpha, beta, cos1, cos2, costheta, dl, dm, dn, dnf, dnorm
    real(kind=dp) :: sintheta, theta
    integer(kind=regint_k) :: i, ii, isymop, j, jj, jsymop, k, l, m, n, ne0, nf
    integer(kind=regint_k) :: noppt, nt, nt0, ntpt, ntra, nv, nvpt

    itrvo = 0
    itreo = 0
    iedo  = 0

    ipflag = 1
    if (ipflag == 1) then
        itsnum = int(4.0d0 * pi * ren**2 / tsare + 0.5d0 )
    end if

    if (itsnum > mxts) then
        write(lvpri, '(a, 2i6)') "Requested number of tesserae in polyhedron exceeds maximum number of tesserae: ", itsnum, mxts
        pedra_error_code = 6
        stop
    end if

! Build the initial tesselation. This depends on the point group and (for C1) on
! the number of tesserae requested.
    nt0 = 1
    nv = 3
    ne0 = 3
    itrvo(1, 1) = 1
    itrvo(1, 2) = 2
    itrvo(1, 3) = 3
    iedo(1, 1) = 1
    iedo(1, 2) = 2
    iedo(2, 1) = 2
    iedo(2, 2) = 3
    iedo(3, 1) = 1
    iedo(3, 2) = 3
    itreo(1, 1) = 1
    itreo(1, 2) = 2
    itreo(1, 3) = 3

    cv(1, 1) = d1
    cv(1, 2) = d0
    cv(1, 3) = d0

    cv(2, 1) = d0
    cv(2, 2) = d1
    cv(2, 3) = d0

    cv(3, 1) = d0
    cv(3, 2) = d0
    cv(3, 3) = d1

! Calculate the subdivision frequency.
    nf = int(sqrt(d1 * itsnum / ( d1 * nt0  * 8 )) + 0.5d0)
    if (nf <= 1) nf = 2
    if (nf <= 2) then
            write(lvpri, '(a/)') "** WARNING  ** A very poor tesselation has " &
            //"been chosen. It is valuable almost only for testing."
    end if
    if (nf >= 8) then
        write(lvpri, '(a)') "** WARNING ** A very fine tesselation has been " &
        //"chosen, it will probably produce a HUGE amount of tesserae."
    end if
    dnf = dble(nf)
! eck WRITE(LVPRI,*) 'dopo polydata',NT0,NE0,NV,NF

    itseff = nt0 * nf**2

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
        end do
        costheta= &
        ( v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)) / &
        (sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2) * &
        sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2))
        theta=acos(costheta)
        sintheta=sin(theta)
        DO l=1,NF-1
            DL = dble(L)
        ! F          m=NF-l
            cos1=cos(theta*Dl/DNF)
            cos2=cos(theta*(DNF-Dl)/DNF)
            alpha=(cos1-costheta*cos2)/sintheta**2
            beta=(cos2-costheta*cos1)/sintheta**2
            DM = 0.0D0
            ALPHA = dble(NF-L)
            BETA  = dble(L)
            DO k=1,3
                v3(k)=alpha*v1(k)+beta*v2(k)
            end do
            dnorm = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
            DO k=1,3
                v3(k) = v3(k) /dnorm
            end do
            DO K = 1,3
                CV(NVPT,K)=V3(K)
            end do
            ednew(j,l+1)=NVPT
            NVPT=NVPT+1
        end do
    end do


!  -nuovi vertici non posti lungo i vecchi spigoli
!  -a partire dai vertici in ednew secondo regole analoghe alle
!  precedenti, vengono memorizzati a seconda del triangolo in trnew
!  trnew(triangolo,fila,n ordine)=nvertice

! f      DO J=1,NT0
    j=1
    ii=1
    jj=3
    DO L =3,NF
        DL = dble(L)
        DO N=1,l-2
            DM = dble(NF - L + 1)
            DL = dble(L - 1 - N)
            DN = dble(N)
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
        end do
    end do
! f       end do
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
        end do

    !  -3 nuovi vertici non lungo i vecchi spigoli

        DO l=3,NF
            DO m=2,l-1
                oldtr(l,m)=trnew(n,l,m)
            end do
        end do

    !  -ora si creano i nuovi triangoli

        DO i=1,NF
            DO j=1,i
                JTR(NTPT,1)=oldtr(i,j)
                JTR(NTPT,2)=oldtr(i+1,j)
                JTR(NTPT,3)=oldtr(i+1,j+1)
                NTPT=NTPT+1
            end do
        end do
        DO i=2,NF
            DO j=1,i-1
                JTR(NTPT,1)=oldtr(i,j)
                JTR(NTPT,2)=oldtr(i,j+1)
                JTR(NTPT,3)=oldtr(i+1,j+1)
                NTPT=NTPT+1
            end do
        end do
    2310 end do
    NV=NVPT-1
    NT=NTPT-1

! Generate the right irreducible part of the tesselation for the current sphere
! under the selected point group.
    call prerep(nv, nt, itseff, cv, jtr, numver, numts)
    do i = 1, nv
        v1 = 0.0d0
        do j = 1,3
            do k = 1,3
                v1(j) = v1(j) + rotcav(k, j) * cv(i, k)
            end do
        end do
        cv(i, :) = v1
    end do

! Replication with local operators, for those cases where the local operator
! transforms the given sphere into another sphere.
    noppt=0
    do isymop = 1, group%maxrep
        ! Get the index of the transformed sphere ntra when applying the
        ! operator with index isymop + 1 to the sphere with index nsfe
        ntra = permutation_table(nsfe, isymop + 1)
        ! The if statment is entered only if the sphere which is equivalent
        ! to the one we are tesselating now under the operation isymop, is
        ! not itself, i.e. if the symmetry operator doesn't bring the current
        ! sphere into itself.
        if (ntra /= nsfe) then
            ! If the replication of this part of the cavity has already been
            ! done we start checking a new simmetry operation
            do jsymop = 1 , isymop-1
                if (permutation_table(nsfe, jsymop + 1) == ntra) go to 100
            end do
            noppt = noppt + 1
            ! Replication of vertices
            do i = 1, nv
                ii = i + noppt * nv
                do k = 1, 3
                    cv(ii, k) = get_pt(iand(group%isymax(k, 1), isymop)) * cv(i,k)
                end do
            end do
            ! Replication of topology
            do i = 1, nt
                ii = i + noppt * nt
                jj = noppt * nv
                do k = 1, 3
                    jtr(ii, k) = jtr(i, k) + jj
                end do
            end do
            100 continue
        end if
    end do
    ! Update indices
    nv = nv * (noppt + 1)
    itseff = itseff * (noppt + 1)

    do i = 1, nv
        cv(i,1) = cv(i,1) * ren + xen
        cv(i,2) = cv(i,2) * ren + yen
        cv(i,3) = cv(i,3) * ren + zen
    end do

    end subroutine polygen

    subroutine repcav(vert, centr, permutation_table, numts)
!
! Reproduce the irreducibile part of the cavity
!

    use pedra_symmetry, only: get_pt

#include "pcm_pcmdef.inc"
#include "pcm_mxcent.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: numts
    real(kind=dp) :: vert(numts, 10, 3), centr(numts, 10, 3)
    integer(kind=regint_k) :: permutation_table(nesf, *)

    integer(kind=regint_k) :: i, ii, ii2, isymop, k, l

    ntsirr = nts
    nts = nts * (group%maxrep + 1)

    if (nts > mxts) then
        write(lvpri, '(a)') "Number of tesserae exceeds maximum."
        pedra_error_code = 7
        stop
    end if
    ! Loop over symmetry operations. The identity is excluded from the loop.
    do isymop = 1, group%maxrep
        do i = 1, ntsirr
            ii        = i + ntsirr * isymop
            as(ii)    = as(i)
            nvert(ii) = nvert(i)
            isphe(ii) = permutation_table(isphe(i), isymop+1)
            do k = 1, nvert(i)
                do l = 1, 3
                  vert(ii, k, l)  = get_pt(iand(group%isymax(l, 1), isymop)) * vert(i, k, l)
                  centr(ii, k, l) = get_pt(iand(group%isymax(l, 1), isymop)) * centr(i, k, l)
                end do
            end do
        end do
    end do

! Moving the normal points in a safe location where they will not be overwritten
    do i = 1, ntsirr
        ii          = i + ntsirr
        ii2         = i + (group%maxrep + 1) * ntsirr
        xtscor(ii2) = xtscor(ii)
        ytscor(ii2) = ytscor(ii)
        ztscor(ii2) = ztscor(ii)
    end do

    do isymop = 1, group%maxrep
        do i = 1, ntsirr
            ii         = i + ntsirr * isymop
            xtscor(ii) = get_pt(iand(group%isymax(1, 1), isymop)) * xtscor(i)
            ytscor(ii) = get_pt(iand(group%isymax(2, 1), isymop)) * ytscor(i)
            ztscor(ii) = get_pt(iand(group%isymax(3, 1), isymop)) * ztscor(i)
        !            ii=i+nts*(isymop+group%maxrep+1)
        !            ii2=i+nts
        !            xtscor(ii)=get_pt(iand(group%isymax(1,1),isymop)) * xtscor(ii2)
        !            ytscor(ii)=get_pt(iand(group%isymax(2,1),isymop)) * ytscor(ii2)
        !            ztscor(ii)=get_pt(iand(group%isymax(3,1),isymop)) * ztscor(ii2)
        end do
    end do

    end subroutine repcav

    subroutine tessera(ns, nv, pts, ccc, pp, pp1, area, intsph, numts)

    use pedra_utils, only: around
    use pedra_print, only: output

#include "pcm_pcmdef.inc"
#include "pcm_mxcent.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: ns, nv, numts
    real(kind=dp) :: area
    real(kind=dp) :: pts(3,10), ccc(3,10), pp(3), pp1(3)
    integer(kind=regint_k) :: intsph(numts,10)
    real(kind=dp) :: p1(3), p2(3), p3(3), p4(3), point(3)
    real(kind=dp) :: pscr(3,10),cccp(3,10),pointl(3,10)
    integer(kind=regint_k) :: ind(10), ltyp(10), intscr(10), ntrhso(10)
    logical :: lan

    real(kind=dp) :: dcheck, de2, delr, delr2, diffdr, dist, dist1, dist2
    real(kind=dp) :: dnorm, rc, rc2, tol
    real(kind=dp) :: x1, x2, y1, y2, z1, z2
    integer(kind=regint_k) :: i, j, ic, icop, icut, idx, idx2, ii, intcas, iprcav
    integer(kind=regint_k) :: iv1, iv2, ivnew, ivold, jj, k, l, n, nsfe1, nvleft, nvnegl


!     Coord. del centro che sottende l`arco tra i vertici
!     n e n+1 (per i primi tre vertici e' sicuramente il centro della
!     sfera) e sfera alla cui intersezione con NS appartiene l'arco (se
!     appartiene alla sfera originaria INTSPH(numts,N)=NS)

    iprcav = 0
    lan = .false.
    area = 0.0d+00
    do j = 1, 3
        ccc(1, j) = xe(ns)
        ccc(2, j) = ye(ns)
        ccc(3, j) = ze(ns)
    end do
    if (iprcav >= 10) then
        call around('Input data in tessera', LVPRI)
        write(lvpri,1000) ns, nv, numts
        write(lvpri, *) '=======pts======='
        call output(pts, 1_regint_k, 3_regint_k, 1_regint_k, 3_regint_k, 3_regint_k, 3_regint_k, 1_regint_k, lvpri)
        write(lvpri, *) '=======ccc======='
        call output(ccc, 1_regint_k, 3_regint_k, 1_regint_k, 3_regint_k, 3_regint_k, 3_regint_k, 1_regint_k, lvpri)
    end if

!     INTSPH viene riferito alla tessera -numts-, e in seguito riceve il
!     numero corretto.

    DO N = 1, 3
        INTSPH(NUMTS,N) = NS
    end do

!     Loop sulle altre sfere

    DO 150 NSFE1=1,NESF
        IF(NSFE1 == NS) go to 150
        IF(IPRCAV >= 10) THEN
            WRITE(LVPRI,1005) NSFE1
            WRITE(LVPRI,1010) XE(NSFE1),YE(NSFE1),ZE(NSFE1),RE(NSFE1)
        end if

    !     Memorizza i vertici e i centri che sottendono gli archi

        DO J =1, NV
            INTSCR(J) = INTSPH(NUMTS,J)
            DO I = 1,3
                PSCR(I,J) = PTS(I,J)
                CCCP(I,J) = CCC(I,J)
            end do
        end do
        IF(IPRCAV >= 10) THEN
            CALL AROUND('ACTUAL TESSERA_ STATUS', lvpri)
            WRITE(LVPRI,1000) NS,NV,NUMTS
            WRITE(LVPRI,*) '=======PSCR======='
            CALL OUTPUT(PSCR,1_regint_k,3_regint_k,1_regint_k,NV,3_regint_k,NV,1_regint_k,LVPRI)
            WRITE(LVPRI,*) '=======CCCP======='
            CALL OUTPUT(CCCP,1_regint_k,3_regint_k,1_regint_k,NV,3_regint_k,NV,1_regint_k,LVPRI)
        end if


        ICOP = 0
        DO J =1, 10
            IND(J) = 0
            LTYP(J) = 0
        end do

    !     Loop sui vertici della tessera considerata

        DO 100 I=1,NV
            DELR2=(PTS(1,I)-XE(NSFE1))**2+(PTS(2,I)-YE(NSFE1))**2+ &
            (PTS(3,I)-ZE(NSFE1))**2
            DELR=SQRT(DELR2)
            IF (IPRPCM >= 10) THEN
                WRITE(LVPRI,1015) DELR,RE(NSFE1)
            end if
            DIFFDR = DELR - RE(NSFE1)
            IF(DIFFDR < 0.0D0) THEN
                IND(I) = 1
                ICOP = ICOP+1
            end if
        100 end do
    !     Se la tessera e' completamente coperta, la trascura
        IF(ICOP == NV) THEN
            IF(IPRCAV >= 10) WRITE(LVPRI,1020) NSFE1
            RETURN
        !           ******
        end if


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
                        end do
                        go to 160
                    end if
                end do
            end if
            160 CONTINUE
            IF(IPRCAV >= 10) THEN
                WRITE(LVPRI,1025) L,TYPLAB(LTYP(L))
            end if
        end do

    !     Se la tessera e' spezzata in due o piu' tronconi, la trascura

        ICUT = 0
        DO L = 1, NV
            IF(LTYP(L) == 1 .OR. LTYP(L) == 2) ICUT = ICUT + 1
            IF(LTYP(L) == 3) ICUT = ICUT + 2
        end do
        ICUT = ICUT / 2
        IF(ICUT > 1) THEN
            IF(IPRPCM >= 10) WRITE(LVPRI,*) 'Tessera cut in pieces and removed.'
            RETURN
        end if

    !     Creazione dei nuovi vertici e lati della tessera
    !     Loop sui lati

        N = 1
        IF(IPRCAV >= 10) WRITE(LVPRI,*) &
        'NOW CREATING NEW VERTICES AND EDGES IF NEED BE....'
        DO 300 L = 1, NV
        !     Se il lato L e' coperto:
            IF(LTYP(L) == 0) go to 300
            IV1 = L
            IV2 = L+1
            IF(L == NV) IV2 = 1
            IF (IPRCAV >= 10) THEN
                WRITE(LVPRI,1030) L,TYPLAB(LTYP(L)),LTYP(L)
                WRITE(LVPRI,1035) IV1,IV2
            end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     Se il lato L e' tagliato (con il I vertice scoperto):
            IF(LTYP(L) == 1) THEN
            ! f store first point in the final set
                DO JJ = 1, 3
                    PTS(JJ,N) = PSCR(JJ,IV1)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                end do
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
                end do
                INTCAS=0
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1040) INTCAS,NSFE1
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                end if
                CALL INTER(P1,P2,P3,P4,NSFE1,INTCAS)
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1050)
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                end if
                DIST1 = SQRT((P4(1)-P1(1))**2 + (P4(2)-P1(2))**2 + &
                (P4(3)-P1(3))**2)
                DIST2 = SQRT((P4(1)-P2(1))**2 + (P4(2)-P2(2))**2 + &
                (P4(3)-P2(3))**2)
            !     Aggiorna i vertici della tessera e il centro dell'arco
                DO JJ = 1,3
                    PTS(JJ,N) = P4(JJ)
                end do

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
            end if

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
                end do
                INTCAS=1
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1040) INTCAS,NSFE1
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                end if
                CALL INTER(P1,P2,P3,P4,NSFE1,INTCAS)
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1050)
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                end if
            !     Aggiorna i vertici della tessera e il centro dell'arco
                DO JJ = 1,3
                    PTS(JJ,N) = P4(JJ)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                end do
                INTSPH(NUMTS,N) = INTSCR(IV1)
                N = N+1
            end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     Se il lato e' intersecato due volte:
            IF(LTYP(L) == 3) THEN
                DO JJ = 1, 3
                    PTS(JJ,N) = PSCR(JJ,IV1)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                end do
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
                end do
                INTCAS=0
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1040) INTCAS,NSFE1
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                end if
                CALL INTER(P1,P2,P3,P4,NSFE1,INTCAS)
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1050)
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                end if
            !     Aggiorna i vertici della tessera e il centro dell'arco
                DO JJ = 1,3
                    PTS(JJ,N) = P4(JJ)
                end do

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
                end do
                INTCAS=1
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1040) INTCAS,NSFE1
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                end if
                CALL INTER(P1,P2,P3,P4,NSFE1,INTCAS)
                IF(IPRCAV >= 10) THEN
                    WRITE(LVPRI,1050)
                    WRITE(LVPRI,1045) 1,(P1(I),I=1,3)
                    WRITE(LVPRI,1045) 2,(P2(I),I=1,3)
                    WRITE(LVPRI,1045) 3,(P3(I),I=1,3)
                    WRITE(LVPRI,1045) 4,(P4(I),I=1,3)
                end if
                DIST1 = SQRT((P4(1)-P1(1))**2 + (P4(2)-P1(2))**2 + &
                (P4(3)-P1(3))**2)
                DIST2 = SQRT((P4(1)-P2(1))**2 + (P4(2)-P2(2))**2 + &
                (P4(3)-P2(3))**2)
            !     Aggiorna il vertice e il centro dell'arco
                DO JJ = 1,3
                    PTS(JJ,N) = P4(JJ)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                end do
                INTSPH(NUMTS,N) = INTSCR(IV1)
                N = N + 1
            end if

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !     Se il lato e' scoperto:
            IF(LTYP(L) == 4) THEN
                DO JJ = 1, 3
                    PTS(JJ,N) = PSCR(JJ,IV1)
                    CCC(JJ,N) = CCCP(JJ,IV1)
                end do
                INTSPH(NUMTS,N) = INTSCR(IV1)
                N = N+1
            end if
        300 end do
        NV = N - 1
        IF(IPRCAV >= 10) THEN
            CALL AROUND('AFTER INTER_SECTION TESSERA_ STATUS', lvpri)
            WRITE(LVPRI,1000) NS,NV,NUMTS
            WRITE(LVPRI,*) '=======PTS======='
            CALL OUTPUT(PTS,1_regint_k,3_regint_k,1_regint_k,NV,3_regint_k,NV,1_regint_k,LVPRI)
            WRITE(LVPRI,*) '=======CCC======='
            CALL OUTPUT(CCC,1_regint_k,3_regint_k,1_regint_k,NV,3_regint_k,NV,1_regint_k,LVPRI)
            DO IDX=1,NV
                IDX2=IDX+1
                IF (IDX2 > NV) IDX2=1
                DIST1=0.0D0
                DIST2=0.0D0
                DO IC=1,3
                    DIST1 = DIST1 + (PTS(IC,IDX ) - CCC(IC,IDX))**2
                    DIST2 = DIST2 + (PTS(IC,IDX2) - CCC(IC,IDX))**2
                end do
                DCHECK=DIST1-DIST2
                WRITE(LVPRI,1055) IDX,DCHECK
            end do
        end if
    !     Controlla che il numero di vertici creati non sia eccessivo
        IF(NV > 10) THEN
            WRITE(LVPRI,*)'TOO MANY VERTICES IN TESSERA_: BYE BYE...'
            pedra_error_code = 8
            STOP
        end if
    150 end do
    IF(IPRCAV >= 10) THEN
        CALL AROUND('FINAL TESSERA_ STATUS', lvpri)
        WRITE(LVPRI,1000) NS,NV,NUMTS
        WRITE(LVPRI,*) '=======PSCR======='
        CALL OUTPUT(PSCR,1_regint_k,3_regint_k,1_regint_k,NV,3_regint_k,NV,1_regint_k,LVPRI)
        WRITE(LVPRI,*) '=======CCCP======='
        CALL OUTPUT(CCCP,1_regint_k,3_regint_k,1_regint_k,NV,3_regint_k,NV,1_regint_k,LVPRI)
    end if

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
        end if
    end do
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
                end do
            end if
        end do
        IF (IVNEW /= NVLEFT) THEN
            WRITE(LVPRI,*) 'SIRCAV: BADLY MADE ALGORITHM!'
            pedra_error_code = 9
            STOP
        end if
        NV = NVLEFT
    end if

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

#include "pcm_mxcent.inc"
#include "pcm_pcmdef.inc"
#include "pcm_pcm.inc"

    real(kind=dp) :: p1(3), p2(3), p3(3), p4(3)
    integer(kind=regint_k) :: ns, i
    logical :: lalow, lblow

    real(kind=dp) :: alphat, delta, diff, diff2, diff2a, diff2b
    real(kind=dp) :: diffa, diffb, dnorm, p1p3, p2p3, r, r2, tol
    integer(kind=regint_k) :: j, jj, m
    integer(kind=regint_k) :: iprcav = 0

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
    end if
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
        end do
        go to 100
    end if
    IF(LBLOW .AND. .NOT. LALOW) THEN
        IF(IPRCAV >= 10) WRITE (LVPRI,*) &
        'INTER_: TAKEN SECOND POINT'
        DO J=1,3
            P4(J)=P2(J)
        end do
        go to 100
    end if
    IF(LALOW .AND. LBLOW) THEN
        IF(DIFFA <= DIFFB) THEN
            IF(IPRCAV >= 10) WRITE (LVPRI,*) &
            'INTER_: TAKEN FIRST POINT'
            DO J=1,3
                P4(J)=P1(J)
            end do
        ELSE
            IF(IPRCAV >= 10) WRITE (LVPRI,*) &
            'INTER_: TAKEN SECOND POINT'
            DO J=1,3
                P4(J)=P2(J)
            end do
        end if
        go to 100
    end if

!     Start iterations

    M = 1
    10 CONTINUE
    ALPHAT = ALPHAT + DELTA
    DNORM = 0.0D+00
    DO JJ = 1,3
        P4(JJ)=P1(JJ)+ALPHAT*(P2(JJ)-P1(JJ))-P3(JJ)
        DNORM = DNORM + P4(JJ)**2
    end do
    DNORM = SQRT(DNORM)
    DO JJ = 1,3
        P4(JJ)= P4(JJ)*R/DNORM + P3(JJ)
    end do
    DIFF2=(P4(1)-XE(NS))**2 + (P4(2)-YE(NS))**2 + (P4(3)-ZE(NS))**2
    DIFF = SQRT(DIFF2) - RE(NS)
    IF(ABS(DIFF) < TOL) go to 100
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
        pedra_error_code = 10
        STOP
    end if
    IF (M > 300)THEN
        WRITE(LVPRI,*)'Too many iterations in INTER_! BYE BYE ...'
        WRITE(LVPRI,*)'P1',P1(1),P1(2),P1(3),DIFFA,LALOW
        WRITE(LVPRI,*)'P2',P2(1),P2(2),P2(3),DIFFB,LBLOW
        WRITE(LVPRI,*)'P3',P3(1),P3(2),P3(3)
        WRITE(LVPRI,*)'P4',P4(1),P4(2),P4(3),DIFF,I
        WRITE(LVPRI,*)'SPHERE',XE(NS),YE(NS),ZE(NS),RE(NS)
        WRITE(*,*) '8'
        STOP
    end if
    go to 10

! Final printing and return

    100 CONTINUE
    IF(IPRPCM >= 10) THEN
        WRITE(LVPRI,1005)
        WRITE(LVPRI,1010) 1,P1(1),P1(2),P1(3)
        WRITE(LVPRI,1010) 2,P2(1),P2(2),P2(3)
        WRITE(LVPRI,1010) 3,P3(1),P3(2),P3(3)
        WRITE(LVPRI,1010) 4,P4(1),P4(2),P4(3)
        WRITE(LVPRI,1015) XE(NS),YE(NS),ZE(NS),RE(NS)
    end if
    RETURN

    1000 FORMAT('INTER_, distance consistency check: ',F14.12)
    1005 FORMAT(/'Final result from INTER_!')
    1010 FORMAT('P',I1,' X = ',F12.9,' Y = ',F12.9,' Z = ',F12.9)
    1015 FORMAT('SPHERE',' X = ',F12.9,' Y = ',F12.9,' Z = ',F12.9, &
    ' R = ',F12.9)

    end subroutine inter

    subroutine gaubon(nv, ns, pts, ccc, pp, pp1, area, intsph, numts)

    use pedra_dblas, only: vector_product

#include "pcm_pcmdef.inc"
#include "pcm_mxcent.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: nv, ns, numts
    real(kind=dp) :: area
    real(kind=dp) :: pts(3, 10), ccc(3, 10), pp(3), pp1(3), beta(10)
    integer(kind=regint_k) :: intsph(numts, 10)
    real(kind=dp) :: p1(3),p2(3),p3(3),u1(3),u2(3),phin(10),weight(0:10)

    real(kind=dp), parameter :: d0 = 0.0d0
    real(kind=dp), parameter :: d1 = 1.0d0
    real(kind=dp), parameter :: pi = acos(-1.0d0)

    real(kind=dp) :: cosphin, costn, sum1, sum2, sumphi, dnorm1, dnorm2, tpi
    real(kind=dp) :: scal, dnorm, dnorm3
    real(kind=dp) :: x1, x2, y1, y2, z1, z2
    integer(kind=regint_k) :: i, jj, n, n0, n1, n2, nsfe1

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

    !         IF (NTRHSO(N) .EQ. 1) go to 100

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
        end if
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
        IF(abs(DNORM1) <= 1.0d-14) DNORM1 = 1.0D+00
        X2 = PTS(1,N) - XE(NS)
        Y2 = PTS(2,N) - YE(NS)
        Z2 = PTS(3,N) - ZE(NS)
        DNORM2 = SQRT(X2*X2 + Y2*Y2 + Z2*Z2)
        COSTN = (X1*X2+Y1*Y2+Z1*Z2)/(DNORM1*DNORM2)
        SUM1 = SUM1 + PHIN(N) * COSTN
    100 end do

    do i = 1, nv
        weight(i) = phin(i)
    end do
    weight(0) = weight(nv)
    do i = nv, 1, -1
        weight(i) = weight(i) + weight(i-1)
    end do

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
        end do
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
        N0 = MOD(NV+N0-1_regint_k,NV)
        IF(N0 == 0) N0 = NV
    ! f         write(lvpri,*) "N0 is", N0, NTRHSO(N0)
    ! f         IF(NTRHSO(N0).EQ.1) go to 220
        N2 = N
    ! f         NCONT = 0
    ! f 230     CONTINUE
        N2 = MOD(N2+1_regint_k,NV)
        IF(N2 == 0) N2 = NV
    ! f         write(lvpri,*) "N2 is", N2, NTRHSO(N1)
    ! f         IF(NTRHSO(N1 + NCONT).EQ.1) THEN
    ! f            NCONT = NCONT + 1
    ! f            go to 230
    ! F         end if

    !     Trova i vettori posizione rispetto ai centri corrispondenti
    !     e i versori tangenti

    !     Lato N0-N1:
        DO JJ = 1, 3
            P1(JJ) = PTS(JJ,N1) - CCC(JJ,N0)
            P2(JJ) = PTS(JJ,N0) - CCC(JJ,N0)
        end do

        CALL vector_product(P1,P2,P3,DNORM3)
        DO JJ = 1, 3
            P2(JJ) = P3(JJ)
        end do
        CALL vector_product(P1,P2,P3,DNORM3)
        DO JJ = 1, 3
            U1(JJ) = P3(JJ)/DNORM3
        end do

    !     Lato N1-N2:
        DO JJ = 1, 3
            P1(JJ) = PTS(JJ,N1) - CCC(JJ,N1)
            P2(JJ) = PTS(JJ,N2) - CCC(JJ,N1)
        end do

        CALL vector_product(P1,P2,P3,DNORM3)
        DO JJ = 1, 3
            P2(JJ) = P3(JJ)
        end do
        CALL vector_product(P1,P2,P3,DNORM3)
        DO JJ = 1, 3
            U2(JJ) = P3(JJ)/DNORM3
        end do

        BETA(N) = ACOS(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
    ! F         SUM2 = SUM2 + (PI - BETAN)
    200 end do
    do I=1,NV
    ! F         II = I - 1
    ! f         IF (II .EQ. 0) II = NV
    ! f         SUM2 = SUM2 - BETA(I) * (D1 - DP5 * DBLE(NTRHSO(I)+NTRHSO(II)))
    ! f     $               + PI * (D1 - DBLE(NTRHSO(I)))
        SUM2 = SUM2 - BETA(I)
    end do
!     Calcola l'area della tessera
    AREA = RE(NS)*RE(NS)*(DBLE(2-NV) * PI + SUM1 - SUM2)
!     Trova il punto rappresentativo (come media dei vertici)
    DO JJ = 1, 3
        PP(JJ) = 0.0D+00
    end do
    DO I = 1, NV
        PP(1) = PP(1) + (PTS(1,I)-XE(NS)) * WEIGHT(I)
        PP(2) = PP(2) + (PTS(2,I)-YE(NS)) * WEIGHT(I)
        PP(3) = PP(3) + (PTS(3,I)-ZE(NS)) * WEIGHT(I)
    end do
    DNORM = 0.0D+00
    DO JJ = 1, 3
        DNORM = DNORM + PP(JJ)*PP(JJ)
    end do
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
    end if

    END SUBROUTINE GAUBON

    subroutine cavspl(icav1, icav2, ncav1, ncav2, some)
!
! Check if the solute is contained into a single cavity or in
! more disjoint cavities. If yes, the order numbers of the spheres
! for the first and second cavity are stored into icav1 and icav2 respectively.
!
#include "pcm_mxcent.inc"
#include "pcm_pcmdef.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k) :: icav1(mxcent), icav2(mxcent)
    integer(kind=regint_k) :: ncav1, ncav2
    logical :: some

    integer(kind=regint_k) :: i, n1, n2, nn, icen, n
    real(kind=dp) :: r, rr, sum, x, y, z, xx, yy, zz

    DO I=1,MXCENT
        ICAV1(I)=0
        ICAV2(I)=0
    end do
    NCAV1=0
    NCAV2=0
    N1=1
    N2=1
    NN=1
    ICAV1(N1)=1
    N1=N1+1
    50 DO I=1,MXCENT
        ICAV2(I)=0
    end do
    N2=1
    N=ICAV1(NN)
    DO 100 ICEN=1,NESF
        IF(ICEN == N) go to 100
        DO I=1,NESF
            IF(ICEN == ICAV1(I)) go to 100
        end do
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
        end if
    100 end do
    NN=NN+1
    IF(ICAV1(NN) /= 0) go to 50
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
        end if
    end if
    RETURN

    200 FORMAT(/10X,'THE SOLUTE IS ENCLOSED IN ONE CAVITY')
    300 FORMAT(/10X,'THE SOLUTE IS ENCLOSED IN TWO DISTINCT CAVITIES'/ &
    10X,'OF',I3,3X,'E',I3,3X,'SPHERE(S),  RESPECTIVELY')
    400 FORMAT(/10X,'THE FIRST CAVITY IS FORMED BY SPHERE(S) :'/)
    500 FORMAT(/10X,'THE SECOND CAVITY IS FORMED BY SPHERE(S) :'/)

    end subroutine cavspl

    subroutine prerep(nv, nt, its, cv, jtr, nvert, numts)
!
! This subroutine identifies the symmetry elements needed for the replication of
! the cavity. The algorithm always starts from the D2h tesselation of a single
! sphere, i.e. it tesselated one eighth of the sphere.
! We then need additional symmetry operators to generate the right unique for
! symmetry portion of the surface, unless the group is D2h.
! That is what happens in the big if-else if statement below, where we construct
! the lsymop(0:7) array.
!
    use pedra_symmetry, only: get_pt

    integer(kind=regint_k) :: nv, nt, its, nvert, numts
    integer(kind=regint_k) :: jtr(numts, *)
    real(kind=dp) :: cv(nvert, *)

    integer(kind=regint_k) :: i, isymop, ii, jj, k, i_tmp
    logical :: lsymop(0:7)
    character(len=3) :: group_name

    group_name = group%group_name

! Symmetry operations in the Abelian groups
!      SYMOP(0) = ' E '
!      SYMOP(1) = 'Oyz'
!      SYMOP(2) = 'Oxz'
!      SYMOP(3) = 'C2z'
!      SYMOP(4) = 'Oxy'
!      SYMOP(5) = 'C2y'
!      SYMOP(6) = 'C2x'
!      SYMOP(7) = ' i '

    do i = 0, 7
        lsymop(i) = .false.
    end do
    ! Every group has the identity
    lsymop(0) = .true.
    if (group_name == 'C1 ') then
            ! C1 has 0 generators. Use all three planes of reflection
            lsymop(1) = .true.
            lsymop(2) = .true.
            lsymop(4) = .true.
    else if (group_name == 'Cs ') then
            ! Cs has 1 generator. Use Oxz plane and inversion
            lsymop(2) = .true.
            lsymop(7) = .true.
    else if (group_name == 'C2 ') then
            ! C2 has 1 generator. Use Oxz and inversion
            lsymop(2) = .true.
            lsymop(7) = .true.
    else if (group_name == 'Ci ') then
            ! Ci has 1 generator. Use Oxz and Oxy planes
            lsymop(2) = .true.
            lsymop(4) = .true.
    else if (group_name == 'C2h') then
            ! C2h has 2 generators. Use Oxz plane
            lsymop(2) = .true.
    else if (group_name == 'D2 ') then
            ! D2 has 2 generators. Use Oyz plane
            lsymop(1) = .true.
    else if (group_name == 'C2v') then
            ! C2v has 2 generators. Use inversion
            lsymop(7) = .true.
    else if (group_name /= 'D2h') then
            ! D2h has 3 generators. No need to select other operations for the replication.
            ! If we get here it  means something went awry before...
            write(lvpri, '(a)') "Check symmetry group."
            pedra_error_code = 11
            stop
    end if

    ! We DO NOT want to include the identity operator in this loop.
    ! That would cause a cavity twice as big as the correct one to be generated!
    do isymop = 1, 7
        if (lsymop(isymop)) then
            ! Replication of vertices
            do i = 1, nv
                ii = i + nv
                do k = 1, 3
                    i_tmp = 2_regint_k**(k - 1_regint_k)
                    cv(ii, k) = get_pt(iand(isymop, i_tmp)) * cv(i, k)
                end do
            end do
            ! Replication of topology
            do i = 1, nt
                ii = i + nt
                jj = nv
                do k = 1, 3
                    jtr(ii, k) = jtr(i, k) + jj
                end do
            end do
            ! Update indices
            nt  = nt  * 2
            nv  = nv  * 2
            its = its * 2
        end if
    end do

    end subroutine prerep

    subroutine pcmtns(vmat, geom, amass, katom)

    use pedra_utils, only: wlkdin

#include "pcm_mxcent.inc"
#include "pcm_pcmdef.inc"
#include "pcm_pcm.inc"

    integer(kind=regint_k), intent(in) :: katom
    real(kind=dp), intent(inout) :: geom(katom, 3)
    real(kind=dp),    intent(in) :: amass(katom)
    real(kind=dp), intent(inout) :: vmat(3, 3)

    real(kind=dp) :: eigval(3), eigvec(3, 3), tinert(3, 3)
    real(kind=dp) :: angmom(3), omegad(3), eiginv(3, 3), scal(3)
    integer(kind=regint_k) :: iax(6)
    logical :: planar, linear
    integer(kind=regint_k) :: i, j, k, jax, nmax
    integer(kind=regint_k) :: nopax, nshift
    real(kind=dp) :: com(3), total_mass
    real(kind=dp) :: dij

    iax = [1, 2, 3, 3, 2, 1]

    angmom = [1.0d0, 1.0d0, 1.0d0]

    ! Calculate center-of-mass (com) and translate
    total_mass = sum(amass)
    do i = 1, 3
      com(i) = 0.0_dp
      do j = 1, katom
        com(i) = com(i) + geom(j, i) * amass(j)
      end do
      com(i) = com(i) / total_mass
    end do
    ! Translation to COM
    do i = 1, katom
      geom(i, :) = geom(i, :) - com(:)
    end do

    call wlkdin(geom, amass, nesfp, angmom, tinert, omegad, eigval, eigvec, .true., planar, linear)

    do i = 1, 3
        do j = 1, 3
            eiginv(i,j) = eigvec(j,i)
        end do
    end do

    nopax = group%nr_rotations + group%nr_reflections
    if (nopax >= 3) then
        do i = 1, 3
            do j = 1, 3
                dij = 0.0d0
                if (i == j) then
                        dij = 1.0d0
                end if
                vmat(j, iax(group%jsop(i))) = dij
            end do
        end do
    elseif (nopax >= 1) then
        jax = iax(group%jsop(1))
        scal(1) = abs(eiginv(1,jax))
        scal(2) = abs(eiginv(2,jax))
        scal(3) = abs(eiginv(3,jax))
        nmax = 1
        do j = 2,3
            if (scal(j) > scal(nmax)) nmax = j
        end do
        nshift = mod(nmax-1_regint_k,3_regint_k)
        do i = 0,2
            k = mod(i + nshift, 3_regint_k) + 1
            do j =1,3
                vmat(i+1,j) = eiginv(k,j)
            end do
        end do
    elseif (nopax == 0) then
        do i = 1,3
            do j = 1,3
                vmat(i,j) = eiginv(i,j)
            end do
        end do
    else
        pedra_error_code = 9
        stop
    end if

    end subroutine pcmtns

    subroutine ordpcm(nts, xtscor, ytscor, ztscor, as)
!
! Performs some kind of check on tesserae
!
    integer(kind=regint_k), intent(in) :: nts
    real(kind=dp), intent(in) :: xtscor(nts), ytscor(nts), ztscor(nts), as(nts)

    logical :: lchk, lswtch
    real(kind=dp) :: xbak, ybak, zbak, abak
    integer(kind=regint_k) :: ibak, i, j, ii
    ! Some scratch space
    real(kind=dp) :: privec(4, nts) ! Tesserae centers and area
    integer(kind=regint_k) :: idxpri(nts)    ! Index of tessera i at i-th position

    lswtch = .false.
    privec = 0.0d0
    idxpri = 0

    do i = 1, nts
        ! Transfer coordinates, areas and index to scratch space.
        privec(1,i) = xtscor(i)
        privec(2,i) = ytscor(i)
        privec(3,i) = ztscor(i)
        privec(4,i) = as(i)
        idxpri(i) = i
    end do

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
        end do
    end if

    write(lvpri, '(a)') "Tess. #      x (AU)              y (AU)             "&
    //" z (AU)             a (AU^2)"
    write(lvpri, '(a)') "--------------------------------------------------"  &
    //"----------------------------------"
    do i = 1, nts
        ! This writes centers and areas of tesserae to lvpri (PEDRA.OUT)
        write(lvpri, '(i4, 4f20.14)') i, (privec(j,i), j=1, 4)
    end do

    end subroutine ordpcm

    logical function chktss(x1, y1, z1, a1, x2, y2, z2, a2)

    real(kind=dp) :: x1, y1, z1, a1, x2, y2, z2, a2

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
    end if

    end function chktss

    character(16) function typlab(i)

    integer(kind=regint_k), intent(in) :: i

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
