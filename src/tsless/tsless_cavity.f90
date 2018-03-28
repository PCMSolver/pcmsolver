!
! PCMSolver, an API for the Polarizable Continuum Model
! Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

!> Originally written by Christian Silvio Pomelli (ca. 2003-2004)
!> Introduced in PCMSolver by Christian Silvio Pomelli (ca. 2014)
!> Wrapped in a F90 module by Roberto Di Remigio (2015)
module tsless_cavity

use, intrinsic :: iso_c_binding
use tsless_precision
use tsless_symmetry
use tsless_weight_function

implicit none

public tsless_driver

private

contains

    !> \brief Interface function to TsLess cavity generation
    !> \author Roberto Di Remigio, Christian S. Pomelli
    !> \year 2014, 2015
    !> \param[in] maxts maximum number of tesserae allowed
    !> \param[in] maxsph maximum number of spheres allowed
    !> \param[in] maxvert maximum number of vertices allowed
    !> \param[out] nesfp number of spheres (original + added)
    !> \param[out] nts number of generated tesserae
    !> \param[out] ntsirr number of generated irreducible tesserae
    !> \param[out] addsph number of added spheres
    !> \param[out] xtscor x-coordinate of tesserae centers
    !> \param[out] ytscor y-coordinate of tesserae centers
    !> \param[out] ztscor z-coordinate of tesserae centers
    !> \param[out] ar area of the tessera
    !> \param[out] xsphcor x-coordinate of the sphere center the tessera belongs to
    !> \param[out] ysphcor y-coordinate of the sphere center the tessera belongs to
    !> \param[out] zsphcor z-coordinate of the sphere center the tessera belongs to
    !> \param[out] rsph radii of the sphere the tessera belongs to, i.e. its curvature
    !> \param[out] xe x-coordinate of the sphere center
    !> \param[out] ye y-coordinate of the sphere center
    !> \param[out] ze z-coordinate of the sphere center
    !> \param[out] rin radius of the spheres
    !> \param[in] masses atomic masses (for inertia tensor formation in TSLESS)
    !> \param[in] nr_gen number of symmetry generators
    !> \param[in] gen1 first generator
    !> \param[in] gen2 second generator
    !> \param[in] gen3 third generator
    !> \param[in] avgArea average tesserae area
    !> \param[in] dmin mininal distance between sampling points
    !> \param[in] nord maximum order of continuous derivative of weight function
    !> \param[in] ifun whether to use the normalized or unnormalized form of the weight function
    !> \param[in] rsolv solvent probe radius
    !> \param[in] work scratch space
    subroutine tsless_driver(maxts, maxsph, maxvert, nesfp, nts, ntsirr, addsph, &
                             xtscor, ytscor, ztscor, ar,         &
                             xsphcor, ysphcor, zsphcor, rsph,    &
                             xe, ye, ze, rin,                    &
                             masses,                             &
                             nr_gen, gen1, gen2, gen3,           &
                             avgArea, dmin, nord, ifun, rsolv, work) &
                             bind(C, name='tsless_driver')

    !> Passed variables
    integer(c_size_t), intent(in) :: maxts, maxsph, maxvert
    integer(c_int),    intent(in) :: nesfp
    real(c_double),   intent(out) :: xtscor(maxts), ytscor(maxts), ztscor(maxts), ar(maxts)
    real(c_double),   intent(out) :: xsphcor(maxts), ysphcor(maxts), zsphcor(maxts), rsph(maxts)
    real(c_double),    intent(in) :: xe(nesfp), ye(nesfp), ze(nesfp), rin(nesfp)
    real(c_double),    intent(in) :: masses(nesfp)
    real(c_double),    intent(in) :: work(maxts*maxsph)
    real(c_double),    intent(in) :: avgArea, rsolv, dmin
    integer(c_int),    intent(in) :: nord, ifun
    integer(c_int),   intent(out) :: nts, ntsirr, addsph
    integer(c_int),    intent(in) :: nr_gen, gen1, gen2, gen3
    !> Parameters
    real(kind=dp), parameter :: pi = acos(-1.0_dp)
    real(kind=dp), parameter :: tpi = 2.0_dp * pi
    real(kind=dp), parameter :: fpi = 4.0_dp * pi
    !> Local variables
    integer :: print_unit
    integer :: i, j, k, ntssp
    logical :: tsless_file_exists, tsless_file_opened
    integer(kind=regint_k) :: nesf_orig
    type(weight_function) :: weights
    real(kind=dp) :: x, y, z, r, dscale(3)
    real(kind=dp) :: dist, dist2, rho, t, w0
    type(point_group) :: pgroup

    print_unit = 121221_regint_k
    tsless_file_exists = .false.
    inquire(file = 'TSLESS.OUT', exist = tsless_file_exists, opened = tsless_file_opened)
    if (.not. tsless_file_opened) then
        open(print_unit,              &
            file = 'TSLESS.OUT',      &
            status = 'replace',       &
            form = 'formatted',       &
            access = 'sequential')
    end if

    !> Save original number of spheres
    nesf_orig = nesfp
    !> Create point group
    pgroup = build_point_group(nr_gen, gen1, gen2, gen3, print_unit)
    write(print_unit, *) 'Passed from host: nord, dmin, avgarea', nord, dmin, avgarea

    ! version 0.1 10/1/2014: plain implementation
    ! added sphere positions, radii and derivatives (TBI)
    ! generation and weighting of sampling points

    ! Generate weight function
    weights = create_weight_function(nderiv=nord, normalization=ifun, dmin=dmin, printer=print_unit)
    nts = 0_regint_k

    ! main loop on spheres (in future will implement linear scaling)
    LoopOnSpheres: do i = 1, nesfp
        ! point generation with weighting according to Leopardi scheme [ref]
        ntssp = int(fpi * rin(i)**2 / avgArea) + 1_regint_k
        !> Original weight factor
        w0 = fpi * rin(i)**2 / ntssp
        !> Interaction radius
        rho = sqrt(w0 * pi)
        write(print_unit, '(a, f12.5, a, f12.5)') 'Original weight factor ', w0, ' and interaction radius ', rho
        x = xe(i)
        y = ye(i)
        z = ze(i)
        r = rin(i)

        ! generate sphere sampling points
        write(print_unit, '(a, i4, a, i4, a)') 'Generating ', ntssp, ' sampling points on sphere', i, ' according to Leopardi'
        write(print_unit, '(a, f12.5, a, f12.5, a, f12.5, a, f12.5)') 'Center = (', x, ', ', y, ', ', z, ') Radius = ', r
        call leopardi_grid(ntssp, nts, x, y, z, r, xtscor, ytscor, ztscor, print_unit)
        write(print_unit, '(a)') 'Sampling points generated, now refining for intersections between spheres'

        ! loop on sphere sampling points
        LoopOnSamplingPoints: do j = nts + 1, nts + ntssp
            ar(j) = w0
            ! loop on intersecting spheres
            LoopOnIntersectingSpheres: do k = 1, nesfp
                if (k .eq. i) then
                    goto 10
                endif
                dist2 = (xe(k) - xtscor(j))**2 + (ye(k) - ytscor(j))**2 + (ze(k) - ztscor(j))**2
                dist = sqrt(dist2)
                t = (dist - rin(k)) / rho
                write(print_unit, *) 't = ', t
                if (t .le. dmin) then
                    ar(j) = -1.0_dp
                    write(print_unit, *) 'point ', j, ' was deleted'
                    goto 20
                else if(t .le. 1.0_dp) then
                    write(print_unit, *) 'point ', j, ' was scaled'
                    write(print_unit, *) 'Original value ', ar(j)
                    dscale = scaling_factor(weights, t, dmin)
                    ar(j) = ar(j) * dscale(1)
                    write(print_unit, *) 'Current value ', ar(j), ' Scaling factor ', dscale(1)
                endif
 10          continue
            end do LoopOnIntersectingSpheres
 20       continue

           ! save points with non zero weights
           SavePoints: if (ar(j) .gt. 0.0_dp) then
               nts = nts + 1_regint_k
               xtscor(nts) = xtscor(j)
               ytscor(nts) = ytscor(j)
               ztscor(nts) = ztscor(j)
               ar(nts) = ar(j)
               xsphcor(nts) = x
               ysphcor(nts) = y
               zsphcor(nts) = z
               rsph(nts) = r
           end if SavePoints
        end do LoopOnSamplingPoints
    end do LoopOnSpheres

    !> Number of irreducible sampling points, hardcoded to nts for the moment
    ntsirr = nts
    !> Number of added spheres
    addsph = nesfp - nesf_orig

    call destroy_weight_function(weights)

    write(print_unit, *) ' Surface = ', cavity_surface(nts, ar), ' volume = ', cavity_volume(nts, rsph, ar)
    close(print_unit)

    call write_visualization_file(nts, nesfp, xtscor, ytscor, ztscor, ar, xe, ye, ze, rin)

    end subroutine tsless_driver

    !> \brief Create and translate Leopardi points
    !> \author Christian S. Pomelli
    !> \param[in] ntssp number of points on the sphere
    !> \param[in] nts   a counter for the total number of points generated
    !> \param[in] x x-coordinate of sphere center
    !> \param[in] y y-coordinate of sphere center
    !> \param[in] z z-coordinate of sphere center
    !> \param[in] r radius of sphere
    !> \param[in] xtscor x-coordinates of points on the sphere
    !> \param[in] ytscor y-coordinates of points on the sphere
    !> \param[in] ztscor z-coordinates of points on the sphere
    !> \param[in] print_unit
    subroutine leopardi_grid(ntssp, nts, x, y, z, r, xtscor, ytscor, ztscor, print_unit)

    use, intrinsic :: iso_fortran_env, only: output_unit

    !> Passed variables
    integer(kind=regint_k), intent(in) :: ntssp, nts
    integer, optional,      intent(in) :: print_unit
    real(kind=dp),          intent(in) :: x, y, z, r
    real(kind=dp),       intent(inout) :: xtscor(*), ytscor(*), ztscor(*)
    !> Parameters
    real(kind=dp), parameter :: pi = acos(-1.0_dp)
    real(kind=dp), parameter :: tpi = 2.0_dp * pi
    real(kind=dp), parameter :: fpi = 4.0_dp * pi
    !> Local variables
    integer :: i, j, k, n, print_out
    real(kind=dp) :: thetaf(ntssp), ypsi(ntssp), theta(ntssp), apsi(ntssp), m(ntssp)
    real(kind=dp) :: vr, thetac, delta1, deltaf
    real(kind=dp) :: thetam, costm, sintm, omegaj

    if (present(print_unit)) then
       print_out = print_unit
    else
       print_out = output_unit
    end if

    thetaf = 0.0_dp
    ypsi = 0.0_dp
    theta = 0.0_dp
    apsi = 0.0_dp
    m = 0.0_dp
    ! Notation consistent with [ref]
    vr = fpi / ntssp
    thetac = 2.0_dp * (asin(sqrt(vr / fpi)))
    delta1 = sqrt(vr)

    ! Polar points
    ! North pole
    xtscor(nts+1) = x
    ytscor(nts+1) = y
    ztscor(nts+1) = z + r
    ! South pole
    xtscor(nts+ntssp) = x
    ytscor(nts+ntssp) = y
    ztscor(nts+ntssp) = z - r

    ! two points case (rare)
    if (ntssp == 2) return

    ! Number of collars
    n = int((pi - 2.0_dp * thetac) / delta1 + 0.5_dp)
    if (n == 0) n = 1
    deltaf = (pi - 2.0_dp * thetac) / (1.0_dp * n)

    ! Colatitudes of collars
    CollarColatitudes: do i = 1, n + 1
        thetaF(i) = thetac + (i - 1) * deltaf
    end do CollarColatitudes

    ! Collars building
    ypsi(1) = fpi*((sin(thetaf(2) / 2.0_dp))**2 - &
                   (sin(Thetaf(1) / 2.0_dp))**2) / vr
    m(1) = int(ypsi(1) + 0.5_dp)
    apsi(1) = ypsi(1) - m(1)
    Collars: do i = 2, n
        ypsi(i) = fpi * ((sin(thetaf(i+1) / 2.0_dp))**2 - &
            (sin(thetaf(i  ) / 2.0_dp))**2) / vr
        m(i) = int(ypsi(i) + apsi(i-1) + 0.5_dp)
        apsi(i) = ypsi(i) - m(i) + apsi(i-1)
    end do Collars

    j = 1
    do i = 1, n + 1
        theta(i) = 2.0_dp * asin(sqrt((vr * j) / fpi))
        j = j + int(m(i) + 0.5_dp)
    end do

    k = 2
    do i = 1, n
        thetam = (theta(i) + theta(i+1)) / 2.0_dp
        costm = cos(thetam)
        sintm = sin(thetam)
        TranslatePoints: do j = 1, int(m(i) + 0.5_dp)
            omegaj = tpi * j / (1.0_dp * m(i))
            xtscor(nts+k) = x + r * sintm * cos(omegaj)
            ytscor(nts+k) = y + r * sintm * sin(omegaj)
            ztscor(nts+k) = z + r * costm
            k=k+1
        end do TranslatePoints
    end do

    end subroutine leopardi_grid

    subroutine write_visualization_file(nts, nesfp, xtscor, ytscor, ztscor, ar, xe, ye, ze, rin)

    !> writes a script for VMD

    !> Passed variables
    integer(kind=regint_k), intent(in) :: nts, nesfp
    real(kind=dp), intent(in) :: xtscor(nts), ytscor(nts), ztscor(nts), ar(nts), xe(nesfp), ye(nesfp), ze(nesfp), rin(nesfp)
    !> Local variables
    integer(kind=regint_k) :: file_unit
    logical :: file_exists, file_opened
    integer :: k

    file_unit = 221221_regint_k
    file_exists = .false.
    inquire(file = 'tsless_visual.vmd', exist = file_exists, opened = file_opened)
    if (.not. file_opened) then
        open(file_unit,                 &
            file = 'tsless_visual.vmd', &
            status = 'replace',         &
            form = 'formatted',         &
            access = 'sequential')
    end if

    write(file_unit, *) 'mol load graphics name'

    do k = 1, nesfp
       write(file_unit,1000)  xe(k), ye(k), ze(k), rin(k)
    enddo

    write(file_unit, *) 'graphics top color 3'

    do k = 1, nts
        write(file_unit,1000)  xtscor(k), ytscor(k), ztscor(k), 0.1_dp
    end do
    write(file_unit, *)

    close(file_unit)

1000 format ('graphics top sphere {', 3F10.4, '} radius', F10.4, ' resolution 120')

    end subroutine write_visualization_file

    pure function cavity_surface(nts, ar) result(surface)

    !> Passed variables
    integer(kind=regint_k), intent(in) :: nts
    real(kind=dp), intent(in) :: ar(nts)
    real(kind=dp) :: surface
    !> Local variables
    integer :: i

    surface = 0.0_dp
    do i = 1, nts
        surface = surface + ar(i)
    end do

    end function cavity_surface

    pure function cavity_volume(nts, rsph, ar) result(volume)

    !> Passed variables
    integer(kind=regint_k), intent(in) :: nts
    real(kind=dp), intent(in) :: rsph(nts), ar(nts)
    real(kind=dp) :: volume
    !> Local variables
    integer :: i

    volume = 0.0_dp
    do i = 1, nts
        volume = volume + ar(i) * rsph(i)
    end do
    volume = volume / 3.0_dp

    end function cavity_volume

end module tsless_cavity
