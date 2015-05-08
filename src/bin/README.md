In all cases the dielectric profile is parametrized by
the profile function:

     epsilon(r) = eR+eL/2 + eR-eL/2 * tanh((r-c)/w)

with eR = 2.0  (outside)
     eL = 80.0 (inside)
     c  = 100.0 (radius of the dielectric sphere)
     w  =  10.0 (width of the interface layer)

Naming:
    gf_spherical_CASEn.log
        - grid points
        - Green's function
        - singular part of the Green's function
        - image part of the Green's function
    gf_uniform_inside_CASEn.log
        - grid points
        - Green's function
    gf_uniform_outside_CASEn.log
        - grid points
        - Green's function

CASE1:
    source point at (4.0, 0.0, 85.0)
    probe point at  (0.0, 0.0, z) with z in 10.0, 300.0

CASE2:
    source point at (4.0, 0.0, 150.0)
    probe point at  (0.0, 0.0, z) with z in 10.0, 300.0

CASE3:
    source point at (4.0, 0.0, 100.0)
    probe point at  (0.0, 0.0, z) with z in 10.0, 300.0

CASE4:
    source point at (0.0, 0.0, 1.0)
    probe point at  (0.0, 0.0, z) with z in 10.0, 300.0
