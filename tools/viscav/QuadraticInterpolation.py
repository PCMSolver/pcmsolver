# Quadratic interpolation in two dimensions through nine 3D points
import math

def interpolate1D(s, p0, p1, p2):
    res = [0.0, 0.0, 0.0]
    for i in xrange(3):
        c = p0[i]
        b = -3.0*p0[i] + 4.0*p1[i] - p2[i]
        a = 2.0*p0[i] - 4.0*p1[i] + 2.0*p2[i]
        res[i] = a*s*s + b*s + c
    return res

def gen2DInterpolationCoefficients(mesh):

    # Assume mesh points are located at u=[0.0, 0.5, 1.0] and v=[0.0, 0.5, 1.0]

    a = [[[0.0, 0.0, 0.0] for i in xrange(3)] for j in xrange(3)]

    for i in xrange(3):
        pi = lambda j: mesh[j/3][j%3][i]
        a[0][0][i] = 4.0*pi(0) - 8.0*pi(1) + 4.0*pi(2) - 8.0*pi(3) +16.0*pi(4) - 8.0*pi(5) + 4.0*pi(6) - 8.0*pi(7) + 4.0*pi(8)
        a[0][1][i] =-6.0*pi(0) + 8.0*pi(1) - 2.0*pi(2) +12.0*pi(3) -16.0*pi(4) + 4.0*pi(5) - 6.0*pi(6) + 8.0*pi(7) - 2.0*pi(8)
        a[0][2][i] = 2.0*pi(0) - 4.0*pi(3) + 2.0*pi(6)

        a[1][0][i] =-6.0*pi(0) +12.0*pi(1) - 6.0*pi(2) + 8.0*pi(3) -16.0*pi(4) + 8.0*pi(5) - 2.0*pi(6) + 4.0*pi(7) - 2.0*pi(8)
        a[1][1][i] = 9.0*pi(0) -12.0*pi(1) + 3.0*pi(2) -12.0*pi(3) +16.0*pi(4) - 4.0*pi(5) + 3.0*pi(6) - 4.0*pi(7) + 1.0*pi(8)
        a[1][2][i] =-3.0*pi(0) + 4.0*pi(3) - 1.0*pi(6)

        a[2][0][i] = 2.0*pi(0) - 4.0*pi(1) + 2.0*pi(2)
        a[2][1][i] =-3.0*pi(0) + 4.0*pi(1) - 1.0*pi(2)
        a[2][2][i] = pi(0)

    return a

def interpolate2D(u, v, coeffs):
    p = [0.0, 0.0, 0.0]
    tmp = [0.0, 0.0, 0.0]
    for i in xrange(3):
        ai = lambda x,y: coeffs[y][x][i]
        tmp[0] = ai(0,0)*v*v + ai(0,1)*v + ai(0,2)
        tmp[1] = ai(1,0)*v*v + ai(1,1)*v + ai(1,2)
        tmp[2] = ai(2,0)*v*v + ai(2,1)*v + ai(2,2)
        p[i] = tmp[0]*u*u + tmp[1]*u + tmp[2]

    return p

def interpolate2DNormal(u, v, coeffs):
    du = [0.0, 0.0, 0.0]
    dv = [0.0, 0.0, 0.0]

    tmpu = [0.0, 0.0, 0.0]
    tmpv = [0.0, 0.0, 0.0]

    for i in xrange(3):
        ai = lambda x,y: coeffs[y][x][i]
        
        tmpu[0] = ai(0,0)*v*v + ai(0,1)*v + ai(0,2)
        tmpu[1] = ai(1,0)*v*v + ai(1,1)*v + ai(1,2)
        #tmpu[2] = ai(2,0)*v*v + ai(2,1)*v + ai(2,2)
        du[i] = tmpu[0]*2.0*u + tmpu[1]

        tmpv[0] = ai(0,0)*2.0*v + ai(0,1)
        tmpv[1] = ai(1,0)*2.0*v + ai(1,1)
        tmpv[2] = ai(2,0)*2.0*v + ai(2,1)
        dv[i] = tmpv[0]*u*u + tmpv[1]*u + tmpv[2]

    normal = [0.0, 0.0, 0.0]
    normal[0] = du[1]*dv[2] - du[2]*dv[1]
    normal[1] = du[2]*dv[0] - du[0]*dv[2]
    normal[2] = du[0]*dv[1] - du[1]*dv[0]

    length = math.sqrt(sum([n*n for n in normal]))
    res = [-n/length for n in normal]

    return res

