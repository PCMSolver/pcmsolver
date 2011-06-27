      function gsfera
     $     (eps,eps2,xs,ys,zs,rs,xi,yi,zi,xj,yj,zj)
      implicit real*8 (a-h,o-z)
      complex*16 eps2, gsfera, coefl
Clf      print *,eps,eps2
Clf      print *, 'sfera', xs,ys,zs,rs
Clf      print *, 'primo punto', xi,yi,zi
Clf      print *, 'secondo punto', xj,yj,zj
      di=sqrt((xi-xs)**2+(yi-ys)**2+(zi-zs)**2)
      dj=sqrt((xj-xs)**2+(yj-ys)**2+(zj-zs)**2)
      dim=rs**2/di
      xim=dim*(xi-xs)/di+xs
      yim=dim*(yi-ys)/di+ys
      zim=dim*(zi-zs)/di+zs
      qim=rs/di
      gc=qim/sqrt((xj-xim)**2+(yj-yim)**2+(zj-zim)**2)
      gc=gc-qim/sqrt((xj-xs)**2+(yj-ys)**2+(zj-zs)**2)
      cost=(xj-xs)*(xi-xs)+(yj-ys)*(yi-ys)+(zj-zs)*(zi-zs)
      cost=cost/(di*dj)
      aa=(rs*rs/(di*dj))**4.
      maxl=int(800./(abs(eps2+eps))**0.4*aa)
      maxl=max(maxl,10)
Clf      write (6,*) "maxl",maxl
      arg=rs*rs/(di*dj)
      argl=rs/(di*dj)
      gsfera2=0.d0
      m=0
      do l=1,maxl
        argl=argl*arg
        coefl=(eps2-eps)*float(l)/((eps2+eps)*float(l)+1.)
        coefl=coefl-(eps2-eps)/(eps2+eps)
        gsfera=gsfera+coefl*argl*plgndr(l,m,cost)
      enddo
      gsfera=gsfera+gc*(eps2-eps)/(eps2+eps)
Clf: WARNING: sign change and divided by epsilon of solvent!!!!!!
      gsfera = -gsfera/eps
C Solo per debugging
C      arg=rs*rs/(di*dj)
C      argl=rs/(di*dj)
C      gsfera=0.d0
C      m=0
C      do l=1,5000
C        argl=argl*arg
C        coefl=(eps2-eps)*float(l)/((eps2+eps)*float(l)+1.)
C        gsfera=gsfera+coefl*argl*plgndr(l,m,cost)
C      enddo 
Clf      write (6,*) "gsfera2",gsfera
      return
      end

      subroutine gsfera_cpp(epssol, epsre, epsim, radius, ps, p1, p2, 
     $     greenre, greenim)
      implicit real*8 (a-h, o-z)
      complex*16 gsfera, green, eps2, ui
      double precision p1(*), p2(*), ps(*)
      ui = (0.0, 1.0)
      eps2 = epsre + epsim * ui
      green = gsfera(epssol, eps2, ps(1), ps(2), ps(3), radius,
     &               p1(1), p1(2), p1(3),
     &               p2(1), p2(2), p2(3))
      greenre = dble(green)
      greenim = dimag(green)
      return
      end
