C
C     simple input reader for cavity generator
C     written by Krzysztof Mozgawa, 2010
C
      subroutine generatecavity_cpp(xtscor,ytscor,ztscor,ar,xsphcor,
     &     ysphcor, zsphcor, rsph, nts, nesfp, xe, ye, ze, rin, avgArea,
     $     work, lwork)

      implicit real*8 (a-h, o-z)
      double precision xtscor(*), ytscor(*), ztscor(*)
      double precision xsphcor(*), ysphcor(*), zsphcor(*), rsph(*)
      double precision ar(*), xe(*), ye(*), ze(*), rin(*)
      integer nts, nesfp
      call makecav_(work,lwork,xtscor,ytscor,ztscor,ar,xsphcor,
     & ysphcor, zsphcor, rsph, nts, nesfp, xe, ye, ze, rin, avgArea)
C      print *, 'cavity created!'
      end

