    module diagonal
    
    implicit none
 
    public diagonal_anisotropic
    public diagonal_ionic

    contains

    subroutine diagonal_anisotropic(vert, centr, euphi, euthe, eupsi, eps1, eps2, eps3) 
    ! This is (almost) the original subroutine for the calculation of
    ! the diagonal elements of S and D for an anisotropic dielectric.
    ! Taken from GAMESS, possibly written by Eric Cancès in 1997-1998
    implicit none

#include "pcm_mxcent.h"
#include "pcm_pcmdef.h"
#include "pcm_pcm.h"
      
    real(8) :: vert(mxts, 10, 3), centr(mxts, 10, 3)
    real(8) :: euphi, euthe, eupsi
    real(8) :: eps1, eps2, eps3
    ! Parameters
    real(8), parameter :: zero = 0.0d0
    real(8), parameter :: two  = 2.0d0
    real(8), parameter :: pi   = acos(-1.0d0)
    real(8), parameter :: tpi  = 2.0d0 * pi
    real(8), parameter :: fpi  = 4.0d0 * pi
    real(8), parameter :: to_radians = acos(-1.0d0) / 180.0d0
    integer, parameter :: ng16pts = 8, ng64pts = 32
    ! Local variables
    real(8) :: theta(10), phi(10), phinumb(10)
    integer :: numb(10) 
    real(8) :: xgp16pts(8), wgp16pts(8), xgp64pts(32), wgp64pts(32)
    real(8) :: rot(3, 3)
    real(8) :: epsxx, epsxy, epsxz, epsyy, epsyz, epszz, epsm
    real(8) :: epsm1xx, epsm1xy, epsm1xz, epsm1yy, epsm1yz, epsm1zz
    real(8) :: csphi, snphi, csthe, snthe, cspsi, snpsi
    real(8) :: aph, bph, ath, pha, phb
    real(8) :: bb, cc, cosph, costh, cs, cotgthmax, dc
    real(8) :: ds, dcs, ocx, ocy, ocz, ph, rcc, rmin
    real(8) :: aa, rth, rtheps, sinph, sinth, th, tha, thb
    real(8) :: thmax, vecx, vecy, vecz, vx, vy, vz
    real(8) :: xcc, xx, xy, xz, ycc, yx, yy, yz, zcc, zx, zy, zz
    real(8) :: sse, dde, ssep, ddep, ssepp, ddepp, deteps
    real(8) :: rth_test, test_S, test_Sp, test_Spp ! Same nomenclature
    real(8) :: accum
    integer :: i, j, k, l, m, li, nedge, nph, nphsgn, nth, nthsgn, nv
    
    ! 16-points Gaussian rule: abscissas 
    xgp16pts = [9.894009350d-01, 9.445750231d-01, &
                8.656312024d-01, 7.554044084d-01, &
                6.178762444d-01, 4.580167777d-01, &
                2.816035508d-01, 9.501250984d-02]
    ! 16-points Gaussian rule: weights
    wgp16pts = [2.715245941d-02, 6.225352394d-02, &
                9.515851168d-02, 1.246289713d-01, &
                1.495959888d-01, 1.691565194d-01, &
                1.826034150d-01, 1.894506105d-01]
    ! 64-points Gaussian rule: abscissas
    xgp64pts = [9.993050417d-01, 9.963401168d-01, &  
                9.910133715d-01, 9.833362539d-01, &
                9.733268278d-01, 9.610087996d-01, &
                9.464113748d-01, 9.295691721d-01, &
                9.105221371d-01, 8.893154460d-01, &
                8.659993982d-01, 8.406292963d-01, & 
                8.132653151d-01, 7.839723589d-01, & 
                7.528199073d-01, 7.198818502d-01, & 
                6.852363131d-01, 6.489654713d-01, &
                6.111553552d-01, 5.718956462d-01, &
                5.312794640d-01, 4.894031457d-01, &
                4.463660173d-01, 4.022701580d-01, &
                3.572201583d-01, 3.113228720d-01, &
                2.646871622d-01, 2.174236437d-01, &
                1.696444204d-01, 1.214628193d-01, &
                7.299312179d-02, 2.435029266d-02]
    ! 64-points Gaussian rule: weights
    wgp64pts = [1.783280722d-03, 4.147033261d-03, &
                6.504457969d-03, 8.846759826d-03, &
                1.116813946d-02, 1.346304790d-02, &
                1.572603048d-02, 1.795171578d-02, &
                2.013482315d-02, 2.227017381d-02, &
                2.435270257d-02, 2.637746972d-02, &
                2.833967261d-02, 3.023465707d-02, &
                3.205792835d-02, 3.380516184d-02, &
                3.547221326d-02, 3.705512854d-02, &
                3.855015318d-02, 3.995374113d-02, &
                4.126256324d-02, 4.247351512d-02, &
                4.358372453d-02, 4.459055816d-02, &
                4.549162793d-02, 4.628479658d-02, &
                4.696818282d-02, 4.754016571d-02, &
                4.799938860d-02, 4.834476223d-02, &
                4.857546744d-02, 4.869095701d-02]
    
    accum = 0.0d0
    euphi = euphi * to_radians
    euthe = euthe * to_radians
    eupsi = eupsi * to_radians
    csphi=cos(euphi)
    snphi=sin(euphi)
    csthe=cos(euthe)
    snthe=sin(euthe)
    cspsi=cos(eupsi)
    snpsi=sin(eupsi)

    rot(1,1)=csphi*cspsi-snphi*csthe*snpsi
    rot(1,2)=-snphi*cspsi-csphi*csthe*snpsi
    rot(1,3)=snthe*snpsi
    rot(2,1)=csphi*snpsi+snphi*csthe*cspsi
    rot(2,2)=-snphi*snpsi+csphi*csthe*cspsi
    rot(2,3)=-snthe*cspsi
    rot(3,1)=snphi*snthe
    rot(3,2)=csphi*snthe
    rot(3,3)=csthe

    ! Cube root of the determinant
    deteps = eps1 * eps2 * eps3
    epsm=(eps1*eps2*eps3)**(1.0d+00/3.0d+00)
    ! Build permittivity tensor in molecule-fixed frame
    epsxx=(eps1*rot(1,1)*rot(1,1)+eps2*rot(2,1)*rot(2,1)+eps3*rot(3,1)*rot(3,1))
    epsxy=(eps1*rot(1,1)*rot(1,2)+eps2*rot(2,1)*rot(2,2)+eps3*rot(3,1)*rot(3,2))
    epsxz=(eps1*rot(1,1)*rot(1,3)+eps2*rot(2,1)*rot(2,3)+eps3*rot(3,1)*rot(3,3))
    epsyy=(eps1*rot(1,2)*rot(1,2)+eps2*rot(2,2)*rot(2,2)+eps3*rot(3,2)*rot(3,2))
    epsyz=(eps1*rot(1,2)*rot(1,3)+eps2*rot(2,2)*rot(2,3)+eps3*rot(3,2)*rot(3,3))
    epszz=(eps1*rot(1,3)*rot(1,3)+eps2*rot(2,3)*rot(2,3)+eps3*rot(3,3)*rot(3,3))
    ! Build inverse of permittivity tensor
    epsm1xx=(eps1**(-1)*rot(1,1)*rot(1,1)+eps2**(-1)*rot(2,1)*rot(2,1)+eps3**(-1)*rot(3,1)*rot(3,1))
    epsm1xy=(eps1**(-1)*rot(1,1)*rot(1,2)+eps2**(-1)*rot(2,1)*rot(2,2)+eps3**(-1)*rot(3,1)*rot(3,2))
    epsm1xz=(eps1**(-1)*rot(1,1)*rot(1,3)+eps2**(-1)*rot(2,1)*rot(2,3)+eps3**(-1)*rot(3,1)*rot(3,3))
    epsm1yy=(eps1**(-1)*rot(1,2)*rot(1,2)+eps2**(-1)*rot(2,2)*rot(2,2)+eps3**(-1)*rot(3,2)*rot(3,2))
    epsm1yz=(eps1**(-1)*rot(1,2)*rot(1,3)+eps2**(-1)*rot(2,2)*rot(2,3)+eps3**(-1)*rot(3,2)*rot(3,3))
    epsm1zz=(eps1**(-1)*rot(1,3)*rot(1,3)+eps2**(-1)*rot(2,3)*rot(2,3)+eps3**(-1)*rot(3,3)*rot(3,3))
    eps=epsm

    sse=zero
    dde=zero
    test_S = 0.0d0

    do i = 1, nts
    li=isphe(i)
    nv=30*(i-1)
!
! -- LOOP ON INTEGRATION POINTS
!
        ssep=0.0d+00
        ddep=0.0d+00
        test_Sp = 0.0d0

        xz=(xtscor(i)-xe(li))/re(li)
        yz=(ytscor(i)-ye(li))/re(li)
        zz=(ztscor(i)-ze(li))/re(li)
        azimuth = atan((ytscor(i) - ye(li))/(xtscor(i)-xe(li)))
        rmin=0.99d+00
        if (abs(xz).le.rmin) then
          rmin=abs(xz)
          xx=0.0d+00
          yx=-zz/sqrt(1.0d+00-xz*xz)
          zx=yz/sqrt(1.0d+00-xz*xz)
        end if
        if (abs(yz).le.rmin) then
          rmin=abs(yz)
          xx=zz/sqrt(1.0d+00-yz*yz)
          yx=0.0d+00
          zx=-xz/sqrt(1.0d+00-yz*yz)
        end if
        if (abs(zz).le.rmin) then
          xx=yz/sqrt(1.0d+00-zz*zz)
          yx=-xz/sqrt(1.0d+00-zz*zz)
          zx=0.0d+00
        end if
        xy=yz*zx-yx*zz
        yy=zz*xx-zx*xz
        zy=xz*yx-xx*yz


        ! Clean-up heap-crap        
        theta = 0.0d0
        phi   = 0.0d0
        numb  = 0
        phinumb = 0.0d0

        do k=1,nvert(i)
          vecx=vert(i,k,1)-xe(li)
          vecy=vert(i,k,2)-ye(li)
          vecz=vert(i,k,3)-ze(li)
          dc=(vecx*xz+vecy*yz+vecz*zz)/re(li)
          if(dc.ge.1.0d+00)dc=1.0d+00
          if(dc.le.-1.0d+00)dc=-1.0d+00
          theta(k)=acos(dc)
          dcs=(vecx*xx+vecy*yx+vecz*zx)/(re(li)*sin(theta(k)))
          if(dcs.ge.1.0d+00)dcs=1.0d+00
          if(dcs.le.-1.0d+00)dcs=-1.0d+00
          phi(k)=acos(dcs)
          snphi=(vecx*xy+vecy*yy+vecz*zy)/(re(li)*sin(theta(k)))
          if (snphi.le.0.0d+00) phi(k)=tpi-phi(k)
        enddo
       
        do k=2,nvert(i)
          phi(k)=phi(k)-phi(1)
          if (phi(k).lt.0.0d+00) phi(k)=tpi+phi(k)
        enddo

        xx=xx*cos(phi(1))+xy*sin(phi(1))
        yx=yx*cos(phi(1))+yy*sin(phi(1))
        zx=zx*cos(phi(1))+zy*sin(phi(1))
        xy=yz*zx-yx*zz
        yy=zz*xx-zx*xz
        zy=xz*yx-xx*yz
        
        phi(1)=0.0d+00
        numb(1)=1
        numb(2)=2
        phinumb(1)=phi(1)
        phinumb(2)=phi(2)
        do 210 k=3,nvert(i)
          do l=2,k-1
            if (phi(k).lt.phinumb(l)) then
              do m=1,k-l
                numb(k-m+1)=numb(k-m)
                phinumb(k-m+1)=phinumb(k-m)
              enddo
              numb(l)=k
              phinumb(l)=phi(k)
              go to 210
            end if
          enddo
          numb(k)=k
          phinumb(k)=phi(k)
  210   continue
        numb(nvert(i)+1)=numb(1)
        phinumb(nvert(i)+1)=tpi
!
! -- LOOP ON GAUSS POINTS
!
        do nedge=1,nvert(i)
          aph=(phinumb(nedge+1)-phinumb(nedge))/2.0d+00
          bph=(phinumb(nedge+1)+phinumb(nedge))/2.0d+00
          tha=theta(numb(nedge))
          thb=theta(numb(nedge+1))
          pha=phinumb(nedge)
          phb=phinumb(nedge+1)
          ocx=(centr(i,nedge,1)-xe(li))/re(li)
          ocy=(centr(i,nedge,2)-ye(li))/re(li)
          ocz=(centr(i,nedge,3)-ze(li))/re(li)
          rcc=ocx**2+ocy**2+ocz**2
          
          do nph=1,ng64pts
          do nphsgn=0,1
           ph=(2*nphsgn-1)*aph*xgp64pts(nph)+bph
           cosph=cos(ph)
           sinph=sin(ph)
          if (rcc.lt.1.0d-07) then
           cotgthmax=(sin(ph-pha)/tan(thb)+sin(phb-ph)/tan(tha))/sin(phb-pha)
           thmax=atan(1.0d+00/cotgthmax)
          else
           xcc=xx*ocx+yx*ocy+zx*ocz
           ycc=xy*ocx+yy*ocy+zy*ocz
           zcc=xz*ocx+yz*ocy+zz*ocz
           aa=(xcc*cosph+ycc*sinph)**2+zcc**2
           bb=-zcc*rcc
           cc=rcc**2-(xcc*cosph+ycc*sinph)**2
           ds=bb**2-aa*cc
           if(ds.lt.zero)ds=zero
           cs=(-bb+sqrt(ds))/aa
           if(cs.gt.1.0d+00)cs=1.0d+00
           if(cs.lt.-1.0d+00)cs=1.0d+00
           thmax=acos(cs)
          end if
           if(thmax.lt.1.0d-08) go to 33
           ath=thmax/two
           ssepp=0.0d+00
           ddepp=0.0d+00
           test_Spp = 0.0d0
           
           do nth=1,ng16pts ! Inner 16-points rule
           do nthsgn=0,1
            th=(2*nthsgn-1)*ath*xgp16pts(nth)+ath
            costh=cos(th)
            sinth=sin(th)
            vx=xx*sinth*cosph+xy*sinth*sinph+xz*(costh-1.0d+00)
            vy=yx*sinth*cosph+yy*sinth*sinph+yz*(costh-1.0d+00)
            vz=zx*sinth*cosph+zy*sinth*sinph+zz*(costh-1.0d+00)
            rth=sqrt(2*(1-costh))
            rtheps=sqrt(epsm1xx*vx*vx+2*epsm1xy*vx*vy+2*epsm1xz*vx*vz &
                                     +  epsm1yy*vy*vy+2*epsm1yz*vy*vz &
                                                     +  epsm1zz*vz*vz)
            ssepp=ssepp+(re(li)/(rtheps))*sinth*ath*wgp16pts(nth)
            ddepp=ddepp-(rth**2/(2*rtheps**3))*sinth*ath*wgp16pts(nth)
           enddo ! Close loop on nthsgn (n-theta-sign)
           enddo ! Close loop on nth (n-theta)
           ssep=ssep+ssepp*aph*wgp64pts(nph)
           ddep=ddep+ddepp*aph*wgp64pts(nph)
           
  33       continue
          enddo ! Close loop on nphsgn (n-phi-sign)
          enddo ! Close loop on nph (n-phi)
        enddo ! Close loop on nedge
        sse=ssep
        dde=ddep
        
        ! Divide by determinant
        sse = sse / deteps
        dde = dde / deteps
      
      end do

      end subroutine diagonal_anisotropic

      subroutine diagonal_ionic(vert, centr, eps, kappa) 
      ! This is (almost) the original subroutine for the calculation of 
      ! the diagonal elements of S and D for an anisotropic dielectric.
      ! Taken from GAMESS, possibly written by Eric Cancès in 1997-1998
      implicit none

#include "pcm_mxcent.h"
#include "pcm_pcmdef.h"
#include "pcm_pcm.h"
      
      real(8) :: vert(mxts, 10, 3), centr(mxts, 10, 3)                   
      real(8) :: eps, kappa
      ! Parameters
      real(8), parameter :: zero = 0.0d0
      real(8), parameter :: two  = 2.0d0
      real(8), parameter :: pi   = acos(-1.0d0)
      real(8), parameter :: tpi  = 2.0d0 * pi
      real(8), parameter :: fpi  = 4.0d0 * pi
      real(8), parameter :: to_radians = acos(-1.0d0) / 180.0d0
      integer, parameter :: ng16pts = 8, ng64pts = 32
      ! Local variables
      real(8) :: theta(10), phi(10), phinumb(10)
      integer :: numb(10) 
      real(8) :: xgp16pts(8), wgp16pts(8), xgp64pts(32), wgp64pts(32)
      real(8) :: aph, bph, ath, pha, phb
      real(8) :: bb, cc, cosph, costh, cs, cotgthmax, dc
      real(8) :: ds, dcs, ocx, ocy, ocz, ph, rcc, rmin
      real(8) :: aa, rth, rtheps, sinph, sinth, th, tha, thb
      real(8) :: thmax, vecx, vecy, vecz, vx, vy, vz
      real(8) :: xcc, xx, xy, xz, ycc, yx, yy, yz, zcc, zx, zy, zz
      real(8) :: sse, dde, ssep, ddep, ssepp, ddepp, deteps
      real(8) :: rth_test, test_S, test_Sp, test_Spp ! Same nomenclature
      real(8) :: accum
      integer :: i, j, k, l, m, li, nedge, nph, nphsgn, nth, nthsgn, nv
      
      ! 16-points Gaussian rule: abscissas 
      xgp16pts = [9.894009350d-01, 9.445750231d-01, &
                  8.656312024d-01, 7.554044084d-01, &
                  6.178762444d-01, 4.580167777d-01, &
                  2.816035508d-01, 9.501250984d-02]
      ! 16-points Gaussian rule: weights
      wgp16pts = [2.715245941d-02, 6.225352394d-02, &
                  9.515851168d-02, 1.246289713d-01, &
                  1.495959888d-01, 1.691565194d-01, &
                  1.826034150d-01, 1.894506105d-01]
      ! 64-points Gaussian rule: abscissas
      xgp64pts = [9.993050417d-01, 9.963401168d-01, &  
                  9.910133715d-01, 9.833362539d-01, &
                  9.733268278d-01, 9.610087996d-01, &
                  9.464113748d-01, 9.295691721d-01, &
                  9.105221371d-01, 8.893154460d-01, &
                  8.659993982d-01, 8.406292963d-01, & 
                  8.132653151d-01, 7.839723589d-01, & 
                  7.528199073d-01, 7.198818502d-01, & 
                  6.852363131d-01, 6.489654713d-01, &
                  6.111553552d-01, 5.718956462d-01, &
                  5.312794640d-01, 4.894031457d-01, &
                  4.463660173d-01, 4.022701580d-01, &
                  3.572201583d-01, 3.113228720d-01, &
                  2.646871622d-01, 2.174236437d-01, &
                  1.696444204d-01, 1.214628193d-01, &
                  7.299312179d-02, 2.435029266d-02]
      ! 64-points Gaussian rule: weights                                                               
      wgp64pts = [1.783280722d-03, 4.147033261d-03, &
                  6.504457969d-03, 8.846759826d-03, &
                  1.116813946d-02, 1.346304790d-02, &
                  1.572603048d-02, 1.795171578d-02, &
                  2.013482315d-02, 2.227017381d-02, &
                  2.435270257d-02, 2.637746972d-02, &
                  2.833967261d-02, 3.023465707d-02, &
                  3.205792835d-02, 3.380516184d-02, &
                  3.547221326d-02, 3.705512854d-02, &
                  3.855015318d-02, 3.995374113d-02, &
                  4.126256324d-02, 4.247351512d-02, &
                  4.358372453d-02, 4.459055816d-02, &
                  4.549162793d-02, 4.628479658d-02, &
                  4.696818282d-02, 4.754016571d-02, &
                  4.799938860d-02, 4.834476223d-02, &
                  4.857546744d-02, 4.869095701d-02]
      
      accum = 0.0d0

      sse=zero                                                     
      dde=zero
      test_S = 0.0d0
                                                                   
      do i = 1, nts
      li=isphe(i)
      nv=30*(i-1)
!                                                                  
! --   LOOP ON INTEGRATION POINTS
!                                                                  
          ssep=0.0d+00
          ddep=0.0d+00
          test_Sp = 0.0d0
                                                                   
          xz=(xtscor(i)-xe(li))/re(li)
          yz=(ytscor(i)-ye(li))/re(li)
          zz=(ztscor(i)-ze(li))/re(li)
          azimuth = atan((ytscor(i) - ye(li))/(xtscor(i)-xe(li)))
          rmin=0.99d+00
          if (abs(xz).le.rmin) then
            rmin=abs(xz)
            xx=0.0d+00
            yx=-zz/sqrt(1.0d+00-xz*xz)
            zx=yz/sqrt(1.0d+00-xz*xz)
          end if
          if (abs(yz).le.rmin) then
            rmin=abs(yz)
            xx=zz/sqrt(1.0d+00-yz*yz)
            yx=0.0d+00
            zx=-xz/sqrt(1.0d+00-yz*yz)
          end if
          if (abs(zz).le.rmin) then
            xx=yz/sqrt(1.0d+00-zz*zz)
            yx=-xz/sqrt(1.0d+00-zz*zz)
            zx=0.0d+00
          end if
          xy=yz*zx-yx*zz
          yy=zz*xx-zx*xz
          zy=xz*yx-xx*yz
                                                                   
                                                                   
          ! Clean-up heap-crap        
          theta = 0.0d0
          phi   = 0.0d0
          numb  = 0
          phinumb = 0.0d0
                                                                   
          do k=1,nvert(i)
            vecx=vert(i,k,1)-xe(li)
            vecy=vert(i,k,2)-ye(li)
            vecz=vert(i,k,3)-ze(li)
            dc=(vecx*xz+vecy*yz+vecz*zz)/re(li)
            if(dc.ge.1.0d+00)dc=1.0d+00
            if(dc.le.-1.0d+00)dc=-1.0d+00
            theta(k)=acos(dc)
            dcs=(vecx*xx+vecy*yx+vecz*zx)/(re(li)*sin(theta(k)))
            if(dcs.ge.1.0d+00)dcs=1.0d+00
            if(dcs.le.-1.0d+00)dcs=-1.0d+00
            phi(k)=acos(dcs)
            snphi=(vecx*xy+vecy*yy+vecz*zy)/(re(li)*sin(theta(k)))
            if (snphi.le.0.0d+00) phi(k)=tpi-phi(k)
          enddo
         
          do k=2,nvert(i)
          phi(k)=phi(k)-phi(1)
          if (phi(k).lt.0.0d+00) phi(k)=tpi+phi(k)
          enddo                                                                 
                                                                                
          xx=xx*cos(phi(1))+xy*sin(phi(1))
          yx=yx*cos(phi(1))+yy*sin(phi(1))
          zx=zx*cos(phi(1))+zy*sin(phi(1))
          xy=yz*zx-yx*zz
          yy=zz*xx-zx*xz
          zy=xz*yx-xx*yz
          
          phi(1)=0.0d+00
          numb(1)=1
          numb(2)=2
          phinumb(1)=phi(1)
          phinumb(2)=phi(2)
          do 210 k=3,nvert(i)
            do l=2,k-1
              if (phi(k).lt.phinumb(l)) then
                do m=1,k-l
                  numb(k-m+1)=numb(k-m)
                  phinumb(k-m+1)=phinumb(k-m)
                enddo
                numb(l)=k
                phinumb(l)=phi(k)
                go to 210
              end if
            enddo
            numb(k)=k
            phinumb(k)=phi(k)
  210     continue
          numb(nvert(i)+1)=numb(1)
          phinumb(nvert(i)+1)=tpi
!                                                                               
! -- LOOP ON GAUSS POINTS
!                                                                               
          do nedge=1,nvert(i)
            aph=(phinumb(nedge+1)-phinumb(nedge))/2.0d+00
            bph=(phinumb(nedge+1)+phinumb(nedge))/2.0d+00
            tha=theta(numb(nedge))
            thb=theta(numb(nedge+1))
            pha=phinumb(nedge)
            phb=phinumb(nedge+1)
            ocx=(centr(i,nedge,1)-xe(li))/re(li)
            ocy=(centr(i,nedge,2)-ye(li))/re(li)
            ocz=(centr(i,nedge,3)-ze(li))/re(li)
            rcc=ocx**2+ocy**2+ocz**2
            
            do nph=1,ng64pts
            do nphsgn=0,1
             ph=(2*nphsgn-1)*aph*xgp64pts(nph)+bph
             cosph=cos(ph)
             sinph=sin(ph)
            if (rcc.lt.1.0d-07) then
             cotgthmax=(sin(ph-pha)/tan(thb)+sin(phb-ph)/tan(tha))/sin(phb-pha)
             thmax=atan(1.0d+00/cotgthmax)
            else
             xcc=xx*ocx+yx*ocy+zx*ocz
             ycc=xy*ocx+yy*ocy+zy*ocz
             zcc=xz*ocx+yz*ocy+zz*ocz
             aa=(xcc*cosph+ycc*sinph)**2+zcc**2
             bb=-zcc*rcc
             cc=rcc**2-(xcc*cosph+ycc*sinph)**2
             ds=bb**2-aa*cc
             if(ds.lt.zero)ds=zero
             cs=(-bb+sqrt(ds))/aa
             if(cs.gt.1.0d+00)cs=1.0d+00
             if(cs.lt.-1.0d+00)cs=1.0d+00
           thmax=acos(cs)
          end if
           if(thmax.lt.1.0d-08) go to 33
           ath=thmax/two
           ssepp=0.0d+00
           ddepp=0.0d+00
           test_Spp = 0.0d0
           
           do nth=1,ng16pts ! Inner 16-points rule
           do nthsgn=0,1
            th=(2*nthsgn-1)*ath*xgp16pts(nth)+ath
            costh=cos(th)
            sinth=sin(th)
            vx=xx*sinth*cosph+xy*sinth*sinph+xz*(costh-1.0d+00)
            vy=yx*sinth*cosph+yy*sinth*sinph+yz*(costh-1.0d+00)
            vz=zx*sinth*cosph+zy*sinth*sinph+zz*(costh-1.0d+00)
            rth=sqrt(2*(1-costh))
            ssepp=ssepp+(re(li)*exp(-kappa*re(li)*rth)/(rth*eps))*sinth*ath*wgp16pts(nth)
            ddepp=ddepp-(rth**2*exp(-kappa*re(li)*rth)*(1.0d+00+kappa*re(li)*rth)/(2*rth**3))*sinth*ath*wgp16pts(nth)
           enddo ! Close loop on nthsgn (n-theta-sign)
           enddo ! Close loop on nth (n-theta)
           ssep=ssep+ssepp*aph*wgp64pts(nph)
           ddep=ddep+ddepp*aph*wgp64pts(nph)
           
  33       continue
          enddo ! Close loop on nphsgn (n-phi-sign)
          enddo ! Close loop on nph (n-phi)
        enddo ! Close loop on nedge
        sse=ssep
        dde=ddep
      end do

      end subroutine diagonal_ionic
      
      end module diagonal
