      SUBROUTINE DIAGAN(I,SSE,DDE,VERT,CENTR,XCTS,YCTS,ZCTS,AS,
     *                  NVERT,XE,YE,ZE,RE,ISPHE)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (MXVRT=20)
C
      DIMENSION VERT(NTS*MXVRT*3),CENTR(NTS*MXVRT*3)
      DIMENSION NUMB(MXVRT),THETA(MXVRT),PHI(MXVRT),PHINUMB(MXVRT),
     *          XGP16PTS(8),WGP16PTS(8), XGP64PTS(32),WGP64PTS(32),
     *          XCTS(*),YCTS(*),ZCTS(*),AS(*),NVERT(*),XE(*),YE(*),
     *          ZE(*),RE(*),ISPHE(*)
C
      COMMON /ANIDAT/ EPSM,EPSXX,EPSXY,EPSXZ,EPSYY,EPSYZ,EPSZZ,
     *                EPSM1XX,EPSM1XY,EPSM1XZ,EPSM1YY,EPSM1YZ,EPSM1ZZ,
     *                EPS1,EPS2,EPS3,ROT(3,3)
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
C
      DATA ZERO,TWO/0.0D+00,2.0D+00/
      DATA TPI,FPI/6.28318531D+00,12.56637061D+00/
C
      DATA NG16PTS/8/
      DATA XGP16PTS/9.894009350D-01,9.445750231D-01,
     *              8.656312024D-01,7.554044084D-01,
     *              6.178762444D-01,4.580167777D-01,
     *              2.816035508D-01,9.501250984D-02/
      DATA WGP16PTS/2.715245941D-02,6.225352394D-02,
     *              9.515851168D-02,1.246289713D-01,
     *              1.495959888D-01,1.691565194D-01,
     *              1.826034150D-01,1.894506105D-01/
C
      DATA NG64PTS/32/
      DATA XGP64PTS/9.993050417D-01,9.963401168D-01,
     *              9.910133715D-01,9.833362539D-01,
     *              9.733268278D-01,9.610087996D-01,
     *              9.464113748D-01,9.295691721D-01,
     *              9.105221371D-01,8.893154460D-01,
     *              8.659993982D-01,8.406292963D-01,
     *              8.132653151D-01,7.839723589D-01,
     *              7.528199073D-01,7.198818502D-01,
     *              6.852363131D-01,6.489654713D-01,
     *              6.111553552D-01,5.718956462D-01,
     *              5.312794640D-01,4.894031457D-01,
     *              4.463660173D-01,4.022701580D-01,
     *              3.572201583D-01,3.113228720D-01,
     *              2.646871622D-01,2.174236437D-01,
     *              1.696444204D-01,1.214628193D-01,
     *              7.299312179D-02,2.435029266D-02/
      DATA WGP64PTS/1.783280722D-03,4.147033261D-03,
     *              6.504457969D-03,8.846759826D-03,
     *              1.116813946D-02,1.346304790D-02,
     *              1.572603048D-02,1.795171578D-02,
     *              2.013482315D-02,2.227017381D-02,
     *              2.435270257D-02,2.637746972D-02,
     *              2.833967261D-02,3.023465707D-02,
     *              3.205792835D-02,3.380516184D-02,
     *              3.547221326D-02,3.705512854D-02,
     *              3.855015318D-02,3.995374113D-02,
     *              4.126256324D-02,4.247351512D-02,
     *              4.358372453D-02,4.459055816D-02,
     *              4.549162793D-02,4.628479658D-02,
     *              4.696818282D-02,4.754016571D-02,
     *              4.799938860D-02,4.834476223D-02,
     *              4.857546744D-02,4.869095701D-02/
C
C
      SSE=ZERO
      DDE=ZERO
      LI=ISPHE(I)
      NV=30*(I-1)
C
C -- LOOP ON INTEGRATION POINTS
C
        SSEP=0.0D+00
        DDEP=0.0D+00
        XZ=(XCTS(I)-XE(LI))/RE(LI)
        YZ=(YCTS(I)-YE(LI))/RE(LI)
        ZZ=(ZCTS(I)-ZE(LI))/RE(LI)
        RMIN=0.99D+00
        IF (ABS(XZ).LE.RMIN) THEN
          RMIN=ABS(XZ)
          XX=0.0D+00
          YX=-ZZ/SQRT(1.0D+00-XZ*XZ)
          ZX=YZ/SQRT(1.0D+00-XZ*XZ)
        END IF
        IF (ABS(YZ).LE.RMIN) THEN
          RMIN=ABS(YZ)
          XX=ZZ/SQRT(1.0D+00-YZ*YZ)
          YX=0.0D+00
          ZX=-XZ/SQRT(1.0D+00-YZ*YZ)
        END IF
        IF (ABS(ZZ).LE.RMIN) THEN
          XX=YZ/SQRT(1.0D+00-ZZ*ZZ)
          YX=-XZ/SQRT(1.0D+00-ZZ*ZZ)
          ZX=0.0D+00
        END IF
        XY=YZ*ZX-YX*ZZ
        YY=ZZ*XX-ZX*XZ
        ZY=XZ*YX-XX*YZ
C
        DO K=1,NVERT(I)
          VECX=VERT(NV+3*(K-1)+1)-XE(LI)
          VECY=VERT(NV+3*(K-1)+2)-YE(LI)
          VECZ=VERT(NV+3*(K-1)+3)-ZE(LI)
          DC=(VECX*XZ+VECY*YZ+VECZ*ZZ)/RE(LI)
          IF(DC.GE.1.0D+00)DC=1.0D+00
          IF(DC.LE.-1.0D+00)DC=-1.0D+00
          THETA(K)=ACOS(DC)
          DCS=(VECX*XX+VECY*YX+VECZ*ZX)
     *        /(RE(LI)*SIN(THETA(K)))
          IF(DCS.GE.1.0D+00)DCS=1.0D+00
          IF(DCS.LE.-1.0D+00)DCS=-1.0D+00
          PHI(K)=ACOS(DCS)
          SNPHI=(VECX*XY+VECY*YY+VECZ*ZY)/(RE(LI)*SIN(THETA(K)))
          IF (SNPHI.LE.0.0D+00) PHI(K)=TPI-PHI(K)
        ENDDO
C
        DO K=2,NVERT(I)
          PHI(K)=PHI(K)-PHI(1)
          IF (PHI(K).LT.0.0D+00) PHI(K)=TPI+PHI(K)
        ENDDO
C
        XX=XX*COS(PHI(1))+XY*SIN(PHI(1))
        YX=YX*COS(PHI(1))+YY*SIN(PHI(1))
        ZX=ZX*COS(PHI(1))+ZY*SIN(PHI(1))
        XY=YZ*ZX-YX*ZZ
        YY=ZZ*XX-ZX*XZ
        ZY=XZ*YX-XX*YZ
C
        PHI(1)=0.0D+00
        NUMB(1)=1
        NUMB(2)=2
        PHINUMB(1)=PHI(1)
        PHINUMB(2)=PHI(2)
C
        DO 210 K=3,NVERT(I)
          DO L=2,K-1
            IF (PHI(K).LT.PHINUMB(L)) THEN
              DO M=1,K-L
                NUMB(K-M+1)=NUMB(K-M)
                PHINUMB(K-M+1)=PHINUMB(K-M)
              ENDDO
              NUMB(L)=K
              PHINUMB(L)=PHI(K)
              GO TO 210
            END IF
          ENDDO
          NUMB(K)=K
          PHINUMB(K)=PHI(K)
  210   CONTINUE
        NUMB(NVERT(I)+1)=NUMB(1)
        PHINUMB(NVERT(I)+1)=TPI
C
C -- LOOP ON GAUSS POINTS
C
        DO NEDGE=1,NVERT(I)
C
          APH=(PHINUMB(NEDGE+1)-PHINUMB(NEDGE))/2.0D+00
          BPH=(PHINUMB(NEDGE+1)+PHINUMB(NEDGE))/2.0D+00
          THA=THETA(NUMB(NEDGE))
          THB=THETA(NUMB(NEDGE+1))
          PHA=PHINUMB(NEDGE)
          PHB=PHINUMB(NEDGE+1)
          OCX=(CENTR(NV+3*(NEDGE-1)+1)-XE(LI))/RE(LI)
          OCY=(CENTR(NV+3*(NEDGE-1)+2)-YE(LI))/RE(LI)
          OCZ=(CENTR(NV+3*(NEDGE-1)+3)-ZE(LI))/RE(LI)
          RCC=OCX**2+OCY**2+OCZ**2
C
C
          DO NPH=1,NG64PTS
          DO NPHSGN=0,1
C
           PH=(2*NPHSGN-1)*APH*XGP64PTS(NPH)+BPH
           COSPH=COS(PH)
           SINPH=SIN(PH)
C
          IF (RCC.LT.1.0D-07) THEN
           COTGTHMAX=(SIN(PH-PHA)/TAN(THB)+SIN(PHB-PH)/TAN(THA))
     *               /SIN(PHB-PHA)
           THMAX=ATAN(1.0D+00/COTGTHMAX)
          ELSE
           XCC=XX*OCX+YX*OCY+ZX*OCZ
           YCC=XY*OCX+YY*OCY+ZY*OCZ
           ZCC=XZ*OCX+YZ*OCY+ZZ*OCZ
           AA=(XCC*COSPH+YCC*SINPH)**2+ZCC**2
           BB=-ZCC*RCC
           CC=RCC**2-(XCC*COSPH+YCC*SINPH)**2
           DS=BB**2-AA*CC
           IF(DS.LT.ZERO)DS=ZERO
           CS=(-BB+SQRT(DS))/AA
           IF(CS.GT.1.0D+00)CS=1.0D+00
           IF(CS.LT.-1.0D+00)CS=1.0D+00
           THMAX=ACOS(CS)
          END IF
           IF(THMAX.LT.1.0D-08) GO TO 33
C
           ATH=THMAX/TWO
C
           SSEPP=0.0D+00
           DDEPP=0.0D+00
C
           DO NTH=1,NG16PTS
           DO NTHSGN=0,1
C
            TH=(2*NTHSGN-1)*ATH*XGP16PTS(NTH)+ATH
            COSTH=COS(TH)
            SINTH=SIN(TH)
            VX=XX*SINTH*COSPH+XY*SINTH*SINPH+XZ*(COSTH-1.0D+00)
            VY=YX*SINTH*COSPH+YY*SINTH*SINPH+YZ*(COSTH-1.0D+00)
            VZ=ZX*SINTH*COSPH+ZY*SINTH*SINPH+ZZ*(COSTH-1.0D+00)
            RTH=SQRT(2*(1-COSTH))
            RTHEPS=SQRT(EPSM1XX*VX*VX+2*EPSM1XY*VX*VY+2*EPSM1XZ*VX*VZ
     *                               +  EPSM1YY*VY*VY+2*EPSM1YZ*VY*VZ
     *                                               +  EPSM1ZZ*VZ*VZ)
            SSEPP=SSEPP+(RE(LI)/(FPI*RTHEPS))*SINTH*ATH*WGP16PTS(NTH)
            DDEPP=DDEPP-(RTH**2/(2*FPI*RTHEPS**3))*SINTH
     *                   *ATH*WGP16PTS(NTH)
C
           ENDDO
           ENDDO
C
           SSEP=SSEP+SSEPP*APH*WGP64PTS(NPH)
           DDEP=DDEP+DDEPP*APH*WGP64PTS(NPH)
  33       CONTINUE
C
          ENDDO
          ENDDO
C
        ENDDO
C
        SSE=SSEP*AS(I)
        DDE=DDEP*AS(I)
C
      RETURN
      END
C*MODULE PCMIEF  *DECK DIAGION
      SUBROUTINE DIAGION(I,SSE,DDE,VERT,CENTR,XCTS,YCTS,ZCTS,AS,
     *                   NVERT,XE,YE,ZE,RE,ISPHE)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (MXVRT=20)
C
      DIMENSION VERT(NTS*MXVRT*3),CENTR(NTS*MXVRT*3)
      DIMENSION NUMB(MXVRT),THETA(MXVRT),PHI(MXVRT),PHINUMB(MXVRT),
     *          XGP16PTS(8),WGP16PTS(8),XGP64PTS(32),WGP64PTS(32),
     *          XCTS(*),YCTS(*),ZCTS(*),AS(*),NVERT(*),XE(*),YE(*),
     *          ZE(*),RE(*),ISPHE(*)
C
      COMMON /IONDAT/ EPSI,DK2,DALP
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
C
      DATA ZERO,TWO/0.0D+00,2.0D+00/
      DATA TPI,FPI/6.28318531D+00,12.56637061D+00/
C
      DATA NG16PTS/8/
      DATA XGP16PTS/9.894009350D-01,9.445750231D-01,
     *              8.656312024D-01,7.554044084D-01,
     *              6.178762444D-01,4.580167777D-01,
     *              2.816035508D-01,9.501250984D-02/
      DATA WGP16PTS/2.715245941D-02,6.225352394D-02,
     *              9.515851168D-02,1.246289713D-01,
     *              1.495959888D-01,1.691565194D-01,
     *              1.826034150D-01,1.894506105D-01/
C
      DATA NG64PTS/32/
      DATA XGP64PTS/9.993050417D-01,9.963401168D-01,
     *              9.910133715D-01,9.833362539D-01,
     *              9.733268278D-01,9.610087996D-01,
     *              9.464113748D-01,9.295691721D-01,
     *              9.105221371D-01,8.893154460D-01,
     *              8.659993982D-01,8.406292963D-01,
     *              8.132653151D-01,7.839723589D-01,
     *              7.528199073D-01,7.198818502D-01,
     *              6.852363131D-01,6.489654713D-01,
     *              6.111553552D-01,5.718956462D-01,
     *              5.312794640D-01,4.894031457D-01,
     *              4.463660173D-01,4.022701580D-01,
     *              3.572201583D-01,3.113228720D-01,
     *              2.646871622D-01,2.174236437D-01,
     *              1.696444204D-01,1.214628193D-01,
     *              7.299312179D-02,2.435029266D-02/
      DATA WGP64PTS/1.783280722D-03,4.147033261D-03,
     *              6.504457969D-03,8.846759826D-03,
     *              1.116813946D-02,1.346304790D-02,
     *              1.572603048D-02,1.795171578D-02,
     *              2.013482315D-02,2.227017381D-02,
     *              2.435270257D-02,2.637746972D-02,
     *              2.833967261D-02,3.023465707D-02,
     *              3.205792835D-02,3.380516184D-02,
     *              3.547221326D-02,3.705512854D-02,
     *              3.855015318D-02,3.995374113D-02,
     *              4.126256324D-02,4.247351512D-02,
     *              4.358372453D-02,4.459055816D-02,
     *              4.549162793D-02,4.628479658D-02,
     *              4.696818282D-02,4.754016571D-02,
     *              4.799938860D-02,4.834476223D-02,
     *              4.857546744D-02,4.869095701D-02/
C
      SSE=ZERO
      DDE=ZERO
      LI=ISPHE(I)
      NV=30*(I-1)
C
C -- LOOP ON INTEGRATION POINTS
C
        SSEP=0.0D+00
        DDEP=0.0D+00
        XZ=(XCTS(I)-XE(LI))/RE(LI)
        YZ=(YCTS(I)-YE(LI))/RE(LI)
        ZZ=(ZCTS(I)-ZE(LI))/RE(LI)
        RMIN=0.99D+00
        IF (ABS(XZ).LE.RMIN) THEN
          RMIN=ABS(XZ)
          XX=0.0D+00
          YX=-ZZ/SQRT(1.0D+00-XZ*XZ)
          ZX=YZ/SQRT(1.0D+00-XZ*XZ)
        END IF
        IF (ABS(YZ).LE.RMIN) THEN
          RMIN=ABS(YZ)
          XX=ZZ/SQRT(1.0D+00-YZ*YZ)
          YX=0.0D+00
          ZX=-XZ/SQRT(1.0D+00-YZ*YZ)
        END IF
        IF (ABS(ZZ).LE.RMIN) THEN
          XX=YZ/SQRT(1.0D+00-ZZ*ZZ)
          YX=-XZ/SQRT(1.0D+00-ZZ*ZZ)
          ZX=0.0D+00
        END IF
        XY=YZ*ZX-YX*ZZ
        YY=ZZ*XX-ZX*XZ
        ZY=XZ*YX-XX*YZ
C
        DO K=1,NVERT(I)
          VECX=VERT(NV+3*(K-1)+1)-XE(LI)
          VECY=VERT(NV+3*(K-1)+2)-YE(LI)
          VECZ=VERT(NV+3*(K-1)+3)-ZE(LI)
          DC=(VECX*XZ+VECY*YZ+VECZ*ZZ)/RE(LI)
          IF(DC.GE.1.0D+00)DC=1.0D+00
          IF(DC.LE.-1.0D+00)DC=-1.0D+00
          THETA(K)=ACOS(DC)
          DCS=(VECX*XX+VECY*YX+VECZ*ZX)
     *        /(RE(LI)*SIN(THETA(K)))
          IF(DCS.GE.1.0D+00)DCS=1.0D+00
          IF(DCS.LE.-1.0D+00)DCS=-1.0D+00
          PHI(K)=ACOS(DCS)
          SNPHI=(VECX*XY+VECY*YY+VECZ*ZY)/(RE(LI)*SIN(THETA(K)))
          IF (SNPHI.LE.0.0D+00) PHI(K)=TPI-PHI(K)
        ENDDO
C
        DO K=2,NVERT(I)
          PHI(K)=PHI(K)-PHI(1)
          IF (PHI(K).LT.0.0D+00) PHI(K)=TPI+PHI(K)
        ENDDO
C
        XX=XX*COS(PHI(1))+XY*SIN(PHI(1))
        YX=YX*COS(PHI(1))+YY*SIN(PHI(1))
        ZX=ZX*COS(PHI(1))+ZY*SIN(PHI(1))
        XY=YZ*ZX-YX*ZZ
        YY=ZZ*XX-ZX*XZ
        ZY=XZ*YX-XX*YZ
C
        PHI(1)=0.0D+00
        NUMB(1)=1
        NUMB(2)=2
        PHINUMB(1)=PHI(1)
        PHINUMB(2)=PHI(2)
C
        DO 210 K=3,NVERT(I)
          DO L=2,K-1
            IF (PHI(K).LT.PHINUMB(L)) THEN
              DO M=1,K-L
                NUMB(K-M+1)=NUMB(K-M)
                PHINUMB(K-M+1)=PHINUMB(K-M)
              ENDDO
              NUMB(L)=K
              PHINUMB(L)=PHI(K)
              GO TO 210
            END IF
          ENDDO
          NUMB(K)=K
          PHINUMB(K)=PHI(K)
  210   CONTINUE
        NUMB(NVERT(I)+1)=NUMB(1)
        PHINUMB(NVERT(I)+1)=TPI
C
C -- LOOP ON GAUSS POINTS
C
        DO NEDGE=1,NVERT(I)
C
          APH=(PHINUMB(NEDGE+1)-PHINUMB(NEDGE))/2.0D+00
          BPH=(PHINUMB(NEDGE+1)+PHINUMB(NEDGE))/2.0D+00
          THA=THETA(NUMB(NEDGE))
          THB=THETA(NUMB(NEDGE+1))
          PHA=PHINUMB(NEDGE)
          PHB=PHINUMB(NEDGE+1)
          OCX=(CENTR(NV+3*(NEDGE-1)+1)-XE(LI))/RE(LI)
          OCY=(CENTR(NV+3*(NEDGE-1)+2)-YE(LI))/RE(LI)
          OCZ=(CENTR(NV+3*(NEDGE-1)+3)-ZE(LI))/RE(LI)
          RCC=OCX**2+OCY**2+OCZ**2
C
C
          DO NPH=1,NG64PTS
          DO NPHSGN=0,1
C
           PH=(2*NPHSGN-1)*APH*XGP64PTS(NPH)+BPH
           COSPH=COS(PH)
           SINPH=SIN(PH)
C
          IF (RCC.LT.1.0D-07) THEN
           COTGTHMAX=(SIN(PH-PHA)/TAN(THB)+SIN(PHB-PH)/TAN(THA))
     *               /SIN(PHB-PHA)
           THMAX=ATAN(1.0D+00/COTGTHMAX)
          ELSE
           XCC=XX*OCX+YX*OCY+ZX*OCZ
           YCC=XY*OCX+YY*OCY+ZY*OCZ
           ZCC=XZ*OCX+YZ*OCY+ZZ*OCZ
           AA=(XCC*COSPH+YCC*SINPH)**2+ZCC**2
           BB=-ZCC*RCC
           CC=RCC**2-(XCC*COSPH+YCC*SINPH)**2
           DS=BB**2-AA*CC
           IF(DS.LT.ZERO)DS=ZERO
           CS=(-BB+SQRT(DS))/AA
           IF(CS.GT.1.0D+00)CS=1.0D+00
           IF(CS.LT.-1.0D+00)CS=1.0D+00
           THMAX=ACOS(CS)
          END IF
           IF(THMAX.LT.1.0D-08) GO TO 33
C
           ATH=THMAX/TWO
C
           SSEPP=0.0D+00
           DDEPP=0.0D+00
C
           DO NTH=1,NG16PTS
           DO NTHSGN=0,1
C
            TH=(2*NTHSGN-1)*ATH*XGP16PTS(NTH)+ATH
            COSTH=COS(TH)
            SINTH=SIN(TH)
            RTH=SQRT(2*(1-COSTH))
            SSEPP=SSEPP+(RE(LI)*EXP(-DALP*RE(LI)*RTH)
     *                   /(FPI*RTH*EPSI))*SINTH*ATH*WGP16PTS(NTH)
            DDEPP=DDEPP-(RTH**2*EXP(-DALP*RE(LI)*RTH)
     *                   *(1.0D+00+DALP*RE(LI)*RTH)
     *                   /(2*FPI*RTH**3))*SINTH*ATH*WGP16PTS(NTH)
C
           ENDDO
           ENDDO
C
           SSEP=SSEP+SSEPP*APH*WGP64PTS(NPH)
           DDEP=DDEP+DDEPP*APH*WGP64PTS(NPH)
  33       CONTINUE
C
          ENDDO
          ENDDO
C
        ENDDO
C
        SSE=SSEP*AS(I)
        DDE=DDEP*AS(I)
C
      RETURN
      END
