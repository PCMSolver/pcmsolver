#if !defined(__CVERSION__)
C
C     Parameters NINFI must be updated after changes (for parallelization)
C
C     NOTE: Reals and logicals should appear at the end.
C
      INTEGER NINFI
      PARAMETER (NINFI = 333)
      INTEGER MULD2H, NRHF,NROHF,NVIR, NFRO,
     *       NISH,NASH,NSSH,NOCC,NORB,NBAS,
     *       NNORB,NNBAS, N2ORB,N2BAS,
     *       IISH,IASH,ISSH,IOCC,IORB,IBAS,
     *       IIISH,IIASH,IIORB,IIBAS,I2ORB,I2BAS,
     *       ICMO, NSYM,
     *       NISHT,NASHT,NSSHT,NOCCT,NORBT,NBAST,NCMOT,NRHFT,NVIRT,
     *       N2ISHX,NNASHX,N2ASHX,NNASHY,NNOCCX,N2OCCX,
     *       NNORBT,NNORBX,N2ORBT,N2ORBX,NNBAST,N2BAST,NNBASX,N2BASX,
     *       NNRHFT,NNRHFX,N2RHFT,N2RHFX,NNVIRT,NNVIRX,N2VIRT,N2VIRX,
     *       NAS1,NAS2,NAS3,NNOCCT,N2OCCT,
     *       NAS1T,NAS2T,NAS3T
      COMMON /INFORB/ MULD2H(8,8), NRHF(8),NROHF(8),NVIR(8), NFRO(8),
     *       NISH(8),NASH(8),NSSH(8),NOCC(8),NORB(8),NBAS(8),
     *       NNORB(8),NNBAS(8), N2ORB(8),N2BAS(8),
     *       IISH(8),IASH(8),ISSH(8),IOCC(8),IORB(8),IBAS(8),
     *       IIISH(8),IIASH(8),IIORB(8),IIBAS(8),I2ORB(8),I2BAS(8),
     *       ICMO(8), NSYM,
     *       NISHT,NASHT,NSSHT,NOCCT,NORBT,NBAST,NCMOT,NRHFT,NVIRT,
     *       N2ISHX,NNASHX,N2ASHX,NNASHY,NNOCCX,N2OCCX,
     *       NNORBT,NNORBX,N2ORBT,N2ORBX,NNBAST,N2BAST,NNBASX,N2BASX,
     *       NNRHFT,NNRHFX,N2RHFT,N2RHFX,NNVIRT,NNVIRX,N2VIRT,N2VIRX,
     *       NAS1(8),NAS2(8),NAS3(8),NNOCCT,N2OCCT,
     *       NAS1T,NAS2T,NAS3T
C     MXSSYM = maximum number of "super symmetries"
      INTEGER MXSSYM
      PARAMETER ( MXSSYM = 100 )
      INTEGER NSSYM, NORBSS, IORBSS,
     *       NINFSS, MXDGSS
      COMMON /INFOSS/ NSSYM, NORBSS(MXSSYM), IORBSS(MXSSYM),
     *       NINFSS(MXSSYM,3), MXDGSS
#else
/*     MXSSYM = maximum number of "super symmetries" */
#define MXSSYM 100 
extern struct common_inforb {
    integer muld2h[8][8], nrhf[8],nrohf[8],nvir[8], nfro[8],
	nish[8],nash[8],nssh[8],nocc[8],norb[8],nbas[8],
	nnorb[8],nnbas[8], n2orb[8],n2bas[8],
	iish[8],iash[8],issh[8],iocc[8],iorb[8],ibas[8],
	iiish[8],iiash[8],iiorb[8],iibas[8],i2orb[8],i2bas[8],
	icmo[8], nsym,
	nisht,nasht,nssht,nocct,norbt,nbast,ncmot,nrhft,nvirt,
	n2ishx,nnashx,n2ashx,nnashy,nnoccx,n2occx,
	nnorbt,nnorbx,n2orbt,n2orbx,nnbast,n2bast,nnbasx,n2basx,
	nnrhft,nnrhfx,n2rhft,n2rhfx,nnvirt,nnvirx,n2virt,n2virx,
	nas1[8],nas2[8],nas3[8],nnocct,n2occt,
	nas1t,nas2t,nas3t;
} inforb_;
extern struct common_infoss {
    /* NOTICE that the index range is inverted for ninfss */
    integer nssym, norbss[MXSSYM], iorbss[MXSSYM],
	ninfss[3][MXSSYM], mxdgss;
} infoss_;
#endif
