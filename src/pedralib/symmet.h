C --- FILE: symmet.h ---
C
C Symmetry information for basis functions, generated in abacus/herrdn.F.
C
C Below: 0 .le. i .and. i .le. MAXREP (irrep index); i1 = i + 1 = irrep no.
C        1 .le. j .and. j .le. 
C
C FMULT(i) =
C PT(i) =
C MAXREP = 2**NSYMOP - 1; NSYM = MAXREP + 1 is number of irreps
C MAXOPR = maximum number of operations necessary to loop over [0..7]
C MULT(i) =
C ISYMAX(q,1) = irrep of q-axis
C ISYMAX(q,2) = irrep of rotation around q-axis
C ISYMAO(,) =
C NPARSU(8)   : offset pointer for symmetry dependent AOs for given irrep
C NAOS(8)     : Number of AO'S for given irrep
C NPARNU(8,8) : offset pointer from non-symmetric operators for given irrep
C ...
C NCOS(8,-1)  : Number of AO's for density fitting for given irrep
C NCOS(8, 0)  : Number of AO's for Huckel for given irrep
C NCOS(8, 1)  : Number of AO's for Large(DIRAC) / Main(LMULBS) for given irrep
C NCOS(8, 2)  : Number of AO's for Small(DIRAC) / Axiliary(LMULBS) for given irrep
C ...
C ICLASS(j) =
C ICNTAO(j) =
C
      REAL*8  FMULT, PT
      INTEGER MAXREP, MAXOPR, MULT, ISYMAX, ISYMAO, NPARSU,
     &        NAOS, NPARNU, IPTSYM, IPTCNT, NCRREP,
     &        IPTCOR, NAXREP, IPTAX, IPTXYZ, IPTNUC, ISOP, 
     &        NROTS, NINVC, NREFL, IXVAL, NCOS, ICLASS, ICNTAO
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &        MAXREP, MAXOPR, MULT(0:7),
     &          ISYMAX(3,2), ISYMAO(MXQN,MXAQN), NPARSU(8),
     &        NAOS(8), NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &          IPTCNT(3*MXCENT,0:7,2), NCRREP(0:7,2),
     &        IPTCOR(3*MXCENT,2), NAXREP(0:7,2), IPTAX(3,2),
     &          IPTXYZ(3,0:7,2), IPTNUC(MXCENT,0:7), ISOP(0:7),
     &        NROTS,NINVC,NREFL,IXVAL(0:7,0:7),NCOS(8,-1:2),
     &          ICLASS(MXCORB), ICNTAO(MXCORB)
C --- end of symmet.h ---
