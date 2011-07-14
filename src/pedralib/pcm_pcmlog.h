C
C Logical variables for pcm calculations
C
C PCM:    True if we perform a PCM calculation
C OUTFLD: Local field correction for pure liquids
      LOGICAL PCM,OUTFLD,NEWMAT,LOCFLD,NONEQ,NEQRSP,
     $        NPCMIN,OLDCEN,NEWQR
      COMMON /PCMLOG/ PCM,OUTFLD,NEWMAT,LOCFLD,NONEQ,NEQRSP,
     $                NPCMIN,OLDCEN,NEWQR
