#if defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_HPUX) || defined (SYS_SUN) || defined (SYS_NEC) || defined (SYS_HAL) || defined (SYS_LINUX)
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = IEOR(I,J)
#elif defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = SHIFTL(I,J)
      IBTSHR(I,J) = SHIFTR(I,J)
      IBTXOR(I,J) = XOR(I,J)
#else
      You must define IBTFUN in comdeck file for this computer.
#endif
