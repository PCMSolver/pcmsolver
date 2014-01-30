    module pedra_ibtfun
    
    implicit none

    public ibtand
    public ibtor
    public ibtshl
    public ibtshr
    public ibtxor

    private

    contains

    function ibtand(i, j)
              
    integer :: i, j
    integer ibtand
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
    ibtand = and(i, j)
#else
    ibtand = iand(i, j)
#endif

    end function
      
    function ibtor(i, j)
              
    integer :: i, j
    integer ibtor
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
    ibtor = or(i, j)
#else
    ibtor = ior(i, j)
#endif

    end function
      
    function ibtshl(i, j)
              
    integer :: i, j
    integer ibtshl
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
    ibtshl = shiftl(i, j)
#else
    ibtshl = ishft(i, j)
#endif

    end function

    function ibtshr(i, j)
              
    integer :: i, j
    integer ibtshr
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
    ibtshr = shiftr(i, j)
#else
    ibtshr = ishft(i, j)
#endif

    end function
      
    function ibtxor(i, j)
                
    integer :: i, j
    integer ibtxor
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
    ibtxor = xor(i, j)
#else
    ibtxor = ieor(i, j)
#endif

    end function
      
    end module pedra_ibtfun
