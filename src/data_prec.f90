module data_prec
        implicit none
#ifdef CPU_DP
        integer, parameter, public :: cputype = KIND(0.0D0)
#else
        integer, parameter, public :: cputype = KIND(0.0)
#endif

end module data_prec
