!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File minres_data_mod.f90
!
! Defines real(kind=dp) and a few constants for use in other modules.
!
! 14 Oct 2007: First version implemented after realizing -r8 is not
!              a standard compiler option.
! 15 Oct 2007: Temporarily used real(8) everywhere.
! 16 Oct 2007: Found that we need
!                 use minres_data_mod
!              at the beginning of modules AND inside interfaces.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module minres_data_mod

    use, intrinsic :: iso_fortran_env, only : dp=>real64
    implicit none
    private

    public :: dp
    public :: one
    public :: zero

    real(dp), parameter :: zero = 0.0_dp 
    real(dp), parameter :: one  = 1.0_dp

    !intrinsic                        ::      selected_real_kind
    !integer,       parameter, public :: dp = selected_real_kind(15)
    !real(kind=dp), parameter, public :: zero = 0.0_dp, one = 1.0_dp

end module minres_data_mod
