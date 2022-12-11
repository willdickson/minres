program minres_ez_example_1
    ! 
    ! Example program demonstrating the use of minres_ez_t to solve 
    ! linear systems of the form Ax = b and (A - shift*I)x = b where
    ! A is symmetric.
    ! 
    use, intrinsic :: iso_fortran_env, only : dp=>real64
    use minres, only : minres_ez_t
    use minres, only : minres_info_t
    implicit none

    integer,  parameter :: n         = 3                               ! size of matrix A (symmetric)
    integer,  parameter :: nnz       = 7                               ! number of nonzero elem. in A 
    integer,  parameter :: irow(nnz) = [1, 2, 3, 2, 3, 1, 2]           ! row indices nonzero elem. of A
    integer,  parameter :: icol(nnz) = [1, 2, 3, 1, 2, 2, 3]           ! col indices nonzero elem. of A
    real(dp), parameter :: a(nnz)    = real([3, 2, 1, 4, 2, 4, 2], dp) ! nonzero elem. of A
    real(dp), parameter :: b(n)      = real([1, 1, 1], dp)             ! the rhs vector b

    type(minres_ez_t)   :: minres_ez   ! the main ez solver class
    type(minres_info_t) :: minres_info ! information on solution 
    real(dp)            :: x(n)        ! the solution vector x

    ! Set some options
    minres_ez % rtol   = 1.0e-13
    minres_ez % checka = .true.  

    ! display settings
    call minres_ez % print()

    ! Solve linear system
    call minres_ez % solve(irow, icol, a, b, x, minres_info)

    ! Print info on result
    call minres_info % print()

end program minres_ez_example_1
