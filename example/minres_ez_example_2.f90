program minres_ez_example_2 
    ! 
    ! Example program demonstrating the use of minres_ez_t to solve 
    ! linear systems of the form Ax = b and (A - shift*I)x = b where
    ! A is symmetric.
    ! 
    use, intrinsic :: iso_fortran_env, only : dp=>real64
    use minres, only : minres_ez_t
    use minres, only : minres_info_t
    implicit none

    integer,  parameter    :: k = 1000  ! size of matrix the matrix A
    integer,  allocatable  :: irow(:)   ! row indices of nonzero elements
    integer,  allocatable  :: icol(:)   ! col indices of nonzero elements
    real(dp), allocatable  :: a(:)      ! value of nonzero elements
    real(dp), allocatable  :: b(:)      ! the right hand side vector b
    real(dp), allocatable  :: xtrue(:)  ! the true solution vector x
    real(dp), allocatable  :: xsolv(:)  ! the solution vector found by minres
    real(dp)               :: max_err   ! maximum abs error 

    type(minres_ez_t)      :: minres_ez
    type(minres_info_t)    :: minres_info

    ! Print minres_ez configuration
    minres_ez % checka = .true.
    minres_ez % precon = .false.
    call minres_ez % print()

    ! Create system of linear equations Ax=b
    call create_linsys(k, irow, icol, a, b, xtrue)
    allocate(xsolv(k), source=0.0_dp)

    ! Solve system of equations and print solver info
    call minres_ez % solve(irow, icol, a, b, xsolv, minres_info)
    call minres_info % print()

    ! Compute the maximum absolute error between true solution and that
    ! found by the solver.
    max_err = maxval(abs(xtrue - xsolv))
    print *, 'max abs error = ', max_err

contains

    ! Generates an example linear system Ax=b with symmetric positive definite matrix A
    subroutine create_linsys(k, irow, icol, a, b, x)
        integer,  intent(in)                 :: k        ! matrix size
        integer,  allocatable, intent(inout) :: irow(:)  ! row indices   (nonzero elem.)
        integer,  allocatable, intent(inout) :: icol(:)  ! col indices   (nonzero elem.)
        real(dp), allocatable, intent(inout) :: a(:)     ! matrix values (nonzero elem.)
        real(dp), allocatable, intent(inout) :: b(:)     ! right hand side vector
        real(dp), allocatable, intent(inout) :: x(:)     ! true solution of Ax=b 

        ! local variables
        integer   :: ind  ! index for elem. in irow, icol, a
        integer   :: nz   ! number of nonzero elements
        integer   :: i    ! loop index

        ! number of nonzero elements
        nz = 3*k - 2
        
        ! Deallocate is necessary
        if (allocated(irow)) then 
            deallocate(irow)
        end if
        allocate(irow(nz))
        irow = 0

        if (allocated(icol)) then 
            deallocate(icol)
        end if
        allocate(icol(nz))
        icol = 0

        if (allocated(a)) then
            deallocate(a)
        end if
        allocate(a(nz))
        a = 0.0_dp

        if (allocated(b)) then 
            deallocate(b)
        end if
        allocate(b(k))

        if (allocated(x)) then
            deallocate(x)
        end if
        allocate(x(k), source=1.0_dp)

        ! Populate irow, icol and a
        ind = 0
        do i = 1, k
            ! diagonal elements
            ind = ind + 1
            icol(ind) = i
            irow(ind) = i
            a(i) = 1.0_dp
        end do 
        do i = 1, k-1
            ! off diagonal elements
            ind = ind + 1
            icol(ind) = i
            irow(ind) = i+1
            a(i) = 2.0_dp
            ind = ind + 1
            icol(ind) = i+1
            irow(ind) = i
            a(i) = 2.0_dp
        end do

        ! Create vector b
        call sparse_mul(irow, icol, a, x, b)
    end subroutine create_linsys


    subroutine sparse_mul(irow, icol, a, x, y) 
        integer,  intent(in)  :: irow(:)
        integer,  intent(in)  :: icol(:)
        real(dp), intent(in)  :: a(:)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        integer               :: i
        y = 0.0_dp
        do i = 1, size(a)
            y(irow(i)) = y(irow(i)) + a(i)*x(icol(i))
        end do
    end subroutine sparse_mul


end program minres_ez_example_2 
