program main

    use, intrinsic :: iso_fortran_env, only : dp=>real64
    use minres, only : minres_ez_t
    use minres, only : minres_info_t

    implicit none


    integer,  parameter    :: n = 80000 
    integer,  parameter    :: sz = 3*n - 2
    integer,  allocatable  :: irow(:)
    integer,  allocatable  :: icol(:)
    real(dp), allocatable  :: a(:)
    real(dp), allocatable  :: b(:)
    real(dp), allocatable  :: x(:)
    real(dp), allocatable  :: y(:)
    type(minres_ez_t)      :: minres_ez
    type(minres_info_t)    :: minres_info

    integer                :: i, j
    integer                :: cnt
    integer                :: nout
    real(dp)               :: err


    allocate(irow(sz))
    allocate(icol(sz))
    allocate(a(sz))

    allocate(b(n))
    allocate(x(n))
    allocate(y(n))

    a = 0.0_dp
    b = 0.0_dp
    x = 0.0_dp
    y = 0.0_dp
    icol = 0
    irow = 0

    cnt = 0
    do i = 1, n
        cnt = cnt + 1
        icol(cnt) = i
        irow(cnt) = i
        a(i) = 1.0_dp
    end do 
    do i = 1, n-1
        cnt = cnt + 1
        icol(cnt) = i
        irow(cnt) = i+1
        a(i) = 2.0_dp
        cnt = cnt + 1
        icol(cnt) = i+1
        irow(cnt) = i
        a(i) = 2.0_dp
    end do

    do i = 1, n
        x(i) = 1.0_dp
    end do

    ! Create vector b
    call sparse_mult(irow, icol, a, x, b)

    call minres_ez % print
    print *, ''

    call minres_ez % solve(irow, icol, a, b, x, minres_info)

    call minres_info % print
    print *, ''

    call sparse_mult(irow, icol, a, x, y)

    err = 0.0_dp
    do i = 1, size(x)
        if (i==1) then
            err = abs(b(i) - y(i))
        else
            err = max(b(i) - y(i), err)
        end if
    end do
    print *, 'err = ', err

contains

    subroutine sparse_mult(irow, icol, a, x, y) 
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
    end subroutine sparse_mult


end program main
