program main

    use, intrinsic :: iso_fortran_env, only : dp=>real64
    use minres, only : minres_ez_t
    use minres, only : minres_info_t

    implicit none


    integer,  parameter    :: n = 5
    integer,  parameter    :: sz = n*3 - 2
    !integer,  parameter    :: sz = n
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

    ! Create A matrix
    cnt = 0
    do i = 1, n 
        cnt = cnt + 1
        a(cnt) = 2.0
        irow(cnt) = i
        icol(cnt) = i
    end do
    do i = 1, n-1
        cnt = cnt + 1
        a(cnt) = 0.5
        irow(cnt) = i
        icol(cnt) = i+1  
        cnt = cnt + 1
        a(cnt) = 0.5
        irow(cnt) = i+1
        icol(cnt) = i
    end do 

    print *, ''
    print *, 'matrix A'
    print *, '-----------------------------------'
    do i = 1, size(a)
        print *, irow(i), icol(i), a(i)
    end do
    print *, ''

    ! Create x vector
    do i = 1, n
        x(i) = 1.0
    end do

    ! Create vector b
    call sparse_mult(irow, icol, a, x, b)

    print *, 'vector x'
    print *, '----------------------------------'
    do i = 1, size(x)
        print *, i, x(i)
    end do
    print *, ''

    print *, 'vector b'
    print *, '----------------------------------'
    do i = 1, size(b)
        print *, i, b(i)
    end do
    print *, ''


    call minres_ez % print
    print *, ''

    call minres_ez % solve(irow, icol, a, b, x, minres_info)

    call minres_info % print
    print *, ''

    call sparse_mult(irow, icol, a, x, y)

    print *, ''
    do i = 1, size(x)
        print *, i, b(i), y(i) 

    end do

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
