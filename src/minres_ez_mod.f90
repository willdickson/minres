module minres_ez_mod

    use  minres_data_mod, only : dp
    use  minres_mod,      only : minres_solver

    implicit none
    private


    ! Default parameters
    integer,  parameter :: DFLT_ITNLIM = 1000
    integer,  parameter :: DFLT_NOUT   = 0
    logical,  parameter :: DFLT_PRECON = .false.
    logical,  parameter :: DFLT_CHECKA = .false.
    real(dp), parameter :: DFLT_RTOL   = 1.0e-12_dp


    type, public :: minres_ez_t
        integer  :: itnlim  = DFLT_ITNLIM ! upper limit on the number of iterations.
        integer  :: nout    = DFLT_NOUT   ! a file number for iter. summarry (if nout > 0) 
        logical  :: precon  = DFLT_PRECON ! whether or not to invoke preconditioning 
        logical  :: checka  = DFLT_CHECKA ! whether or not to check if matrix A is symmetric 
        real(dp) :: rtol    = DFLT_RTOL   ! user-specified residual tolerance 
    contains
        private
        procedure, public :: solve => minres_ez_solve
        procedure, public :: print => minres_ez_print
    end type minres_ez_t


    interface minres_ez_t
        procedure :: minres_ez_constructor
    end interface minres_ez_t


    ! Information from solver.  See minresModule.f90 for further details.
    type, public :: minres_info_t   
        integer  :: istop  = 0      ! integer giving the reason for termination
        integer  :: itn    = 0      ! number of iterations performed
        real(dp) :: anorm  = 0.0_dp ! estimate of the norm of the matrix operator
        real(dp) :: acond  = 0.0_dp ! estimate of the condition of Abar 
        real(dp) :: rnorm  = 0.0_dp ! estimate of norm of residual vector
        real(dp) :: ynorm  = 0.0_dp ! estimate of the norm of xbar
    contains
        private
        !procedure, public :: print => minres_info_print
    end type minres_info_t          




contains

    function minres_ez_constructor(itnlim, nout, precon, checka, rtol) result(minres_ez)
        integer,  optional, intent(in) :: itnlim
        integer,  optional, intent(in) :: nout
        integer,  optional, intent(in) :: precon
        integer,  optional, intent(in) :: checka
        real(dp), optional, intent(in) :: rtol
        type(minres_ez_t) :: minres_ez
        if (present(itnlim)) then
            minres_ez % itnlim = itnlim
        end if

        if (present(nout)) then
            minres_ez % nout = nout
        end if

        if (present(precon)) then 
            minres_ez % precon = precon
        end if

        if (present(checka)) then
            minres_ez % checka = checka
        end if

        if (present(rtol)) then 
            minres_ez % rtol = rtol
        end if
    end function minres_ez_constructor

    
    subroutine minres_ez_solve(this, irow, icol, a, b, num, x, info)
        class(minres_ez_t), intent(in)   :: this
        integer, intent(in)              :: irow(:) ! row indices for nonzero items
        integer, intent(in)              :: icol(:) ! col indices for nonzero items
        real(dp), intent(in)             :: a(:)    ! values for nonzero items
        real(dp), intent(in)             :: b(:)    ! the rhs vector b
        integer, intent(in), optional    :: num     ! the number of nonzero items
        real(dp), intent(out)            :: x       ! the computed solution 
        type(minres_info_t), intent(out) :: info    ! information regarding solution

        integer  :: num_nz ! number of nonzero elements in a

        if (present(num)) then
            num_nz = num
        else
            num_nz = size(a)
        end if


    contains

        ! y = A*x
        subroutine aprod(n, x, y)
            integer, intent(in)   :: n
            real(dp), intent(in)  :: x(n)
            real(dp), intent(out) :: y(n)
            integer :: i
            do i  = 1, num_nz
                y(irow(i)) = y(irow(i)) + a(i)*x(icol(i))
            end do
        end subroutine aprod
        
       ! Solve M*y = x
       subroutine msolve(n,x,y)   
         integer,  intent(in)    :: n
         real(dp), intent(in)    :: x(n)
         real(dp), intent(out)   :: y(n)
       end subroutine Msolve

    end subroutine minres_ez_solve 


    subroutine minres_ez_print(this)
        class(minres_ez_t), intent(in) :: this
        print *, 'itnlim = ', this % itnlim
        print *, 'nout   = ', this % nout  
        print *, 'precon = ', this % precon
        print *, 'checka = ', this % checka
        print *, 'rtol   = ', this % rtol  
    end subroutine minres_ez_print


end module minres_ez_mod
