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
        integer  :: itnlim  = DFLT_ITNLIM ! upper limit on the nnzber of iterations.
        integer  :: nout    = DFLT_NOUT   ! a file nnzber for iter. summarry (if nout > 0) 
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


    ! Information from solver.  See minres_module.f90 for further details.
    type, public :: minres_info_t   
        integer  :: istop  = 0      ! integer giving the reason for termination
        integer  :: itn    = 0      ! nnzber of iterations performed
        real(dp) :: anorm  = 0.0_dp ! estimate of the norm of the matrix operator
        real(dp) :: acond  = 0.0_dp ! estimate of the condition of Abar 
        real(dp) :: rnorm  = 0.0_dp ! estimate of norm of residual vector
        real(dp) :: arnorm = 0.0_dp ! recognize singular systems ||Ar||
        real(dp) :: ynorm  = 0.0_dp ! estimate of the norm of xbar
    contains
        private
        procedure, public :: print => minres_info_print
    end type minres_info_t          


contains

    ! minres_ex_t contructor + methods
    ! -----------------------------------------------------------------------------------

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

    
    subroutine minres_ez_solve(this, irow, icol, a, b, x, info, nnz_in, shift_in)
        class(minres_ez_t), intent(in)   :: this
        integer, intent(in)              :: irow(:)  ! row indices for nonzero items
        integer, intent(in)              :: icol(:)  ! col indices for nonzero items
        real(dp), intent(in)             :: a(:)     ! values for nonzero items
        real(dp), intent(in)             :: b(:)     ! the rhs vector b
        real(dp), intent(out)            :: x(:)     ! the computed solution 
        type(minres_info_t), intent(out) :: info     ! information regarding solution
        integer, intent(in), optional    :: nnz_in   ! number of nonzero items
        real(dp), intent(in), optional   :: shift_in ! shift value (A - shift*I) x = b

        integer                          :: nnz   ! number of nonzero elements in a
        real(dp)                         :: shift ! offset shift value
        real(dp), allocatable            :: m(:) 

        if (present(nnz_in)) then 
            nnz = nnz_in
        else
            nnz = size(a)
        end if

        if (present(shift_in)) then
            shift = shift_in
        else
            shift = 0.0_dp
        end if

        if (this % precon) then
            call create_precon(irow, icol, a, nnz, size(b), m)
        end if

        call minres_solver(      &
                size(b),         & 
                aprod,           & 
                msolve,          & 
                b,               & 
                shift,           & 
                this % checka,   &
                this % precon,   &
                x,               &
                this % itnlim,   &  
                this % nout,     &
                this % rtol,     &
                info % istop,    &
                info % itn,      &
                info % anorm,    &
                info % acond,    &
                info % rnorm,    &
                info % arnorm,   &
                info % ynorm     &
                )

    contains

        ! y = A*x
        subroutine aprod(n, x, y)
            integer, intent(in)   :: n
            real(dp), intent(in)  :: x(n)
            real(dp), intent(out) :: y(n)
            integer               :: i
            y = 0.0_dp
            do i  = 1, nnz 
                y(irow(i)) = y(irow(i)) + a(i)*x(icol(i))
            end do
        end subroutine aprod
        
        ! Solve M*y = x
        subroutine msolve(n,x,y)   
          integer,  intent(in)    :: n
          real(dp), intent(in)    :: x(n)
          real(dp), intent(out)   :: y(n)
          integer                 :: i
          ! NOT DONE .. currently prconditioning matrix is I !!!
          y = 0.0_dp
          do i = 1, n 
              y(i) = x(i)/m(i)
          end do
        end subroutine Msolve 

    end subroutine minres_ez_solve 


    subroutine minres_ez_print(this)
        class(minres_ez_t), intent(in) :: this
        print *, ''
        print *, 'minres_ez'
        print *, '---------------------------------------------------'
        print *, 'itnlim = ', this % itnlim
        print *, 'nout   = ', this % nout  
        print *, 'precon = ', this % precon
        print *, 'checka = ', this % checka
        print *, 'rtol   = ', this % rtol  
        print *, ''
    end subroutine minres_ez_print


    ! minres_info_t methods
    ! -----------------------------------------------------------------------------------

    subroutine minres_info_print(this)
        class(minres_info_t), intent(in) :: this
        print *, ''
        print *, 'minres_info'
        print *, '---------------------------------------------------'
        print *, 'istop  =', this % istop  
        print *, 'itn    =', this % itn    
        print *, 'anorm  =', this % anorm  
        print *, 'acond  =', this % acond  
        print *, 'rnorm  =', this % rnorm  
        print *, 'arnorm =', this % arnorm 
        print *, 'ynorm  =', this % ynorm  
        print *, ''
    end subroutine minres_info_print

    ! utility functions
    ! -----------------------------------------------------------------------------------

    ! Creates a diagonal preconditioner where diagonal elements are from same as A matrix
    subroutine create_precon(irow, icol, a, nnz, sz, m)
        integer,  intent(in)                 :: irow(:)
        integer,  intent(in)                 :: icol(:)
        real(dp), intent(in)                 :: a(:)
        integer,  intent(in)                 :: nnz 
        integer,  intent(in)                 :: sz
        real(dp), intent(inout), allocatable :: m(:)
        integer                              :: i
        allocate(m(sz), source=1.0_dp)
        do i = 1, nnz
            if ((irow(i) == icol(i)) .and. (a(i) > epsilon(1.0_dp))) then
                m(irow(i)) = a(i)
            end if
        end do
    end subroutine create_precon

end module minres_ez_mod
