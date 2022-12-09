!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File minres_test_mod.f90
!
! This file illustrates how minres_solver can call Aprod with a
! short fixed parameter list, even if it needs arbitrary other data.
!
! 11 Oct 2007: First version of minres_test_mod.f90.
! 15 Oct 2007: Use real(8) everywhere.
!              minrestest2 added.
! 16 Oct 2007: Use minresDataModule to define dp = selected_real_kind(15).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module minres_test_mod

  use  minres_data_mod, only : dp
  use  minres_mod,      only : minres_solver

  implicit none
  public   :: minrestest, minrestest2
  private  :: Aprod, Msolve, Aprod2

  ! DYNAMIC WORKSPACE DEFINED HERE.
  ! It is allocated in minrestest and used by Aprod.

  real(dp), allocatable :: d(:)     ! Defines diagonal matrix D.
  real(dp)              :: Ashift   ! Global version of shift.  A = D - Ashift*I.
  real(dp)              :: Mpert    ! Global version of pertM.  Perturbation to D
                                    ! in Msolve to avoid having an exact preconditioner.

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod (n,x,y)

    integer,  intent(in)    :: n
    real(dp), intent(in)    :: x(n)
    real(dp), intent(out)   :: y(n)

    !-------------------------------------------------------------------
    ! Aprod  computes y = A*x for some matrix A.
    ! This is a simple example for testing minres_solver.
    !-------------------------------------------------------------------

    integer  :: i

    do i = 1, n
       y(i) = d(i)*x(i)
    end do

  end subroutine Aprod
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod2 (n,x,y)

    integer,  intent(in)    :: n
    real(dp), intent(in)    :: x(n)
    real(dp), intent(out)   :: y(n)

    !-------------------------------------------------------------------
    ! Aprod  computes y = A*x for some matrix A.
    ! This is a simple example for testing minres_solver.
    !-------------------------------------------------------------------

    integer  :: i

    do i = 1, n
       y(i) = d(i)*x(i)
    end do
    if (n >= 3) then
       y(n-1:n) = 0.0_dp
    end if

  end subroutine Aprod2

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Msolve(n,x,y)
      
    integer,  intent(in)    :: n
    real(dp), intent(in)    :: x(n)
    real(dp), intent(out)   :: y(n)

    !-------------------------------------------------------------------
    ! Msolve solves M*y = x for some symmetric positive-definite matrix M.
    ! This is a simple example for testing minres_solver.
    ! Ashift will be the same as shift in minres_solver.
    !
    ! If Mpert = 0, the preconditioner will be exact, so
    ! minres_solver should require either one or two iterations,
    ! depending on whether (A - shift*I) is positive definite or not.
    !
    ! If Mpert is nonzero, somewhat more iterations will be required.
    !
    ! 04 Dec 2007: Safeguard singular M.
    !-------------------------------------------------------------------
      
    intrinsic     :: abs, mod      
    integer       :: i
    real(dp)      :: di

    do i = 1, n
       di   = abs( d(i) - Ashift )
       if (mod(i,10) == 0  ) di = di + Mpert
       if (abs(di) > 0.0_dp) then
          y(i) = x(i) / di
       else
          y(i) = x(i)
       end if
    end do

  end subroutine Msolve

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine minrestest( n, precon, shift, pertM, nout )

    integer,  intent(in)    :: n, nout
    logical,  intent(in)    :: precon  
    real(dp), intent(in)    :: shift, pertM

    !-------------------------------------------------------------------
    ! minrestest solves sets up and solves a system (A - shift*I)x = b,
    ! using Aprod to define A and Msolve to define a preconditioner.
    !-------------------------------------------------------------------

    intrinsic :: abs, dot_product, real, sqrt

    ! Local arrays and variables
    real(dp)  :: b(n), r1(n), w(n), x(n), xtrue(n), y(n)
    logical   :: checkA 
    integer   :: j, itnlim, istop, itn, nprint
    real(dp)  :: Anorm, Acond, Arnorm, rnorm, rtol, r1norm, ynorm
    real(dp)  :: enorm, etol, wnorm, xnorm

      
    write(nout, 1000) shift, pertM

    allocate( d(n) )         ! Array used in Aprod and Msolve to define A and M.
    Ashift = shift           ! Global version of shift
    Mpert  = pertM           ! Global version of pertM

    do j = 1, n              ! Set d(*) to define A.
       d(j) = real(j)        ! Does this give a full-precision integer?
       d(j) = d(j)/n         ! We don't want exact integers.
       xtrue(j) = real(n+1-j)! Set the true solution and the rhs
    end do                   ! so that (A - shift*I)*xtrue = b.

    call Aprod (n,xtrue,b)   ! b = A*xtrue
    b      = b - shift*xtrue ! Now b = (A - shift*I)*xtrue

    checkA = .true.          ! Set other parameters and solve.
    itnlim = n*2
    rtol   = 1.0e-12_dp

    call minres_solver( n, Aprod, Msolve, b, shift, checkA, precon, &
                 x, itnlim, nout, rtol,                      &
                 istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm )

    call Aprod (n,x,y)       ! y = A*x
    r1     = b - y + shift*x ! Final residual r1 = b - (A - shift*I)*x.    
    r1norm = sqrt(dot_product(r1,r1))
    write(nout, 2000) r1norm

    nprint = min(n,20)
    write(nout,2500) (j, x(j), j=1,nprint)  ! Print some of the solution

    w      = x - xtrue                      ! Print a clue about whether
    wnorm  = sqrt(dot_product(w,w))         ! the solution looks OK.
    xnorm  = sqrt(dot_product(xtrue,xtrue))
    enorm  = wnorm/xnorm
    etol   = 1.0e-5_dp
    if (enorm <= etol) then
       write(nout, 3000) enorm
    else
       write(nout, 3100) enorm
    end if

    deallocate(d)                           ! Free work array
    return

 1000 format(//' -----------------------------------------------------' &
             / ' Test of  minres_solver.'                                      &
             / ' -----------------------------------------------------' &
             / ' shift =', f12.4, 6x, 'pertM =', f12.4)
 2000 format(/ ' Final residual =', 1p, e8.1)
 2500 format(/ ' Solution  x' / 1p, 4(i6, e14.6))
 3000 format(/ ' minres_solver  appears to be successful.  Relative error in x =', 1p, e8.1)
 3100 format(/ ' minres_solver  appears to have failed.    Relative error in x =', 1p, e8.1)

  end subroutine minrestest

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine minrestest2( n, precon, shift, pertM, nout, ifLS)

    integer,  intent(in)     :: n, nout
    logical,  intent(in)     :: precon  
    real(dp), intent(in)     :: shift, pertM

    !-------------------------------------------------------------------
    ! minrestest solves sets up and solves a system (A - shift*I)x = b,
    ! using Aprod2 to define A and Msolve to define a preconditioner.
    !-------------------------------------------------------------------

    intrinsic  :: abs, dot_product, real, sqrt

    ! Local arrays and variables
    real(dp)   :: b(n), r1(n), w(n), x(n), xtrue(n), y(n)
    logical    :: checkA, ifLS
    integer    :: j, itnlim, istop, itn, nprint
    real(dp)   :: Anorm, Acond, Arnorm, rnorm, rtol, r1norm, ynorm
    real(dp)   :: enorm, etol, wnorm, xnorm

      
    write(nout, 1000) shift, pertM

    allocate( d(n) )         ! Array used in Aprod2 and Msolve to define A and M.
    Ashift = shift           ! Global version of shift
    Mpert  = pertM           ! Global version of pertM

    do j = 1, n              ! Set d(*) to define A.
       d(j) = real(j)        ! Does this give a full-precision integer?
       d(j) = d(j)/n         ! We don't want exact integers.
       xtrue(j) = real(n+1-j)! Set the true solution and the rhs
    end do                   ! so that (A - shift*I)*xtrue = b.

    if (n >= 3) then
       xtrue(n-1:n) = 0.0_dp
       d(n-1:n)     = 0.0_dp
       x(n-1:n)     = 0.0_dp
    end if
    
    call Aprod2 (n,xtrue,b)   ! b = A*xtrue
    b      = b - shift*xtrue  ! Now b = (A - shift*I)*xtrue

    if (n >= 3 .and. ifLS) then
       b(n-1:n) = 1.0_dp      ! b not in the range of A
    end if
    
    checkA = .true.           ! Set other parameters and solve.
    itnlim = n*2
    rtol   = 1.0e-12_dp

  ! write(nout,2501) (j, d(j), j=1,n)     ! print A 
  ! write(nout,2502) (j, xtrue(j), j=1,n) ! print xtrue
  ! write(nout,2503) (j, b(j), j=1,n)     ! print b
      
    call minres_solver( n, Aprod2, Msolve, b, shift, checkA, precon, &
                 x, itnlim, nout, rtol,                       &
                 istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm )

    call Aprod2 (n,x,y)       ! y = A*x
    r1     = b - y + shift*x  ! Final residual r1 = b - (A - shift*I)*x.    
    r1norm = sqrt(dot_product(r1,r1))
    write(nout, 2000) r1norm

    call Aprod2 (n,r1,w)      ! w = A*r1  (should be small)
    Arnorm = sqrt(dot_product(w,w))


    nprint = min(n,50)
    
    write(nout,2500) (j, x(j), j=1,nprint)  ! Print some of the solution

    w      = x - xtrue                      ! Print a clue about whether
    wnorm  = sqrt(dot_product(w,w))         ! the solution looks OK.
    xnorm  = sqrt(dot_product(xtrue,xtrue))
    enorm  = wnorm/xnorm
    etol   = 1.0e-5_dp
    if (enorm <= etol) then
       write(nout, 3000) enorm
    else if (Arnorm > etol) then
       write(nout, 3100) enorm
    else
       write(nout, 3200)
    end if

    deallocate(d)                           ! Free work array
    return

 1000 format(//' -----------------------------------------------------' &
             / ' Test of  minres_solver.'                                      &
             / ' -----------------------------------------------------' &
             / ' shift =', f12.4, 6x, 'pertM =', f12.4)
 2000 format(/ ' Final residual =', 1p, e8.1)
!2501 format(/ ' A          ' / 1p, 4(i6, e14.6))
!2502 format(/ ' xtrue      ' / 1p, 4(i6, e14.6))
!2503 format(/ ' b          ' / 1p, 4(i6, e14.6))
 2500 format(/ ' Solution  x' / 1p, 4(i6, e14.6))
 3000 format(/ ' minres_solver  appears to be successful.  Relative error in x =', 1p, e8.1)
 3100 format(/ ' minres_solver  appears to have failed.    Relative error in x =', 1p, e8.1)
 3200 format(/ ' minres_solver  appears to have found a least-squares solution')

  end subroutine minrestest2

end module minres_test_mod
