!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File minres_test.f90
!
!    minres_test
!
! Main program for testing minres_solver via subroutine minrestest in 
! minres_test_mod.
!
! Maintained by Michael Saunders <saunders@stanford.edu>
!             and Sou-Cheng Choi <scchoi@stanford.edu>.
!
! 11 Oct 2007: Initially used compiler option -r8.
! 15 Oct 2007: Use real(kind=8) everywhere.
! 16 Oct 2007: Use minresDataModule to define dp = selected_real_kind(15).
! 04 Dec 2007: Debugged singular compatible and incompatible examples.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program minres_test

  use minres_data_mod, only : dp
  use minres_test_mod, only : minrestest, minrestest2

  implicit none

  !---------------------------------------------------------------------
  ! This program calls minrestest(...) to generate a series of test problems
  ! Ax = b or Ax ~= b and solve them with MINRES.
  ! The matrix A is n x n.  It is defined by routines in minres_test_mod.
  !
  ! 11 Oct 2007: First version of minres_test.f90.
  ! 15 Oct 2007: minrestest2 added.
  !---------------------------------------------------------------------

  ! Local variables
  logical                   :: normal, precon, ifLS
  integer                   :: n, nout
  real(dp)                  :: shift, pertM
  character(:), allocatable :: filename


  ! Local constants    
  real(dp), parameter :: zero = 0.0_dp

  filename = 'MINRES_tests.txt'
  open(newunit=nout,file=filename,status='unknown')

  normal = .false.
  precon = .true.
  shift  = 0.25_dp
  pertM  = 0.1_dp
  
  ! Test the unlikely tiny cases that often trip us up.

  n      = 1
  call minrestest( n, normal, zero , zero, nout  )
  call minrestest( n, normal, shift, zero, nout  )

  n      = 2
  call minrestest( n, normal, zero , zero, nout  )
  call minrestest( n, normal, shift, zero, nout  )

  ! Minrestest small positive-definite and indefinite systems
  ! without preconditioners.  MINRES should take n iterations.

  n      = 50
  call minrestest( n, normal, zero , zero, nout  )
  call minrestest( n, normal, shift, zero, nout  )

  ! if (n <= 5) stop    ! WHILE TESTING

  ! Test small positive-definite and indefinite systems with
  ! exact preconditioners.  MINRES should take 1 and 2 iterations
  ! respectively.

  n      = 1
  call minrestest( n, precon, zero , zero, nout  )
  call minrestest( n, precon, shift, zero, nout  )

  n      = 2
  call minrestest( n, precon, zero , zero, nout  )
  call minrestest( n, precon, shift, zero, nout  )

  n      = 50
  call minrestest( n, precon, zero , zero, nout  )
  call minrestest( n, precon, shift, zero, nout  )

  ! pertM makes the preconditioners incorrect in n/10 entries.
  ! MINRES should take about n/10 iterations.

  call minrestest( n, precon, zero , pertM, nout )
  call minrestest( n, precon, shift, pertM, nout )


  ! Singular but compatible test cases.  No precon
  n      = 1
  ifLS   = .false.
  call minrestest2( n, normal, zero , zero, nout, ifLS )

  n      = 2
  ifLS   = .false.
  call minrestest2( n, normal, zero , zero, nout, ifLS )

  n      = 50
  ifLS   = .false.
  call minrestest2( n, normal, zero , zero, nout, ifLS )

  
  ! Singular and incompatible test cases.  No precon
  n      = 1
  ifLS   = .true.
  call minrestest2( n, normal, zero , zero, nout, ifLS )

  n      = 2
  ifLS   = .true.
  call minrestest2( n, normal, zero , zero, nout, ifLS )

  n      = 50
  ifLS   = .true.
  call minrestest2( n, normal, zero , zero, nout, ifLS )


  ! Singular but compatible test cases with precon
  n      = 1
  ifLS   = .false.
  call minrestest2( n, precon, zero , zero, nout, ifLS )

  n      = 2
  ifLS   = .false.
  call minrestest2( n, precon, zero , zero, nout, ifLS )

  n      = 50
  ifLS   = .false.
  call minrestest2( n, precon, zero , pertM, nout, ifLS )

  print *, 'MINRES tests complete. Test data written to ', filename

end program minres_test
