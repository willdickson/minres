![MINRES](media/logo_smaller.png)
=================================

An implementation of the Minimum Residual Method (MINRES) based on the original
fortran 90 code by Chris Paige, Sou-Cheng Choi, and Michael Saunders. 

### Brief Description

MINRES is designed to solve the system of linear equations

   Ax = b

or the least-squares problem

   min ||Ax - b||\_2,

where A is an n by n symmetric matrix and b is a given vector.  The matrix A
may be indefinite and/or singular.

More generally, MINRES is designed to solve the system

   (A - shift\*I) x = b
or
   min ||(A - shift\*I) x - b||\_2,

where  shift  is a specified scalar value.  Again, the matrix (A - shift\*I)
may be indefinite and/or singular.  The work per iteration is very slightly
less if  shift = 0.

I've added two things: 

1. made it into a package which can be built using the Fortan Package Manager fpm. 
2. created the minres\_ez\_t class which provides a simplified interface to the method.

### Usage

There are two different ways to use this library.  The first is via the
`minres_ez_t` class which provides a simplified interface to the method and the
second is via the `minres_solver` subroutine which gives you full control, but
requires a more work on the part of the user. 

#### minres\_ez\_t

To use the `minres_ez_t` class, you have to provide the matrix `A` in sparse
form, using three arrays: the row indices, column indices, and the nonzero
elements.  Here is an example:


```fortran
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
    minres_ez % itnlim  = 100          ! Upper limit on the nnzber of iterations.
    minres_ez % precon  = .false.      ! Whether or not to invoke preconditioning 
    minres_ez % checka  = .true.       ! Whether or not to check if matrix A is symmetric 
    minres_ez % rtol    = 1.0e-13_dp   ! User-specified residual tolerance 

    ! display settings
    call minres_ez % print()

    ! Solve linear system
    call minres_ez % solve(irow, icol, a, b, x, minres_info)

    ! Print info on result
    call minres_info % print()

end program minres_ez_example_1
```

Two additional optional arguments can be passed to the solve method nnz and shift. 

```fortran
integer,  intent(in), optional   :: nnz      ! number of nonzero items
real(dp), intent(in), optional   :: shift    ! shift value (A - shift*I) x = b
```
If nnz is not specified the minres\_ez\_t class assumes the number of nonzero
elements in A is equal to size(a).  If this is not the case, e.g. the number of
nonzero elements is less then size(a), then nnz can be specified. 

The minres_info_t type contains the following members. 

```fortran
integer  :: istop  ! integer giving the reason for termination
integer  :: itn    ! nnzber of iterations performed
real(dp) :: anorm  ! estimate of the norm of the matrix operator
real(dp) :: acond  ! estimate of the condition of Abar 
real(dp) :: rnorm  ! estimate of norm of residual vector
real(dp) :: arnorm ! recognize singular systems ||Ar||
real(dp) :: ynorm  ! estimate of the norm of xbar
```

#### minres\_solver

The `minres_solver` subroutine  (To do ...)


### Compiling

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the library and test cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `minres` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
minres.git = "https://github.com/willdickson/minres.git"
```

### Documentation

(To do..)

### License

The original fortran 90 implementation of MINRES, everything in the src/minres sub-directory, is distributed under the OSI Common Public License (CPL):
http://www.opensource.org/licenses/cpl1.0.php

My additions, everything else, are distributed under the MIT license https://opensource.org/licenses/MIT. 


### References

1. C. C. Paige and M. A. Saunders (1975). Solution of sparse indefinite systems of linear equations, SIAM J. Numerical Analysis 12, 617-629.
2. S.-C. Choi (2006). Iterative Methods for Singular Linear Equations and Least-Squares Problems, PhD thesis, Stanford University.
3. S.-C. T. Choi, C. C. Paige and M. A. Saunders. MINRES-QLP: A Krylov subspace method for indefinite or singular symmetric systems, SIAM J. Sci. Comput. 33:4, 1810-1836, published electronically Aug 4, 2011.
4. David Chin-Lung Fong and Michael Saunders. CG versus MINRES: An empirical comparison, SQU Journal for Science, 17:1 (2012), 44-62.https://stanford.edu/group/SOL/reports/SOL-2011-2R.pdf


