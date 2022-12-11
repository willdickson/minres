# minres

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
minre\_ez\_t class which provides a simplified interface to the method and the
second is via the minres_solver function which give you more control, but
requires a bit more work on the part of the user. 

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


