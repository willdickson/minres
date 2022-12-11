# minres

An implementation of the MINRES algorithm which based on the original fortran 90
code by Chris Paige, Sou-Cheng Choi, and Michael Saunders. 

I've added two things: 

1. I've made it into a package which can be built using the fortan package manager fpm. 
2. I've created an "easy mode" wrapper type so it can be used in a simplified manner. 

## The MINRES Algorithm 

MINRES is designed to solve the system of linear equations

   Ax = b

or the least-squares problem

   min ||Ax - b||,

where A is an n by n symmetric matrix and b is a given vector.  The matrix A
may be indefinite and/or singular.

    
(In development).

