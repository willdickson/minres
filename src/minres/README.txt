README for f90 minres

The software for MINRES (f90 version) is provided by SOL, Stanford University
under the terms of the OSI Common Public License (CPL):
http://www.opensource.org/licenses/cpl1.0.php

Contributors: Chris Paige
              Michael Saunders
              Sou-Cheng Choi

11 Oct 2007: First set of files available for download from SOL.
04 Dec 2007: Debugged singular cases.  MINRES returns *a* least-squares
             solution, but it may not be the minimum-length solution.

Maintained by Michael Saunders, SOL, Stanford University
              saunders@stanford.edu  650-723-1875
and           Sou-Cheng Choi <scchoi@stanford.edu>
-----------------------------------------------------------------------------

The f90 version of MINRES involves the following files:

   minresModule.f90
   minresTestModule.f90
   minresTestProgram.f90
   MINRES.txt       (example output file from an Intel Xeon system
                     compiled with g95 on Linux Redhat 9)
   Makefile

To compile the code and run the test program on Linux or Unix,
proceed as follows:

   Check that Makefile selects an f90 or f95 compiler
   that you actually have.  Edit if necessary.

   Make                 (creates executable TestProgram)
   ./TestProgram
   grep appears MINRES.txt

  "MINRES  appears to be successful" should occur 22 out of 23 times.
The 20th case should say
  "MINRES  appears to have found a least-squares solution"
after minresttest2 is called to solve a singular incompatible system.
For that case, the quantity ||Ar|| should be small, where r = b - Ax.
