# QEPCAD Monotone

This is a modified version of QEPCAD which is able to compute a CAD, compatible with a family of semialgebraic sets of
dimension at most two and satisfying the frontier condition, based on the work of Basu, Gabrielov and Vorobjov (2015)
"Triangulation of Monotone Families"

```
============================================================
                     Q E P C A D   B

                     A program for

    Quantifier Elimination and Formula Simplification

               in Elementary Real Algebra via

             Cylindrical Algebraic Decomposition
============================================================

Please refer to: https://www.usna.edu/CS/qepcadweb/B/QEPCAD.html

Making QEPCAD
-------------
1. Install saclib, and make sure that the "saclib" environment
   variable gives the full path to your saclib root directory.
2. Set the environment "qe" to the full path to your qepcad
   installation root directory - i.e. to this directory!
3. Type "make" in this directory.  You may adjust flags, etc.
   in the file "Makefile" in this directory.
At this point, "$qe/bin/qepcad" should run QEPCAD.  You might
want to add $qe/bin to your path.

Using Singular as a Computer Algebra Server
-------------------------------------------
QEPCAD B allows you to run Singular to farm out some computations
that might not be as fast in Saclib, like multivariate polynomial
factorization, and to provide Groebner Basis computation, which
is not available in Saclib.  When Singular is installed on your
system, add the line

  SINGULAR <full-path-for-singular-binary-directory>

to the default.qepcadrc file.  For example:

  SINGULAR /usr/lib/Singular/3-0-0/ix86-Linux

might be right for you.

Making CAD2D
-------------
The program "cad2d" is a variant of QEPCAD with some
optimizations for constructing 2D CADs.  Cd to the adj2d
directory and type "gmake", AFTER INSTALLING QEPCAD!  After
this, "$qe/bin/cad2d" should launch the CAD2D program.  Note:
you must use "gmake".

Making ADJ2D_plot
----------------
The ADJ2D_plot program is the rendering engine called up by
QEPCADB's plot-2d-cad command.  Just change to the "plot2d"
directory and type "make" to create it.  You may need to
adjust some entries in the Makefile to match OpenGL and GLUT
locations on your system.  Note that rendering remotely (i.e.
via ssh tunelling) may not work, while local rendering does.
For me, setting the environment variable LIBGL_ALWAYS_INDIRECT
if remote sessions fixed the problem.  E.g.

  export LIBGL_ALWAYS_INDIRECT=1

If you're having the problem, give this fix a try!

Directories
-----------
bin - links to the qepcad and cad2d programs (as well as some
      other files for QEPCAD)

source - The "original" QEPCAD.  In fact, much has been modified
      and added since the original, but we distinguish this
      from the "extensions".

extensions - Extension packages to QEPCAD:
      sfext: solution formula extensions
      adj2d: 2D adjacency extensions
      rend : 2D CAD plotting extensions
      lift2D:Extensions for fast lifting via floating-point for
             2D CADs.
      newadj:Experimental extensions for modified adjacency
             computations.

cad2d - This contains the program "cad2d" is a variant of
      QEPCAD with some optimizations for constructing 2D CADs.

plot2d - The ADJ2D_plot program.
```
