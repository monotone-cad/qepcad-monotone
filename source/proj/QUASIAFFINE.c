/*======================================================================
  QUASIAFFINE(A, r, *A_);

considers the critical points of projections of the smooth two-dimensional locus of input set, onto one and two
dimensional coordinate subspaces.

\Input
  \parm{A} set of projection factors (A_1,...,A_r), each A_i is a list of polynomials in Z[x_1,...,x_i]
  \parm{r} positive integer, number of variables

Output
  \parm{*A}: set of projection factors, modified to include first partial derivatives of each element with respect to
  each vaciable.

SideEffect
  \parm{A} is also modified.

======================================================================*/
#include "qepcad.h"

void QepcadCls::QUASIAFFINE(Word A, Word r, Word *A_)
{
    // strat methos:
    // - construct all sign sets of polynomials. 3^n of them. some will be empty
    // - for each sign set, call stratify.
    // - get back a list of strata.
    // - for each 2-dimensional stratum, add jacobi det for each of the sets of variables, excluding two.

    *A_ = A;
}

