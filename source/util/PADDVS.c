/*======================================================================
                    Q <- PADDVS(P, k)

Adds additional variables to a saclib polynomial.

\Input
    P : polynomial in K[x_1,...,x_r]
    k : non-negative integer, number of variables to add
Output
    Q : P written in K[x_1,...,x_r,x_r+1,...,x_r+k]

======================================================================*/
#include "qepcad.h"

Word PADDVS(Word P, Word k)
{
    Word i = 0;

    while (i < k) {
        P = LIST2(0, P);
        ++i;
    }

    return P;
}

