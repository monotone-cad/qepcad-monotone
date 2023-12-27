/*======================================================================
                    Q <- PPREPVS(P, k)

Adds additional variables to a saclib polynomial.

\Input
    P : polynomial in K[x_1,...,x_r]
    k : non-negative integer, number of variables to add
Output
    Q : P written in K[y_1,...,y_k,x_1,...,x_r]

======================================================================*/
#include "qepcad.h"

// helper function to construct a zero polynomial on top of which the original will be written
Word PrepHelper(Word k, Word P)
{
    // base case: polynomial in one variable
    if (k == 1) return LIST2(0, P);

    // recursive case
    return LIST2(0, PrepHelper(k - 1, P));
}

Word PPREPVS(Word P, Word k)
{
    if (k == 0) return P;

    Word P1 = NIL;
    Word e, Q;
    while (P != NIL) {
        ADV2(P, &e, &Q, &P);

        P1 = COMP2(PrepHelper(k, Q), e, P1);
    }

    return INV(P1);
}

