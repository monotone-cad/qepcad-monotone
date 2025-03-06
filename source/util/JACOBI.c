/*======================================================================
 * J <- JACOBI(r, f, i, Hs, Is)
 * Partial differential operator, similar to a Jacobi Matrix. returns its determinant
 *
 * Input:
 *     r  : positive integer
 *     f  : polynomial in Z[x_1,...,x_r]
 *     i  : positive integer, i <= r
 *     Ps : list of polynomial (f1,...,fk) in Z[x_1,...,x_r]
 *     Is : list of positive integers (i1,...,ik), ij <= r, 1 <= j <= k
 * Output:
 *     J  : determinant of jacobi matrix of partial derivatives
 *              [ d h1 / d x_i1, ..., d h1 / d x_ik, d h1 / x_i ]
 *              [ d h2 / d x_i1, ..., d h2 / d x_ik, d h2 / x_i ]
 *          det [                ...                            ]
 *              [ d hk / d x_i1, ..., d h2 / d x_ik, d hk / x_i ]
 *              [ d f  / d x_i1, ..., d f  / d x_ik, d f  / x_i ]
 *
 *
 *====================================================================*/
#include "qepcad.h"

// construct one row
Word JacobiRow(Word r, Word P, Word Is, Word i)
{
    bool zeroRow = true;
    // rightmost d f / d xi
    Word Mi = NIL, D = 0;

    if (i > 0) {
        D = IPDER(r, P, i);
        if (D != 0) zeroRow = false;

        Mi = COMP(D, Mi);
    }

    // remaining elements
    Word j;
    while (Is != NIL) {
        ADV(Is, &j, &Is);

        Word D = IPDER(r, P, j);
        if (D != 0) zeroRow = false;

        Mi = COMP(D, Mi);
    }

    return zeroRow ? 0 : Mi;
}

Word JACOBI(Word r, Word f, Word i, Word Hs, Word Is)
{
    // base case: ordinary derivative
    if (Hs == NIL) return IPDER(r, f, i);

    // build the jacobian
    Word M = NIL;

    // first k rows in the order INV(h1,...,hk),f
    // order doesnt really matter
    Word h, Ri;
    while (Hs != NIL) {
        ADV(Hs, &h, &Hs);

        Ri = JacobiRow(r, h, Is, i);
        if (Ri == 0) return 0; // zero row => zero determinant
        M = COMP(Ri, M);
    }

    if (i > 0) {
        Ri = JacobiRow(r, f, Is, i);

        if (Ri == 0) return 0; // zero row => zero determinant
        M = COMP(Ri, M);
    }

    return MAIPDE(r, M);
}

