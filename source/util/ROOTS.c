/*======================================================================
 * (J, M, ...) <- ROOTS(Ps)
 *
 & Real root isolation for a univariate polynomal with integer coefficients
 *
 *====================================================================*/
#include "qepcad.h"

Word ROOTS(Word Ps)
{
    Word s, T, B, E;

    // construct list of similar rational polynomials
    IPLSRP(Ps, &s, &T);

    // compute basis B
    IPFSBM(1, T, &B, &E);

    // finally, represent as a list (J, M, ...) where J is the isolating interval for the unique root of polynomial M
    return IPLRRI(B);
}

