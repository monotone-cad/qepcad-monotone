/*======================================================================
 * (J_1, M_1, ..., J_l, M_l) <- ROOTS(Ps)
 *
 * Real root isolation for a list of univariate polynomials with rational coefficients
 *
 * Input
 *   Ps : a list of univariate polynomials with rational coefficients
 * Output
 *   M_i,J_i : represents a real root of a polynomial in P. M_i is a univariate polynomial and J_i is an isolating
 *             interval for the root
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

