/*======================================================================
 * (J_1, M_1, ..., J_l, M_l) <- ROOTS(Ps, I)
 *
 * Real root isolation for a list of univariate polynomials with rational coefficients
 *
 * Input
 *   Ps : a list of univariate polynomials with rational coefficients
 *   I  : interval in which roots should be contained.
 * Output
 *   M_i,J_i : represents a real root of a polynomial in P. M_i is a univariate polynomial and J_i is an isolating
 *             interval for the root
 *
 *====================================================================*/
#include "qepcad.h"

Word Estimate(Word J, Word M)
{
    if (PDEG(M) == 1) {
        return IUPRLP(M);
    }

    // algebraic
    return RNQ(RNSUM(FIRST(J), SECOND(J)), RNINT(2));
}

Word ROOTS(Word Ps, Word I)
{
    Word B, s, T, Rs, Rs1, M, J, l, r;
    // compute the list of "similar polynomials", to satisfy the condition for finest squarefree basis
    IPLSRP(Ps, &s, &T);

    // compute basis B
    B = IPFSFB(1, T);

    // represent as a list (J, M, ...) where J is the isolating interval for the unique root of polynomial M
    Rs = IPLRRI(B);

    // no roots
    if (Rs == NIL) return NIL;

    // finally, throw out roots outside the desired interval
    FIRST2(I, &l, &r);

    Rs1 = Rs;
    if (l != NIL) {
        Word c;
        do {
            ADV2(Rs1, &J, &M, &Rs1);
            c = Estimate(J, M);

            // root is IN if c > l
            if (RNCOMP(c,l) > 0) break;
            Rs = Rs1;
        } while(Rs != NIL);
    }

    // Rs1 now contains all roots whose isolating interval is greater than l

    if (r == NIL || Rs == NIL) return Rs;
    Rs1 = NIL;
    while (Rs != NIL) {
        ADV2(Rs, &J, &M, &Rs);
        Word c = Estimate(J, M);
        // root is IN if c < r
        if (RNCOMP(c,r) >= 0) break;

        Rs1 = COMP2(M, J, Rs1);
    }

    return INV(Rs1);
}

