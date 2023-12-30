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

Word ROOTS(Word Ps, Word I)
{
    Word B, s, T, Rs, Rs1, Rs2, M, J, l, r;
    LWRITE(I); SWRITE("\n");

    // compute the list of "similar polynomials", to satisfy the condition for finest squarefree basis
    IPLSRP(Ps, &s, &T);

    // compute basis B
    B = IPFSFB(1, T);

    // represent as a list (J, M, ...) where J is the isolating interval for the unique root of polynomial M
    Rs = IPLRRI(B);
    LWRITE(Rs); SWRITE("\n");

    // no roots
    if (Rs == NIL) return NIL;

    // finally, throw out roots outside the desired interval
    FIRST2(I, &l, &r);

    printf("total roots %d\n", LENGTH(Rs));
    if (l != NIL) {
        do {
            ADV2(Rs, &J, &M, &Rs);
            RNWRITE(FIRST(J)); RNWRITE(l); printf(" comp = %d \n", RNCOMP(FIRST(J),l));
        } while(RNCOMP(FIRST(J),l) > 0 && Rs != NIL);
    }

    printf("after l %d\n", LENGTH(Rs));
    // Rs1 now contains all roots whose isolating interval is greater than l

    if (r == NIL) return Rs;
    // TODO if rational don't we need to compute it? test pls!
    Rs1 = Rs;
    while (Rs1 != NIL) {
        Rs2 = Rs1;
        ADV2(Rs1, &J, &M, &Rs1);

        if (RNCOMP(SECOND(J), r) < 0) break;
    }
    // Rs1 now points to the las root before r
    printf("before r %d\n", LENGTH(Rs));

    if (Rs2 != NIL) SRED(Rs1, NIL);

    return Rs;
}

