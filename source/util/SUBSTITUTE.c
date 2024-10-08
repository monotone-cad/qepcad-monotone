/*======================================================================
 * Q <- SUBSTITUTE(r, P, S, rc)
 *
 & Substitute sample point S into polynomial P in Z[x_1,...,x_r]
 *
 * Input:
 *     r  : positive integer
 *     P  : polynomial in Z[x_1,...,x_r]
 *     S  : sample point
 *     rc : true to return the polynomial with rational coefficients. false returns integer coefficients
 * Output:
 *     Q  : polynomial in Q[x_{k+1},...,x_r] such that Q = P(s_1,...,s_k,x_{k+1},...,x_r)
 *
 *====================================================================*/
#include "qepcad.h"

// if polynomial deos not depend on variables x_1,..,x_k, then return a polynomial in x_(k+1),...,x_r
Word PNotDependVs(Word r, Word P, Word k)
{
    // base case r == 0, P is an integer
    if (r == 0) {
        return P;
    }

    Word P1, e, Q, Q1;

    // k < r: check for dependent variable
    if (r < k) {
        FIRST2(P, &e, &Q);

        if (e > 0) { // depends on x_i
            return NIL;
        } else {
            return PNotDependVs(r - 1, Q, k);
        }
    }

    // k >= r, reconstruct polynomial
    P1 = NIL;
    while (P != NIL) {
        ADV2(P, &e, &Q, &P);

        Q1 = PNotDependVs(r - 1, Q, k);

        // error case
        if (Q1 == NIL) return NIL;

        // otherwise construct the polynomial
        P1 = COMP2(Q1, e, P1);
    }

    return INV(P1);
}


Word SUBSTITUTE(Word r, Word P, Word S, bool rc)
{
    Word Q, J, c;
    if (LENGTH(S) == 5) {
        // convert sample
        Word SQ, SJ, SM, SI, Sb;
        FIRST5(S, &SQ, &SJ, &SM, &SI, &Sb);
        ConvertToPrimitive(SQ, SJ, SM, SI, Sb, &Q, &J, &c);
    } else {
        FIRST3(S, &Q, &J, &c);
    }

    Word i0 = LENGTH(c);

    // sample is a subset of R^0, nothing to do
    if (i0 == 0) {
        if (rc) {
            return RPFIP(r, P);
        }

        return P;
    }

    // number of coordinates in sample point
    Word k = r - i0;

    // if P does not depend on all x1,...,kk, then no substitution needed
    Word P1 = PNotDependVs(r, P, i0 + 1);
    if (P1 != NIL) {
        if (rc) {
            return RPFIP(k, P1);
        }

        return P1;
    }

    // otherwise perform substitution
    if (PDEG(Q) > 1) { // algebraic sample
        P1 = IPAFME(r, Q, P, c);
        if (k == 0 || P1 == 0) return P1;

        P1 = AFPNORM(k, Q, P1);

        if (rc) {
            return RPFIP(k, P1);
        }

        return P1;
    }

    // rational sample
    P1 = IPRNME(r, P, RCFAFC(c));
    if (rc) {
        return P1;
    }

    return IPFRPmod(k, P1);
}

