/*======================================================================
 * B <- IPFRP(r,A)
 *
 * Integral polynomial from rational polynomial (modified to handle ratinoal coefficients).
 *
 * Inputs
 *     r : a BETA-digit, r >= 0, the number of variables.
 *     A : in Q[X1,...,Xr], The base coefficients of A are not assumed to be integers.
 *
 * Outputs
 *     B : in Z[X1,...,Xr], B is A converted to integral representation.
======================================================================*/
#include "qepcad.h"

// rational univariate polynomial convert to integer coefficients
inline Word RUPConvert(Word A)
{
    Word lcm = 1, Es = NIL, As = NIL, Bs = NIL;

    while (A != NIL) {
        Word a, e, a1, a2;
        ADV2(A, &e, &a, &A);

        // store the exponents and coefficients for reconstruction later
        Es = COMP(e, Es);
        if (a != 0) {
            FIRST2(a, &a1, &a2);
        } else {
            a1 = 0;
            a2 = 1;
        }

        As = COMP(a1, As);
        Bs = COMP(a2, Bs);

        // update lcm for integer conversion
        if (lcm == 1) { // set lcm
            lcm = a2;
        } else if (a2 != 1) { // compute lcm (if needed)
            lcm = ILCM(lcm, a2);
        }
    }

    // reconstruct (backwards for return), multiplying the coefficients by the lcm of their denominators
    Word B = NIL;
    while (As != NIL && Es != NIL) {
        Word a, b, e;
        ADV(As, &a, &As);
        ADV(Bs, &b, &Bs);
        ADV(Es, &e, &Es);

        Word c = IQ(IPROD(lcm, a), b);
        B = COMP2(e, c, B);
    }

    return B;
}

Word IPFRPmod(Word r, Word A)
{
    // zero polynomial
    if (A == 0) {
        return 0;
    }

    // base case: r = 1 (convert coefficients)
    if (r == 1) {
        return RUPConvert(A);
    }

    // recursive case: r > 0 (reconstruct polynomial)
    Word As = A;
    Word rp = r - 1;
    Word B = NIL;

    while (As != NIL) {
        Word e, a, b;
        ADV2(As,&e,&a,&As);

        b = IPFRPmod(rp, a);
        B = COMP2(b,e,B);
    }

    return INV(B);
}
