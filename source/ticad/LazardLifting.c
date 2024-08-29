/*======================================================================
 * L = LazardLifting(Ps, S, As)
 *
 * Perform Lazard Lifting. Ps := [g_1,...,g_k,f],
 * where p = { g_1 = 0, ..., g_k = 0 } defines a 0-dimensional cell
 *       s = p is a k-dimensional sample point
 * and f is a function containing a blow-up above p.
 *====================================================================*/
#include "qepcad.h"
#include "db/CAPolicy.h"
using namespace std;

// from the list of polynomials (Grobner basis) in r+1 variables, find polynomials defining refinement points.
// TODO discard variable r+1 and any polynomials depending on it.
// keep if it has nonzero degree in variable r
// eval at p = f1,f2 and find roots.
Word SectionCells(Word r, Word S, Word Ps)
{
    Word Ps1 = NIL;
    while (Ps != NIL) {
        Word P, e, A, A1;
        ADV(Ps, &P, &Ps);
        FIRST2(P, &e, &A);

        // if P depends on variable r, discard it
        if (e > 0) continue;

        A1 = SUBSTITUTE(r, A, S, false); // w/ integer coefficients

        if (PDEG(A1) == 0) continue;

        Ps1 = COMP(A1, Ps1);
    }

    return Ps1;
}

// saturation of Ps w.r.t Q
Word LazardHelper(Word r, Word S, Word Ps, Word Rs, Word Q, Word k)
{
    // construct 1 - z * Q
    Word Q1 = IPSUM(
        r + 1,
        PPREPVS(PMON(1,0), r), // 1
        LIST2(1, PADDVS(IPNEG(k, Q), r - k)) // z * - Q
    );

    return SectionCells(
        r,
        S,
        GVCAP->GROEBNER(
            COMP(Q1, Ps),
            COMP(r+1, Rs),
            r + 1
        )
    );
}

Word LazardLifting(Word k, Word S, Word As, Word IPs, Word i, Word j)
{
    // locate polynomials
    Word Rs = NIL, Ps = NIL, Q1 = NIL, r1 = NIL, Q2 = NIL, r2 = NIL;
    Word level = 0, index, P, A1s;
    while (IPs != NIL) {
        ++level;
        ADV(IPs, &index, &IPs);
        ADV(As, &A1s, &As);

        P = LELTI(LELTI(A1s, index), PO_POLY);

        if (level == i) {
            Q1 = P;
            r1 = level;
        } else if (level == j) {
            Q2 = P;
            r2 = level;
        } else {
            Ps = COMP(P, Ps);
            Rs = COMP(level, Rs);
        }
    }

    return CONC(
        LazardHelper(k, S, COMP(Q2, Ps), COMP(r2, Rs), Q1, r1),
        LazardHelper(k, S, COMP(Q1, Ps), COMP(r1, Rs), Q2, r2)
    );
}

