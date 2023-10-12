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

        // refactor evaluation code from MONOTONE, eval P at sample x and append.
        A1 = SUBSTITUTE(r, A, S, false); // w/ integer coefficients
        // LWRITE(A); SWRITE("\n");

        if (PDEG(A1) == 0) continue;

        Ps1 = COMP(A1, Ps1);
        // SWRITE("  -> "); LWRITE(SUBSTITUTE(r, A, S, false)); SWRITE("\n");
    }

    return Ps1;
}

void LazardHelper(Word r, Word Ps, Word S, Word As, Word i, Word k, Word R, Word L, Word* Hs_)
{
    // base case, blow-up polynomial f
    if (i == r) {
        Word f = FIRST(Ps);
        printf("base case: f = "); LWRITE(f); SWRITE("\n");

        *Hs_ = CONC(*Hs_, SectionCells(
            r,
            S,
            GVCAP->GROEBNER(
                COMP(f, L),
                COMP(r, R),
                r + 1
            )
        ));

        return;
    }

    // saturation polynomial h, 1 - z * h
    if (i == k) {
        Word g;
        ADV(Ps, &g, &Ps);

        // construct 1 - z * h
        Word h = IPSUM(
            r + 1,
            PPREPVS(PMON(1,0), r), // 1
            LIST2(1, PADDVS(IPNEG(k, g), r - k)) // z * - h
        );

        printf("%d, 1 - zh case = ", i); LWRITE(h);
        LazardHelper(
            r,
            Ps,
            S,
            RED(As),
            i + 1,
            k,
            COMP(r + 1, R),
            COMP(h, L),
            Hs_
        );

        return;
    }

    // rest of elements are projection factors.
    Word A1s, A11;
    ADV(As, &A1s, &As);

    while (A1s != NIL) {
        ADV(A1s, &A11, &A1s);

        printf("%d, proj factor case = ", i); LWRITE(LELTI(A11, PO_POLY)); SWRITE("\n");
        LazardHelper(
            r,
            RED(Ps),
            S,
            As,
            i + 1,
            k,
            COMP(i, R),
            COMP(LELTI(A11, PO_POLY), L),
            Hs_
        );
    }
}

Word LazardLifting(Word k, Word Ps, Word S, Word As)
{
    Word Hs = NIL;

    Word i = 0;
    while (i < k) {
        ++i;

        LazardHelper(k + 1, Ps, S, As, 1, i, NIL, NIL, &Hs);
    }

    return Hs;
}

