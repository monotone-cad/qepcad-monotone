#include "qepcad.h"

// TODO this is not used anywhere
// it is old code for permutation of coordinates and project to obtain monotone cells.

// integer range, from a to b inclusive
Word IRANGE(Word a, Word b);
// swap item at index i with item at index j in the list L (L is modified)
Word LSWAP(Word i, Word j, Word L);
// PPERMV for a list of r-variate polynomials
Word LPERMV(Word r, Word L, Word R);

Word WRAP(Word r, Word P)
{
    return LIST2(0, P);
}

// list of polynomials in r, wrap in brackets to get to r+1
Word LWRAP(Word r, Word L)
{
    Word A1, LL;
    LL = NIL;
    while (L != NIL) {
        ADV(L, &A1, &L);
        Word P = LELTI(A1, PO_POLY);
        LL = COMP(MPOLY(WRAP(r, P), NIL, NIL, PO_OTHER, PO_KEEP), LL);
    }

    return INV(LL);
}

// given list L, return only polynomials with positive degree in main (r-th) variable
Word LPPDEG(Word L)
{
    Word LL, P;

    LL = NIL;
    while (L != NIL) {
        ADV(L, &P, &L);

        if (PDEG(LELTI(P, PO_POLY)) == 0) continue;

        LL = COMP(P, LL);
    }

    return INV(LL);
}

// permute and project
//   A: proj factor structure
//   r: dimension
//   k: dimension to project to
// A is modified to add projections onto k-dimensional subspaces (and returned)
Word QepcadCls::PermuteAndProject(Word A, Word r, Word k)
{
    Word L = NIL;
    // TODO r > 3
    Word A1 = CONC(LCOPY(LELTI(A, r)), LWRAP(r-1, LELTI(A, r-1)));
    Word Seq = IRANGE(1, r);
    Word Seq1 = IRANGE(1, r);

    for (int i = 1; i <= r-1; i++) { /* first coordinate */
        Seq = LSWAP(1, i, Seq);
        for (int j = MAX(3, i+1); j <= r; j++) { /* second coordinate */
            Seq = LSWAP(2, j, Seq);
            // undo
            Seq1 = LSWAP(2, j, Seq1);
            Seq1 = LSWAP(1, i, Seq1);

            LWRITE(Seq);SWRITE("\n");
            Word Pp = LPERMV(r, A1, Seq);
            Pp = LPPDEG(Pp);
            Word Pp1 = PROJ(r, Pp);
            Word Pp2 = LWRAP(r-1, Pp1);
            LWRITE(Seq1);SWRITE("\n");
            Word Pp3 = LPERMV(r, Pp2, Seq1);
            printf("before:%d, not constant: %d, aftor proj: %d\n", LENGTH(A1), LENGTH(Pp), LENGTH(Pp3));

            L = CONC(L, Pp3);

            Seq1 = LSWAP(1, i, Seq1);
            Seq1 = LSWAP(2, j, Seq1);

            Seq = LSWAP(2, j, Seq);
        } /* second coordinate */
        Seq = LSWAP(1, i, Seq);
    } /* first coordinate */

    return L;
}

Word IRANGE(Word a, Word b)
{
    if (a > b) return NIL;

    Word L = NIL;

    int i = b;
    while (i >= a) { /* [do it backwards to avoid INV(L)] */
       L = COMP(i, L);

       i--;
    }

    return L;
}

Word LSWAP(Word i, Word j, Word L)
{
    Word tmp = LELTI(L, i);
    SLELTI(L, i, LELTI(L, j));
    SLELTI(L, j, tmp);

    return L;
}

Word LPERMV(Word r, Word L, Word R)
{
    Word LL, P;

    LL = NIL;
    while (L != NIL) {
        ADV(L, &P, &L);

        IPWRITE(r, LELTI(P, PO_POLY), LIST3(LFS("x1"), LFS("x2"), LFS("x3")));SWRITE("\n -> ");
        Word P1 = PPERMV(r, LELTI(P, PO_POLY), R);
        IPWRITE(r, P1, LIST3(LFS("x1"), LFS("x2"), LFS("x3")));SWRITE("\n");
        LL = COMP(MPOLY(P1, NIL, NIL, PO_OTHER, PO_KEEP), LL);
    }

    return INV(LL);
}
