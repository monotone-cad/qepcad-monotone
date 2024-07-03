/*======================================================================
                      D <- FRONTIER(r,C,P)

Frontier condition (lifting above bad points, lazard-style)

\Input
  \parm{r} positive integer, dimension
  \parm{C} is the original cad with missing signs of proj factors (may be recursive i.e., level of C may be > 1)
  \parm{P} is the list~$(P_1,\ldots,P_r)$,
           where $P_i$ is the list of
           $i$--level projection factors.

Output
  \parm{D} is a truth--invariant partial CAD of $f$--space
           in which every leaf cell has a determined truth value.
======================================================================*/
#include "qepcad.h"
#include <iostream>

Word IndexOfFirstZero(Word S1)
{
    Word i = 1;
    while (S1 != NIL) {
        Word s;
        ADV(S1, &s, &S1);
        if (s == 0) {
            return i;
        }

        ++i;
    }

    return -1;
}

// returns (j1,...,jk) where each ji is the index of first zero for lezel-i projection factors.
Word SignatureIndex(Word S)
{
    Word J, S1;

    J = NIL;
    while (S != NIL) {
        ADV(S, &S1, &S);
        J = COMP(IndexOfFirstZero(S1), J);
    }

    return J;
}

// unique comp, with == check, O(n)
Word UCOMP(Word x, Word L)
{
    Word LL = L;
    while (LL != NIL) {
        Word y;
        ADV(LL, &y, &LL);

        if (x == y) return L;
    }

    return COMP(x, L);
}

void ProcessBadCells(Word r, Word C, Word As, Word i, Word j, Word S)
{
    if (C == NIL) return; // base case, nothing to do

    Word s;
    ADV(S, &s, &S);

    Word Ch = LELTI(C, CHILD);
    if (Ch == NIL) return;

    Word sample = LELTI(C, SAMPLE);
    Ch = RED(Ch);
    bool section = false;
    while (Ch != NIL) {
        Word C1, level, s1, SC1;
        section = !section;
        ADV(Ch, &C1, &Ch);
        level = LELTI(C1, LEVEL);
        SC1 = LELTI(C1, SIGNPF);
        s1 = IndexOfFirstZero(FIRST(SC1));

        if (!section && level > j && s == s1) {
            printf("possible bad cell "); LWRITE(LELTI(C1, INDX)); SWRITE("\n");
            Word Rp = LazardLifting(
                level,
                sample,
                As,
                COMP(s1, SignatureIndex(RED(SC1))),
                i,
                j
            );
            LWRITE(Rp); SWRITE("\n");
        }

        if (section && (level == j || s == s1)) {
            ProcessBadCells(r, C1, As, i, j, S);
        }
    }
}

Word QepcadCls::FRONTIER(Word r, Word k, Word D, Word As, Word* A_, Word* J_, Word* RPs_)
{
    Word Ch, TrueCells, junk, C1, C1_B, C1_T, C;
    Ch = LELTI(D, CHILD);

    // if r < 3, frontier condition is obtaoined automatically, if no children then nothing to do.
    if (r < 3 || Ch == NIL) return D;

    ADV(Ch, &C1, &Ch);

    // only one sector, bad cells are not possible.
    if (Ch != NIL);

    C1_B = NIL, C1_T = NIL;
    while (Ch != NIL) {
        // Ch = (top, next sector, ...)
        // cells will be taken in pairs.
        ADV(Ch, &C1_T, &Ch);
        FRONTIER(r, k+1, C1_T, As, A_, J_, RPs_);

        // get true (1,...)-cells
        LISTOFCWTV(C1, &TrueCells, &junk);
        while (TrueCells != NIL) {
            Word d, j, SI;
            ADV(TrueCells, &C, &TrueCells);

            d = TwoDimIndex(LELTI(C, INDX), &junk, &j);
            if (d != 2 || j == r) continue;

            printf("- true two-dimensional section cell "); LWRITE(LELTI(C, INDX)); SWRITE("\n");

            // find indices of polynomials which are zero on C.
            SI = REDI(SignatureIndex(LELTI(C, SIGNPF)), k);
            LWRITE(SI); SWRITE("\n");

            ProcessBadCells(r, C1_B, As, k, j, SI);
            ProcessBadCells(r, C1_T, As, k, j, SI);
        }

        // next sector.
        ADV(Ch, &C1, &Ch);
        C1_B = C1_T; // one sector's top is the next sector's bottom.
    }

    return D;
}

