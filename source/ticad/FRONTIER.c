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

void ProcessBadCells(Word r, Word C, Word As, Word i, Word j, Word S, Word *RefinedCells_, Word *A_, Word *J_, Word* RPs_)
{
    if (C == NIL) return; // base case, nothing to do

    Word s;
    ADV(S, &s, &S);

    Word Ch = LELTI(C, CHILD);
    if (Ch == NIL) return;

    Word C1 = NIL, C2 = NIL, JT = NIL, JB = NIL;
    Word sample = LELTI(C, SAMPLE);
    bool section = false;
    while (true) {
        Word level, s1, SC1;
        C1 = C2;
        JB = JT;
        section = !section;

        // loop exit check. reached end of list.
        if (Ch == NIL && C1 == NIL) break;

        if (Ch != NIL) {
            ADV(Ch, &C2, &Ch);
            level = LELTI(C2, LEVEL);
        } else {
            C2 = NIL;
            level = LELTI(C1, LEVEL);
        }

        if (!section && C2 != NIL) {
            Word SM, SJ;
            GETSAMPLEK(level, LELTI(C2, SAMPLE), &SM, &SJ);
            JT = RNQ(RNSUM(FIRST(SJ), SECOND(SJ)), RNINT(2));
        } else if (!section) {
            JT = NIL;
        }

        if (C1 == NIL) continue;

        SC1 = LELTI(C1, SIGNPF);
        s1 = IndexOfFirstZero(FIRST(SC1));

        Word Idx = LELTI(C1, INDX);
        Word I1x = COMP(s1, Idx); // bit like a hash, polynomial plus cell index to indicate that a cell has been refined

        // a bad cell is a (0,...,0,1)-cell of level greater than J, with matching sign which has not yet been refined
        if (!section && level > j && s == s1 && LSRCH(I1x, *RefinedCells_) == 0) {
            printf("possible bad cell, polynomial = %d, ", s1); LWRITE(LELTI(C1, INDX));
            JB == NIL ? SWRITE("-infty") : RNWRITE(JB); SWRITE(" ");
            JT == NIL ? SWRITE("infty") : RNWRITE(JT); SWRITE("\n");

            Word RP = LazardLifting(
                level,
                sample,
                As,
                COMP(s1, SignatureIndex(RED(SC1))),
                i,
                j
            );

            ADDREFINEMENTPOINTS(Idx, sample, RP, LIST2(JB, JT), A_, J_, RPs_);
            *RefinedCells_ = COMP(I1x, *RefinedCells_);
            LWRITE(RP); SWRITE("\n");
        }

        if (section && (level == j || s == s1)) {
            ProcessBadCells(r, C1, As, i, j, S, RefinedCells_, A_, J_, RPs_);
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

    Word RefinedCells = NIL;
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

            // find indices of polynomials which are zero on C.
            SI = REDI(SignatureIndex(LELTI(C, SIGNPF)), k);

            ProcessBadCells(r, C1_B, As, k, j, SI, &RefinedCells, A_, J_, RPs_);
            ProcessBadCells(r, C1_T, As, k, j, SI, &RefinedCells, A_, J_, RPs_);
        }

        // next sector.
        ADV(Ch, &C1, &Ch);
        C1_B = C1_T; // one sector's top is the next sector's bottom.
    }

    return D;
}

