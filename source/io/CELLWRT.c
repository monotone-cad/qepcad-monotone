/*======================================================================
                      CELLWRT(c)

Cell Write (Tarski).

Write out the cell C as a tarski formula
======================================================================*/
#include "qepcad.h"

// return the cell in list L with (partial) index I.
// purpose is to find parent cells.
// L : list of cells, level k
// I : partial index to match, level j
// return : a list of cells C s.t. C has index I
Word FindInRcad(Word L, Word I, Word j, Word k);

void QepcadCls::CELLWRT(Word c)
{
    Word S,S1,k,t,i,r,D,I;
    /* hide t; */

    /* Heading. */
    k = LELTI(c,LEVEL);
    r = LENGTH(GVVL);

    I = LELTI(c, INDX);
    SWRITE("---------- Information about the cell "); LWRITE(I);

    /* Dimension at each level (which is why we don't use CELLDIM) */
    SWRITE(" (Dimension ");

    int d = 0, o = 0, j = 0;
    i = NIL;
    // write the dimension at each level, calculate d = CELLDIM(c)
    while (I != NIL) {
        // pretty printing - opening bracket for the first cell, comma for later cells except the last one
        SWRITE(i == NIL ? "(" : ",");
        j++;

        ADV(I, &i, &I);
        // odd indices are 1-dimensional
        o = ODD(i);
        d += o;

        IWRITE(o);
    }

    // cells with index length < k implicitly have one-dimensional cells for remaining indices
    while (j < r) {
        j++;
        d += 1;

        SWRITE(",1");
    }

    SWRITE(") (");
    IWRITE(d);

    SWRITE("))\n\n");
    SWRITE("Signs of projection factors ------------------------\n\n");

    S = LELTI(c,SIGNPF);
    CELLIPLLDWR(GVVL, GVPF, S, k); SWRITE("\n");

    /* Signs of Projection Factors. */
    // To ensure the CAD is projection definable, we will use the RCAD instead.
    Word D1, P1;
    if (GVTD == NIL) { // initialise
        SWRITE("*** Initialising the RCAD. ***\n\n");

        Word D0 = GVPC, P0 = LCOPY(GVPF), J0 = LCOPY(GVPJ);
        for(i = GVNFV - LENGTH(J0); i > 0; i--) {
            J0 = INV(COMP(NIL,INV(J0)));
        }

        STRIPPED_BIGLOOP(J0,P0,P0,D0,GVNFV,&P1,&D1, 1);

        // cache for later
        GVTD = LIST2(P1, D1);
    } else { // computed before, use cached.
        FIRST2(GVTD, &P1, &D1);
    }

    SWRITE("Signs of projection factors (for defining formula) -\n\n");
    Word Cs = FindInRcad(LELTI(D1, CHILD), LELTI(c, INDX), k, 1);
    Word ncs = LENGTH(Cs);

    if (ncs > 1) {
        SWRITE("This cell consists of ");
        IWRITE(ncs);
        SWRITE("cells in the RCAD.\n\n");
    }

    while (Cs != NIL) {
        Word c1;
        ADV(Cs, &c1, &Cs);

        SWRITE("Index in RCAD: "); LWRITE(LELTI(c1, INDX)); SWRITE("\n");
        S = LELTI(c1,SIGNPF);
        CELLIPLLDWR(GVVL, P1, S, k); SWRITE("\n");
    }

    /* Write out the sample point. */
    SWRITE("Sample point ----------------------------------------\n\n");
    SAMPLEWR(c);

    /* Finish. */
    SWRITE("\n----------------------------------------------------\n");

    return;
}

Word FindInRcad(Word L, Word I, Word j, Word k)
{
    Word J, C, C1, Cs;

    Cs = NIL;
    while (L != NIL) {
        ADV(L, &C, &L);

        if (LENGTH(C) < 11) {
            J = LELTI(C, INDX);
        } else { // look at the cell of original cad containing this RCell
            C1 = LELTI(C, 11); // 11: INCELL
            J = LELTI(C1, INDX);
        }

        // no partial match
        if (FIRST(I) != LELTI(J, k)) continue;

        // otherwise, partial match. if j == k then we're done
        if (j == k) {
            Cs = COMP(C, Cs);
        } else {
            // if j < k, then we have a partial but incomplete match, continue searching recursively
            Cs = CONC(FindInRcad(LELTI(C, CHILD), RED(I), j, k + 1), Cs);
        }
    }

    // no partial match.
    return Cs;
}

