/*======================================================================
                      CELLWRT(c)

Cell Write (Tarski).

Write out the cell C as a tarski formula
======================================================================*/
#include "qepcad.h"

void QepcadCls::CELLWRT(Word c)
{
    Word S,S1,k,t,i,r,D,I;
    /* hide t; */

    /* Heading. */
    k = LELTI(c,LEVEL);
    r = LENGTH(GVVL);

    SWRITE("---------- Information about the cell "); LWRITE(LELTI(c,INDX)); SWRITE("\n\n");

    /* Dimension at each level (which is why we don't use CELLDIM) */
    SWRITE("Dimension ");

    int d = 0, o = 0, j = 0;
    I = LELTI(c, INDX);
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

    SWRITE(") ");
    IWRITE(d);

    SWRITE("----------\n\n");
    SWRITE("\nSigns of projection factors ------------------------\n\n");

    S = LELTI(c,SIGNPF);
    CELLIPLLDWR(GVVL, GVPF, S, k); SWRITE("\n");



    SWRITE("\nSigns of (guaranteed definable ) projection factors \n\n");

    /* Signs of Projection Factors. */
    // To ensure the CAD is projection definable, we will use the ESPCAD instead.
    Word D1, P1;
    if (GVTD == NIL) { // initialise
        SWRITE("*** Initialising the ESPCAD. ***\n\n");

        Word D0 = GVPC, P0 = LCOPY(GVPF), J0 = LCOPY(GVPJ);
        for(i = GVNFV - LENGTH(J0); i > 0; i--) {
            J0 = INV(COMP(NIL,INV(J0)));
        }

        STRIPPED_BIGLOOP(J0,P0,P0,D0,GVNFV,&P1,&D1);

        // cache for later
        GVTD = LIST2(P1, D1);
    } else { // computed before, use cached.
        FIRST2(GVTD, &P1, &D1);
    }

    // TODO find by index from D1.

    Word c1 = FindByIndex(LELTI(D1, CHILD), LELTI(c, INDX), k, 1);
    S = LELTI(c1,SIGNPF);
    CELLIPLLDWR(GVVL, P1, S, k); SWRITE("\n");

    /* Write out the sample point. */
    SWRITE("\nSample point ----------------------------------------\n\n");
    SAMPLEWR(c);

    /* Finish. */
    SWRITE("\n----------------------------------------------------\n");

    return;
}

