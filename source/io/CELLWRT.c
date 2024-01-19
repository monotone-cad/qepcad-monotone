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

    /* Signs of Projection Factors. */
    // first check if the CAD is projection definable.
    Word t1, D1;
    if (GVTD == NIL) { // initialise
        SWRITE("*** Determining whether the CAD is projection definable. ***\n\n");

        if (DOPFSUFF(GVPF, LIST1(GVPC)) == NIL) {
            // not projection definable. construct EPSCAD
            t1 = 0;
        } else {
            t1 = 1;
        }
        D1 = NIL;

        GVTD = LIST2(t1, NIL);
    } else {
        FIRST2(GVTD, &t1, &D1);
    }

    if (t1 == 1) {
        S = LELTI(c,SIGNPF);
        CELLIPLLDWR(GVVL, GVPF, S, k); SWRITE("\n");
    } else {
        SWRITE("*** CAD is not projection definable, using ESPCAD. ***\n\n");
    }

    /* Write out the sample point. */
    SWRITE("\nSample point ----------------------------------------\n\n");
    SAMPLEWR(c);

    /* Finish. */
    SWRITE("\n----------------------------------------------------\n");

    return;
}

