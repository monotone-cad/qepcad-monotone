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

Step1: /* Heading. */
    k = LELTI(c,LEVEL);
    r = LENGTH(GVVL);

    SWRITE("---------- Information about the cell "); LWRITE(LELTI(c,INDX)); SWRITE("\n\n");

Step2: /* Dimension at each level (which is why we don't use CELLDIM) */
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

Step3: /* Signs of Projection Factors. */
    S = LELTI(c,SIGNPF);
    CELLIPLLDWR(GVVL, GVPF, S, k); SWRITE("\n"); goto Return;

Step4: /* Finish. */
    SWRITE("\n----------------------------------------------------\n");
    goto Return;

Return: /* Prepare for return. */
    return;
}

