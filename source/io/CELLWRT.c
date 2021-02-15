/*======================================================================
                      CELLWRT(c)

Cell Write (Tarski).

Write out the cell C as a tarski formula
======================================================================*/
#include "qepcad.h"

void QepcadCls::CELLWRT(Word c)
{
        Word S,S1,k,t,i,D,I;
        /* hide t; */

Step1: /* Heading. */
        k = LELTI(c,LEVEL);
        SWRITE("---------- Information about the cell "); LWRITE(LELTI(c,INDX));
        SWRITE(" ----------\n\n");

Step2: /* Dimension */
        SWRITE("Dimension ");

        int d = 0, o = 0;
        I = LELTI(c,INDX);
        i = NIL;
        while (I != NIL) {
            SWRITE(i == NIL ? "(" : ",");
            ADV(I, &i, &I);
            o = ODD(i);

            d += o;

            IWRITE(o);
        }
        SWRITE(") ");
        IWRITE(d);
        SWRITE("\n\n");

Step3: /* Signs of Projection Factors. */
        S = LELTI(c,SIGNPF);
        CELLIPLLDWR(GVVL, GVPF, S, k); SWRITE("\n"); goto Return;

Step4: /* Finish. */
        SWRITE("\n----------------------------------------------------\n");
        goto Return;

Return: /* Prepare for return. */
        return;
}

