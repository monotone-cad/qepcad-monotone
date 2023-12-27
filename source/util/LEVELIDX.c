/*======================================================================
                    L <- LEVELIDX(C)

Level index -- lest of coordinate indices at which cell C is onedimensional

\Input
  \parm{C} cad cell
Output
  \parm{L} list [i_1, ..., i_m] of indices it which cell C is one-dimensional

======================================================================*/
#include "qepcad.h"

Word LEVELIDX(Word C)
{
    Word I, k, j, L;

    L = NIL;
    k = 0;
    I = LELTI(C, INDX);
    while (I != NIL) {
        ADV(I, &j, &I);
        k++;

        if (ODD(j)) {
            L = COMP(k, L);
        }
    }

    return INV(L);
}

