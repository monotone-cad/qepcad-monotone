/*======================================================================
                    F -< MONOTONE(FF)

Turns a quasi-affine CAD into a monotone one.

\Input
  \parm{DD} is a CAD (which has quasi-affine cells)
    \parm{r} is the space in which D lives

  Output
  \parm{D} is the CAD modified so all cells are monotone

======================================================================*/
#include "qepcad.h"

Word QepcadCls::MONOTONE(Word DD, Word r)
{
    Word T, F, C, I;
    printf("%d\n", r);
    if (r != 3) goto Return;

Step1: /* compute */
    LISTOFCWTV(DD, &T, &F);
    while (T != NIL) {
        ADV(T, &C, &T);
        I = LELTI(C, INDX);
        printf("%d\n", CELLDIM(C));

    }

    /*
     * loop through true cells
     * check it's (1,0,1)
     * take up to level 2 polynomials
     * find the "top" and "bottom", cell index (x,y,z+1) and (x,y,z-1), unless it's top or bottom of a stack
     * take level <= 2 polynomials, find first derivs = 0
     * check the degree - can they end up more than degree 1?
     * split cell, then figure out how to split all cells in the stack
     */

Return:
    return DD;
}
