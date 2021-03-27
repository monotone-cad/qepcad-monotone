/*======================================================================
                    F -< QUASIAFFINE(FF)

Modify the input formula to allow CAD to construct quasi-affine cells
by adding critical and singular points.

\Input
  \parm{FF} is a normalised formula

  Output
  \parm{F} is F with added first derivatives $=F(x_1,\ldots,x_r)$

======================================================================*/
#include "qepcad.h"

Word QepcadCls::QUASIAFFINE(Word FF)
{
    Word F;

Step1: /* Initialise */

Step2: /* figure out what type of formula we have */
    if (FF == TRUE || FF == FALSE) { F = FF; goto Return; }

    if (FIRST(FF) == OROP || FIRST(FF) == ANDOP) {
        Word F1, F2, op;
        op = FIRST(FF);
        F = LIST1(op);
        FF = RED(FF);

        while (FF != NIL) {
            ADV(FF, &F1, &FF);
            F2 = QUASIAFFINE(F1);
            F = COMP(F2, F);
        }

        F = INV(F);
        goto Return;
    }

Step3: /* atomic formula */
    Word t, P, r, I, Pp, D, S;
    FIRST4(FF, &t, &P, &r, &I);
    // calculate first derivatives
    Pp = IPALLPARTIALS(r, P, 1, 1);

    F = LIST2(FF, OROP);
    S = LIST1(ANDOP);

    while (Pp != NIL) {
        ADV(Pp, &D, &Pp);
        if (IPCONST(r, D)) continue;

        F = COMP(LIST4(EQOP, D, r, NIL), F);
        S = COMP(LIST4(EQOP, D, r, NIL), S);
    }

    if (LENGTH(S) == 1) {
        // no derivatives - return original
        F = FF;
        goto Return;
    }
    if (LENGTH(S) > 2) {
        S = INV(S);
        F = COMP(S, F);
    }
    F = LIST3(ANDOP, FF, INV(F));

Return:
    return F;
}
