/*======================================================================
                    F -< ADDFIRSTDERIVATIVES(FF)

Add disjunction of first derivatives of polynomials to formula.

\Input
  \parm{FF} is a normalised formula

  Output
  \parm{F} is F with added first derivatives $=F(x_1,\ldots,x_r)$

======================================================================*/
#include "qepcad.h"

Word QepcadCls::ADDFIRSTDERIVATIVES(Word FF)
{
    Word F;

Step1: /* Initialise */
    printf("In addfirstderivatives\n"); fflush(0);

Step2: /* figure out what type of formula we have */
    if (FF == TRUE || FF == FALSE) { F = FF; goto Return; }

    if (FIRST(FF) == OROP || FIRST(FF) == ANDOP) {
        Word F1, F2, op, Fcp;
        op = FIRST(FF);
        Fcp = LCOPY(FF);
        F = LIST1(OROP);
        FF = RED(FF);

        while (FF != NIL) {
            ADV(FF, &F1, &FF);
            F2 = ADDFIRSTDERIVATIVES(F1);
            F = CONC(F, RED(F2));

        }
        F = LIST3(op, Fcp, F);

        goto Return;
    }
Step3: /* atomic formula */
    Word t, P, r, I, Pp, D, S;
    FIRST4(FF, &t, &P, &r, &I);
    SWRITE("we have an atomic formula!");
    // calculate first derivatives
    Pp = IPALLPARTIALS(r, P, 1, 1);

    F = LIST1(FF);
    S = LIST1(FF);

    while (Pp != NIL) {
        ADV(Pp, &D, &Pp);
        if (IPCONST(r, D)) continue;

        F = COMP(LIST4(EQOP, D, r, NIL), F);
        S = COMP(LIST4(EQOP, D, r, NIL), S);
    }

    S = COMP(ANDOP, S);
    if (LENGTH(S) > 2) F = COMP(S, F);
    F = COMP(OROP, F);

Return:
  printf("about to return... %p", F);fflush(0);
    return F;

}
