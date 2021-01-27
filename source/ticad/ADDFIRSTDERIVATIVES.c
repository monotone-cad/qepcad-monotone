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
    F = FF;
    SWRITE("In addFirstDerivatives.");

Step2: /* figure out what type of formula we have */
    if (F == TRUE || F == FALSE) goto Return;
    if (FIRST(F) == ANDOP) {
        // TODO
        SWRITE("Conjunctions not implemented yet!");
        goto Return;
    }

    if (FIRST(F) == OROP) {
        Word F1 = NIL;
     SWRITE("Disjunction");
        ADV(FF, &F1, &FF);
        while (FF != NULL) {
            F1 = ADDFIRSTDERIVATIVES(F1);

            ADV(FF, &F1, &FF);
            goto Return;
        }
    }

Step3: /* atomic formula */
    Word t, P, r, I, Pp, D, S;
    FIRST4(F, &t, &P, &r, &I);
    SWRITE("we have an atomic formula!");
    // calculate first derivatives
    Pp = IPALLPARTIALS(r, P, 1, 1);
    printf("number of derivs: %d\n",LENGTH(Pp));

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
    printf("list length: %d\n",LENGTH(F));

Return:
    return F;

}
