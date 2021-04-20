/*======================================================================
                      ADDPOL(Ps,k,Z;A)

Add (unfactored) polynomials to "A" with given prefix Z

\Input
  P s: a list of k-variate saclib polynomial
  k : the level of P (as well as its "variate-ness")
  A : the list of polynomials extracted from the input formula thus far
  Z : Label prefix

\SideEffect
  A : is modified if P does not already appear in it
  The label of each polynomial in Ps is set to the label it was assigned in A
======================================================================*/
#include "qepcad.h"

void ADDPOLS(Word Ps, BDigit k, Word Z, Word *A_)
{
    Word P, P1, k1, Label;

    while (Ps != NIL) {
        ADV(Ps, &P, &Ps);
        PSIMREP(k, LELTI(P, PO_POLY), &k1, &P1);

        ADDPOL(P1, k1, Z, A_, &Label);

        Label = COMP(Z, Label);
        SLELTI(P, PO_LABEL, Label);
    }
}

