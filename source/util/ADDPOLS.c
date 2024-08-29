/*======================================================================
                      ADDPOL(Ps,k,Z;A)

Add (unfactored) polynomials to "A" with given prefix Z

\Input
 Ps : a list of k-variate saclib polynomial
  k : positive integer, level and number of variables in P
  A : set of [input/projection] polynomials.
  Z : Label prefix

\SideEffect
  A : if A does not contain each P in Ps, A is modified so that it contains P

======================================================================*/
#include "qepcad.h"

void ADDPOLS(Word Ps, BDigit k, Word Z, Word *A_)
{
    Word P, P1, k1, Label;

    while (Ps != NIL) {
        ADV(Ps, &P, &Ps);
        PSIMREP(k, LELTI(P, PO_POLY), &k1, &P1);

        ADDPOL(P1, LELTI(P, PO_PARENT), k1, Z, A_, &Label);

        Label = COMP(Z, Label);
        SLELTI(P, PO_LABEL, Label);
    }
}

