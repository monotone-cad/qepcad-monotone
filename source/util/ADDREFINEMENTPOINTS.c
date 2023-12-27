/*======================================================================
  ADDREFINEMENTPOINTS(I, S, R1s; A, J, RPs)

Add refinement points based on real roots of a list of polynomials.

\Input
    I : (i_1,...,i_k), cell index
    S : (primitive?!) sample point
  R1s : list of univariate polynomials with integer coefficients, defining the refinement of subcad above cell with
        index I
Output
   As : projection factors structure
   Js : input polynomials structure
  RPs : Refinement points structure

======================================================================*/
#include "qepcad.h"

void ADDREFINEMENTPOINTS(Word I, Word S, Word R1s, Word* A_, Word* J_, Word* RPs_)
{
    const Word Z1 = LFS("K");
    const Word Z2 = LFS("M");

    Word k = LENGTH(I); // level of base cell
    Word k1 = k + 1; // level of current polynomials
    printf("adding refinement points for ");LWRITE(I);SWRITE("\n");

    // factorise the list of polynomials, add to J and store as rationals for root finding
    Word Ps = NIL;
    while (R1s != NIL) {
        Word P, s, c, L, Label;
        ADV(R1s, &P, &R1s);
        if (P == 0) continue;

        ADDPOL(PPREPVS(P, k), NIL, k1, Z1, J_, &Label);
        Label = COMP(Z1, Label);
        Word W = MPOLY(P, Label, NIL, PO_POLY, PO_KEEP);

        IPFACDB(1,P,&s,&c,&L);
        while (L != NIL) {
            Word P1, e, Q;
            ADV(L,&P1,&L);
            FIRST2(P1,&e,&Q);

            ADDPOL(PPREPVS(Q, k), LIST1(LIST3(PO_FAC,e,W)), k1, Z2, A_, &Label);
            Ps = COMP(RPFIP(1, Q), Ps);
        }
    }

    if (Ps == NIL) {
        return;
    }

    // compute the list of roots of the polynomials
    Word B = ROOTS(Ps);

    // now construct level k+1 sample points from the basis
    // TODO is S always primitive? should be since there are other coordinates above it
    Word SM, SI, Sb;
    FIRST3(S, &SM, &SI, &Sb);
    while (B != NIL) {
        Word M, J, a, bp, S1;
        ADV2(B, &J, &M, &B);
        Word b = LCOPY(Sb);

        // construct sample point using minimal polynomial M and isolating interval J
        if (PDEG(M) == 1) { // M is linear, easy!
            a = AFFRN(IUPRLP(M));
            bp = CONC(b,LIST1(a));
            S1 = LIST3(SM,SI,bp);
        } else if (PDEG(SM) == 1 || EQUAL(M, SM)) { // nonlinear, construct as an algebraic -- same M and I as in S
            a = AFGEN();
            bp = CONC(b,LIST1(a));

            S1 = LIST3(M,J,bp);
        } else { // algebraic, but with a new M and I. extended form with normalised polynomial.
            S1 = LIST5(AFPFIP(1, M), J ,SM, SI, b);
        }

        // add the sample points to set RPs
        Word RP1 = LELTI(*RPs_, k1);
        Word W = MPOLY(LIST2(S1, J), LIST3(Z2,k1,LENGTH(RP1) + 1), NIL, PO_REFINEMENT, PO_KEEP);
        SLELTI(W, PO_REFINEMENT, I);
        SLELTI(*RPs_, k1, COMP(W, RP1));
    }
}

