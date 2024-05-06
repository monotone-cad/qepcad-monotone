// utility function to perform McCallum Projection
// use proj McCallum without leading coeffs, removes constants, see PROJMCx
// resultants and discriminants of input polynomials
// it is assumed that no polynomial in the list is equal to 0. otherwise the function will crash
// A : a set of polynomials in r variables
// r : positive integer
// return : projection onto R^{r-1} by McCallum's operator
#include "qepcad.h"

Word ProjMcxUtil(Word r, Word A)
{
    Word A1,Ap,Ap1,Ap2,App,D,L,Lh,P,R,W,i,t;

    // set of polynomials to project on this round
    A1 = NIL;

    // set of projected polynomials
    P = NIL;

    // construct the list of polynomials to project now and the list to project later
    while (A != NIL) {
        ADV(A, &Ap1, &A);

        // Ap1 = x_r^e Aq1
        Word e, Aq1;
        FIRST2(Ap1, &e, &Aq1);

        // if nonzero degree in x_r
        if (e > 0) {
            // project on this round
            A1 = COMP(Ap1, A1);
        } else {
            // defer to  next round
            P = COMP(Aq1, P);
        }
    }

    // discriminants
    Ap = A1;
    while (Ap != NIL) {
        ADV(Ap,&Ap1,&Ap);

        if (PDEG(Ap1) < 2) continue;

        D = IPDSCRQE(r,Ap1);
        if (!IPCONST(r-1, D)) {
            P = COMP(D,P);
        }
    }

    // resultants
    Ap = A1;
    while (Ap != NIL) {
        ADV(Ap,&Ap1,&Ap);

        App = Ap;
        while (App != NIL) {
            ADV(App,&Ap2,&App);

            R = IPRESQE(r,Ap1,Ap2);
            if (!IPCONST(r-1, P)) {
                P = COMP(R,P);
            }
        }
    }

    return P;
}


