/*======================================================================
                      ADDPOL(P,PP,k, Z;A,L)

Add (unfactored) polynomial to "A" with given label prefix Z

\Input
  P : a k-variate saclib polynomial
 PP : parent of P (or null if input)
  k : the level of P (as well as its "variate-ness")
  A : the list of polynomials extracted from the input formula thus far
  Z : Label prefix

\SideEffect
  A : is modified if P does not already appear in it
  L : L is set to the label (i.e. pair (i,j) ) of P in A
======================================================================*/
#include "qepcad.h"

void ADDPOL(Word P, Word PP, BDigit k, Word Z, Word *A_, Word *L_)
{
    Word A,A_k,As_k,n,i,As_k_i,Ap,L,Q,Ap_k;

Step1: /* Initialize */
    A = *A_;
    A_k = LELTI(A,k);

Step2: /* Search A for P */
    As_k = A_k; n = LENGTH(A_k);
    for (i = 1; i <= n; i++)
    {
        ADV(As_k,&As_k_i,&As_k);
        if (EQUAL(P,LELTI(As_k_i,PO_POLY)))
        {
            L = LIST2(k,i);
            Ap = A;
            goto Return;
        }
    }

Step3: /* P not found in A */
    L = LIST2(k,n + 1);
    Q = MPOLY(P,LIST3(Z,k,n+1),PP,PO_OTHER,PO_KEEP);
    Ap_k = SUFFIX(A_k,Q);
    SLELTI(A,k,Ap_k);
    Ap = A;
    goto Return;

Return: /* Prepare to return */
    *A_ = Ap;
    *L_ = L;

    return;
}
