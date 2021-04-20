/*======================================================================
                      EXTRACTS(F,A; Ap)

Extract polynomials from a quantifier--free formula, Subalgorithm.

\Input
  \parm{F} is a non--constant normalized quantifier-free formula.
  \parm{A} is a list~$(A_1,\ldots,A_r)$ where each $A_i$ is
           the list of all the $i$--level polynomials found
           util the call to this algorithm was made.
\Output
  \parm{A'} is the list~$(A'_1,\ldots,A'_r)$ where each $A'_i$ is
           obtained by appending to $A_i$ all the $i$--level
           polynomials in $F$ if there are not there yet.

\SideEffect
   \parm{F} is modified in that each polynomial is assigned to
           an index.
   \parm{A} is modified.
======================================================================*/
#include "qepcad.h"

void EXTRACTS(Word F, Word A, Word *Ap_)
{
        Word T,P,k,I,L,j,Pp,Ap,Fp,F_i,X;

Step1: /* Non-constant normalized atomic formula. */
        if (!ISATOMF(F)) goto Step3;
	if (FIRST(F) == IROOT) goto Step2;
        FIRST4(F,&T,&P,&k,&I);
	if (k == 0) { Ap = A; goto Return; }
	ADDPOL(P,k,LFS("A"),&A,&L);
	SLELTI(F,4,L);
        Ap = A;
        goto Return;

Step2: /* Normalized atomic Extended Tarski formula */
	FIRST6(F,&X,&T,&j,&P,&k,&I);
	I = NIL;
	for(Pp = CINV(P); Pp != NIL; Pp = RED(Pp))
	{
	  ADDPOL(FIRST(Pp),k,LFS("A"),&A,&L);
	  I = COMP(L,I);
	}
	SLELTI(F,6,I);
	Ap = A;
	goto Return;

Step3: /* Other. */
        Ap = A; Fp = RED(F);
        while (Fp != NIL) { ADV(Fp,&F_i,&Fp); EXTRACTS(F_i,Ap,&Ap); }
        goto Return;

Return: /* Prepare for return. */
       *Ap_ = Ap;
       return;
}

