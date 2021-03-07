/*======================================================================
                      D <- TICAD(Q,F,f,P,A)

Truth Invariant CAD Construction Phase. Recompute based on new (linear) proj factors.

\Input
  \parm{DD} is the original cad with missing sings of proj factors
  \parm{Q} is the list of quantifiers
           in the input formula.
  \parm{F} $=F(x_1,\ldots,x_r)$
           is the normalized input quantifier-free formula.
  \parm{f} is the number of free variables in the input formula.
  \parm{P} is the list~$(P_1,\ldots,P_r)$,
           where $P_i$ is the list of
           $i$--level projection factors.
  \parm{A} is the list~$(A_1,\ldots,A_r)$,
           where $A_i$ is the list of all
           the distinct $i$--level normalized input polynomials.

Output
  \parm{D} is a truth--invariant partial CAD of $f$--space
           in which every leaf cell has a determined truth value.
======================================================================*/
#include "qepcad.h"

Word QepcadCls::RECOMPUTE(Word DD, Word Q, Word F, Word f, Word P, Word A)
{
    return DD;
}

