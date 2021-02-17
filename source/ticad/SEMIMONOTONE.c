/*======================================================================
                    F -< SEMIMONOTONE(FF)

Adds extra polynomials to projection factor set to ensure semi-monotone cells will be produced

\Input
  \parm{A} is a list a_1, ..., a_r where a_ is a list of i-level projection factors
  \parm{r} is the space in which D lives

  Output
  \parm{P} is a list p_1, ..., p_r where p_i is the list a_i with partial derivatives with respect to x_1 of a_i added

======================================================================*/
#include "qepcad.h"

Word QepcadCls::SEMIMONOTONE(Word A, Word r)
{
    printf("in semimonotone!!");
    return A;
}
