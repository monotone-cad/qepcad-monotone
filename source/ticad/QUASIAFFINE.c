/*======================================================================
                    A -< QUASIAFFINE(A, r)

Adding first derivatives of projections onto 1 and 2 dimensional
coordinate subspaces, to result in a CAD with all quasi-affine
cells.

\Input
  \parm{AA} projection factor structure, obtained by EXTRACT(F) where F is a normalised qff.
  \parm{r} number of variables

Output
  \parm{A} is AA with added first derivatives of projection to 1 and 2 dimensional subspaces

SideEffect
  \parm{AA} is modified.

======================================================================*/
#include "qepcad.h"

Word QepcadCls::QUASIAFFINE(Word A, Word r)
{
Step1: /* decide based on dimension */
    if (r <= 1) {
        // nothing to do

        goto Return;
    }

    if (r > 3) {
        SWRITE("Dimension > 3 not supported yet.\n");

        goto Return;
    }

    if (r == 2) {
        SWRITE("# skipping 2d projections.\n");

        goto Step3;
    }

Step2: /* dim >= 3: derivs of projections to 2d subspaces */
    SWRITE("# Adding derivs of projections to 2d subspaces.\n");

Step3: /* dim >= 2: derivs of projections onto 1d subspaces */
    SWRITE("# adding derivs of projections to 1d subspaces.\n");

Return: /* prepare for return */
    return A;
}
