/*======================================================================
                    A -< QUASIAFFINE(A, J, r; P, J)

Adding first derivatives of projections onto 1 and 2 dimensional
coordinate subspaces, to result in a CAD with all quasi-affine
cells.

\Input
  \parm{AA} projection factor structure. ([A_1, ..., A_r], with A_i being i-level proj factors
  \parm{AJ} projection polynomial (unfactored proj factor) structure. ([A_1, ..., A_r], with A_i being i-level proj factors
  \parm{r} number of variables

Output
  \parm{A}: modified projection factors
  \parm{J}: modified projection polys

SideEffect
  \parm{AA} is modified.
  \parm{AJ} is modified.

======================================================================*/
#include "qepcad.h"

// all (non-constant) partials of P (level k), factorised and formatted ready to append
Word PARTIALS(Word P, Word k);

// wrapper for IPLFAC to unset the parents. for when you don't want them referenced
Word LFAC(Word k, Word L);

void QepcadCls::QUASIAFFINE(Word A, Word J, Word r, Word* A_, Word* J_)
{
    Word AA, A1, A11, i, L;

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

Step2: /* dim >= 3: derivs of projections to 1 and  2d subspaces (add projections) */
    SWRITE("# adding projections to (other [than x1,x2]) 2d subspaces.\n");

    AA = PermuteAndProject(A, r, 2);
    printf("%d\n", LENGTH(AA));
    L = NIL;
    while (AA != NIL) {
        ADV(AA, &A11, &AA);

        L = CONC(L, PARTIALS(A11, r));
    }
    A = APPEND(A, r, LFAC(r, L));

Step3: /* dim >= 2: derivs of projection factors */
    SWRITE("# adding derivs of projections to 1d subspaces.\n");

    /* levels. */
    AA = A; i = 0;
    L = NIL;
    while (AA != NIL) {
        ADV(AA, &A1, &AA);
        i++;

        /* polynomials. */
        while (A1 != NIL) {
            ADV(A1, &A11, &A1);

            L = CONC(L, PARTIALS(A11, i));
        } /* END polynomials. */
        A = APPEND(A, i, LFAC(i, L));
    } /* END level. */

Return: /* prepare for return */
    *A_ = A;
    *J_ = J;
}

Word PARTIALS(Word P, Word k)
{
    Word L = NIL;
    Word Pp = IPALLPARTIALS(k, LELTI(P, PO_POLY), 1, 1);
    IPWRITE(k, LELTI(P, PO_POLY), LIST3(LFS("x"), LFS("y"), LFS("z")));SWRITE("\n");

    Word D = NIL;
    /* derivatives. */
    while (Pp != NIL) {
        ADV(Pp, &D, &Pp);

        if (IPCONST(k, D)) continue;

        IPWRITE(k, D, LIST3(LFS("x"), LFS("y"), LFS("z")));SWRITE("\n");
        L = COMP(MPOLY(D, NIL, NIL, PO_OTHER, PO_KEEP), L);
    } /* END derivatives. */

    return L;
}

Word LFAC(Word k, Word L)
{
    Word LL, P;

    L = IPLFAC(k, L);

    LL = L;
    while (L != NIL) {
        ADV(L, &P, &L);
        SLELTI(P, PO_PARENT, NIL);
    }

    return LL;
}
