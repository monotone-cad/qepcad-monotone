/*======================================================================
                    F -< MONOTONE(FF)

Adds extra polynomials to projection factor set to ensure monotone cells will be produced
Does not recompute the cad

\Input
  \parm{A} is a list a_1, ..., a_r where a_ is a list of i-level projection factors
  \parm{AA} is a list a_1, ..., a_r where a_ is a list of i-level input polynomials
  \parm{D} is a CAD
  \parm{r} is the space in which D lives

======================================================================*/
#include "qepcad.h"

// projection factor index of first zero (if not, then -1)
int PFIDXOFZERO(Word S, Word i, Word r)
{
    Word i1 = r + 1 - i;
    Word S1 = LELTI(S, i1);

    Word s = NIL;
    int idx = 0;
    while (S1 != NIL) {
        ADV(S1, &s, &S1);
        idx++;

        if (s == 0) {
            return idx;
        }
    }

    return -1;
}

// top and bottom of cell in decomposition D with index x, at level j.
// sets to c_top and c_bottom. in case of error, return 0, otherwise return 1
int FINDTOPANDBOTTOM(Word D, Word x, Word j, Word *C_top, Word *C_bottom)
{
    Word t1, t2;
    Word x1 = LCOPY(x);
    Word xj = LELTI(x, j);

    // top: index level j + 1
    SLELTI(x1, j, xj + 1);
    CELLFIDX(x1, D, C_top, &t1);

    // bottom: index level j - 1
    SLELTI(x1, j, xj - 1);
    CELLFIDX(x1, D, C_bottom, &t2);

    if (t1 == 0) {
        *C_top = NIL;
    }
    if (t2 == 0) {
        *C_bottom = NIL;
    }

    return t1 || t2;
}

// total derivative of P with respect to variable i

Word IPTDER(Word r, Word P, Word i)
{
    Word D = IPDER(r, P, i);

    // representation of D0 contains variable r (degree 0), this will remove it
    Word D1 = NIL, A = NIL, e = NIL;
    bool error = false; // indicates whether the derivative has positive degree in variable r
    while (D != NIL) {
        ADV2(D, &e, &A, &D);
        if (e > 0) error = true;

        if (D1 != NIL) D1 = IPSUM(r-1, D1, A);
        else D1 = A;
    }

    if (!error) return D1;

    // in case of error (?)
    Word D0 = IPDER(r, P, r);
    if (!IPCONST(r, D0)) {
        // exclude any roots of D0, these will be the points with "infinite gradient", i.e., when the line is vertical.
        // it should never happen because of quasi-affine cells (splitting based on critical points)
        // TODO
        printf("Non-constant denominator, and non-zero degree in variable r of numerator\n");
    }

    return D1;
}

// returns a polynomial in v_1 indicating critical points on the Surface along the Curve.
Word LAGRANGE(Word Surface, Word Curve)
{
    Word Vs = LIST3(LFS("x"), LFS("y"), LFS("z"));
    Word f = LELTI(Surface, PO_POLY);
    Word g = LELTI(Curve, PO_POLY);

    // compute jacobian matrix of f, g with v1 and v2 (we find derivs with total derivs)
    // solve (f'x * g'y) - (f'y * g'x) = 0 AND g = 0
    // TODO total derivs
    Word f1 = IPTDER(3, f, 1);
    Word f2 = IPTDER(3, f, 2);

    Word g1 = IPDER(2, g, 1);
    Word g2 = IPDER(2, g, 2);
    Word jacobiDet = IPSUM(
        2,
        IPPROD(2, f1, g2),
        IPNEG(2, IPPROD(2, f2, g1))
    );

    IPWRITE(2, jacobiDet, Vs); SWRITE("\n");
    IPWRITE(2, g, Vs); SWRITE("\n");
    // check if we can compute resultant - i.e., they both have degree > 0 in variable v2
    if (PDEG(jacobiDet) == 0) return NIL;

    // since we have 2 polynomials in 2 variables, we can solve for x by computing the resultant
    Word res = IPRESQE(2, jacobiDet, g);

    return MPOLY(res, NIL, NIL, PO_OTHER, PO_KEEP);
}

void QepcadCls::MONOTONE(Word A, Word J, Word D, Word r)
{
    Word Ct, Cf, C, I, i ,j, A1, P, P1, L;

Step1: /* calculate: looping through true cells with CELLDIM == 2 */
    LISTOFCWTV(D, &Ct, &Cf);
    L = NIL;
    while (Ct != NIL) {
        ADV(Ct, &C, &Ct);
        I = LEVELIDX(C);

        // choose an action based on dimension of cell
        if (LENGTH(I) > 2) {
            SWRITE("*** ERROR doesn't handle cells of dimension > 2");

            return;
        } else if (LENGTH(I) < 2) {
            SWRITE("# cell dim < 2 - continue\n");

            continue;
        }

        /* We have a 2-dimensional cell */
        FIRST2(I, &i, &j);

        printf("# two dimensional index: %d %d\n", i, j);

        // algorithm:
        // given our 2d cell, find the curves defining it's tob and bottom.
        // the cell will "fold over" at any critical points along these curves
        // use the lagrange multiplier method to find the critical points of curve (top / bottom of cell)
        // on the surface (function defining the graph, section cell)

        // find g_T and g_B - functions defining the top and bottom of 2d cell X
        Word C_top, C_bottom;
        A1 = LELTI(A, j);

        if (!FINDTOPANDBOTTOM(D, LELTI(C, INDX), j, &C_top, &C_bottom)) continue;
        int k = r; // TODO for r > 3
        // finding f (surface)
        Word S = LELTI(C, SIGNPF);
        int Ik = PFIDXOFZERO(S, k, r);
        P = LELTI(LELTI(A, k), Ik);

        // critical points of curve along the top of the cell
        if (C_top != NIL) {
            Word S_top = LELTI(C_top, SIGNPF);
            int Ij_top = PFIDXOFZERO(S_top, j, r);

            P1 = LAGRANGE(P, LELTI(A1, Ij_top));
            if (P1 != NIL) L = COMP(P1, L);
        }

        // same for the bottom
        if (C_bottom != NIL) {
            Word S_bottom = LELTI(C_bottom, SIGNPF);
            int Ij_bottom = PFIDXOFZERO(S_bottom, j, r);

            P1 = LAGRANGE(P, LELTI(A1, Ij_bottom));
            if (P1 != NIL) L = COMP(P1, L);
        }
    }

    // construct and append list of new polynomials from L.
    printf("length %d\n", LENGTH(L));
    ADDPOLS(L, 1, LFS("D"), &J);
    A = APPEND(A, 1, IPLFAC(1, L));
}

