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
// sets to c_top and c_bottom. both are NIL in case of error
int FINDTOPANDBOTTOM(Word D, Word x, Word j, Word *C_top, Word *C_bottom)
{
    Word t1, t2;
    Word x1 = LCOPY(x);
    Word xj = LELTI(x, j);

    // top
    SLELTI(x1, j, xj + 1);
    CELLFIDX(x1, D, C_top, &t1);

    // bottom
    SLELTI(x1, j, xj - 1);
    CELLFIDX(x1, D, C_bottom, &t2);

    if (t1 == 0 || t2 == 0) {
        printf("ERROR: cannot find top and bottom cells.\n");

        return 0;
    }

    return 1;
}

// total derivative of P with respect to variable i
Word IPTDER(Word r, Word P, Word i)
{
    Word D = IPDER(r, P, i);
    Word D1 = NIL; Word A = NIL; Word e = NIL;
    while (D != NIL) {
        ADV2(D, &e, &A, &D);
        if (D1 != NIL) D1 = IPSUM(r-1, D1, A);
        else D1 = A;
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
    SWRITE("Solve:\n");
    IPWRITE(2, jacobiDet, Vs); SWRITE("\n");
    IPWRITE(2, g, Vs); SWRITE("\n");
    // check if we can compute resultant - i.e., they both have degree > 0 in variable v2
    if (PDEG(jacobiDet) == 0) return NIL;

    Word res = IPRESQE(2, jacobiDet, g);
    IPWRITE(1, res, Vs); SWRITE("\n");


    return MPOLY(res, NIL, NIL, PO_OTHER, PO_KEEP);
}

void QepcadCls::MONOTONE(Word A, Word J, Word D, Word r)
{
    Word Ct, Cf, C, I, i ,j, A1, P, P1, L;

Step1: /* calculate: looping through true cells */
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
            printf("# cell dim < 2 - continue\n");

            continue;
        }

        /* We have a 2-dimensional cell */
        FIRST2(I, &i, &j);

        printf("# two dimensional index: %d %d\n", i, j);
        // ALGORITHM
        // find top and bottom of 2d part (part with i,j)
        // X := (I_i, I_j), X_B = (I_i, I_{j-1}), X_T = (I_i, I_{j+1})
        // X_T = g_T(x,y) = 0, X_B = g_B(x,y) = 0
        // C := f(X)
        // C_T := f(X_T), C_B := f(X_B)
        // f := C -> R | f(C) = 0
        // use lagrange method to find critical points of f along curve G_T and g_B

        // find g_T and g_B - functions defining the top and bottom of 2d cell X
        Word C_top, C_bottom;
        A1 = LELTI(A, j);

        if (FINDTOPANDBOTTOM(D, LELTI(C, INDX), j, &C_top, &C_bottom) == 0) return;

        Word S_top = LELTI(C_top, SIGNPF);
        int Ij_top = PFIDXOFZERO(S_top, j, r);

        Word S_bottom = LELTI(C_bottom, SIGNPF);
        int Ij_bottom = PFIDXOFZERO(S_bottom, j, r);

        int k = r; // TODO for r > 3
        Word S = LELTI(C, SIGNPF);
        int Ik = PFIDXOFZERO(S, k, r);

        P = LELTI(LELTI(A, k), Ik);
        P1 = LAGRANGE(P, LELTI(A1, Ij_top));
        if (P1 != NIL) L = COMP(P1, L);

        P1 = LAGRANGE(P, LELTI(A1, Ij_bottom));
        if (P1 != NIL) L = COMP(P1, L);

        printf("top: (%d,%d), bottom: (%d,%d), k: (%d,%d)\n", j, Ij_top, j, Ij_bottom, k, Ik);
    }

    // construct and append list of new polynomials from L.
    printf("length %d\n", LENGTH(L));
    ADDPOLS(L, 1, LFS("D"), &J);
    A = APPEND(A, 1, IPLFAC(1, L));
}

