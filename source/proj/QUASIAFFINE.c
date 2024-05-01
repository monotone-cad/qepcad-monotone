/*======================================================================
  QUASIAFFINE(A, r, *A_);

considers the critical points of projections of the smooth two-dimensional locus of input set, onto one and two
dimensional coordinate subspaces.

\Input
  \parm{r} positive integer, number of variables
  \parm{V} variable list
  \parm{F} quantifier-free Boolean formula defining F.
  \parm{A} list of input polynomials extracted from F

Output
  \parm{*A}: set of projection factors, modified to include first partial derivatives of each element with respect to
  each vaciable.

SideEffect
  \parm{A} is also modified.

======================================================================*/
#include "qepcad.h"

// return the dimension of cell C
Word Dim(Word C)
{
    Word d = 0, j, I = LELTI(C, INDX);

    while (I != NIL) {
        ADV(I, &j, &I);

        if (ODD(j)) d++;
    }

    return d;
}

// return a list of CAD cells, which is the smooth two-dimensional locus of V (defined by formula F)
void SmoothOneTwoDim(Word r, Word V, Word F, Word *C1s_, Word *C2s_, Word *PIs_, Word *PPs_, Word *PFs_)
{
    Word C1s = NIL, C2s = NIL, C, Ct, Cf;

    QepcadCls Q;
	INITSYS();

    // set input formula
    Q.SETINPUTFORMULA(V,LIST4(r, r, NIL, F));
    Q.PRDQFF();
    Q.CADautoConst();

    // special case: trivially false
    if (Q.GVPC == 0) {
        *C1s_ = NIL, *C2s_ = NIL, *PPs_ = NIL, *PFs_ = NIL;

        return;
    }

    LISTOFCWTV(Q.GVPC, &Ct, &Cf);

    // collect all the 2d cells. Error out if the input set has dimension > 2.
    while (Ct != NIL) {
        ADV(Ct, &C, &Ct);

        Word d = Dim(C);

        if (d > 2) {
            FAIL("source/ticad/MONOTONE", "cell dimension greater than 2 not supported.");
        } else if (d == 1) {
            C1s = COMP(C, C1s);
        } else if (d == 2) {
            C2s = COMP(C, C2s);
        }
    }

    // prepare to return
    *C1s_ = C1s;
    *C2s_ = C2s;
    *PIs_ = Q.GVNIP;
    *PFs_ = Q.GVPF;
    *PPs_ = Q.GVPJ;
}

// list of polynomials (first in list) which have sign zero on cell C.
Word ZeroPols(Word A, Word r, Word C)
{
    Word i, L, A1, P, S, S1, s;
    S = LELTI(C, SIGNPF);
    L = NIL;
    i = 0;

    while (A != NIL && S != NIL) {
        ADV(A, &A1, &A);
        ADV(S, &S1, &S);

        // check signs of level k polynomials
        while (A1 != NIL && S1 != NIL) {
            ADV(A1, &P, &A1);
            ADV(S1, &s, &S1);

            if (s == ZERO) {
                Word Q = PADDVS(LELTI(P, PO_POLY), i);
                L = COMP(Q, L);

                break;
            }
        }

        ++i;
    }

    return L;
}

// returns list (1,...,n) without i,j.
Word GenerateIndex(Word r, Word i, Word j)
{
    Word L = NIL;
    while (r > 0) {
        if (r != i && r != j) {
            L = COMP(r, L);
        }

        --r;
    }

    return L;
}

void AddJacobis(Word r, Word Ps, Word *P_, Word *J_)
{
    // sets of level k jacobis
    Word K = NIL; for (int i=0; i<r; i++) K = COMP(NIL,K);

    // determine level of P.
    while (Ps != NIL) {
        Word PP, P1, e1;
        ADV(Ps, &P1, &Ps);

        Word k = r + 1;
        do {
            --k;
            PP = P1;

            FIRST2(PP, &e1, &P1);
        } while (e1 == 0);

        // append at correct level
        Word K1 = LELTI(K, k);
        K1 = COMP(MPOLY(PP, NIL, NIL, PO_OTHER, PO_KEEP), K1);
        SLELTI(K, k, K1);
    }

    // P is of level k. add to P and J
    Word j = 1;
    while (K != NIL) {
        Word R, K1;
        ADV(K, &K1, &K);

        printf("j = %d, K1 length = %d\n", j, LENGTH(K1));
        ADDPOLS(K1, j, LFS("Q"), J_);

        R = IPLFAC(j, K1);
        ADDPOLS(R, j, LFS("Q"), P_);
        ++j;
    }
}

void QepcadCls::QUASIAFFINE(Word r, Word V, Word F, Word* A_, Word* P_, Word* J_)
{
    Word C, Ps, C1s, C2s, Js = NIL, PIs, PFs, PPs;
    SmoothOneTwoDim(r, V, F, &C1s, &C2s, &PIs, &PPs, &PFs);
    Word PF1 = CINV(PFs); // reverse order of proj factors to match cells

    // one-dimensional cells
    SWRITE("1d cells \n");
    while (C1s != NIL) {
        ADV(C1s, &C, &C1s);
        Word Ps = ZeroPols(PF1, r, C);

        // n-1 polynomials Ps will appear in rows of the jacobi, now we generate the lists of indices
        for (int i = 1; i <= r; ++i) {
            Word Is = GenerateIndex(r, i,0);
            Word J = JACOBI(r, NIL, 0, Ps, Is);

            if (!IPCONST(r, J)) {
                LWRITE(Is); SWRITE(" "); LWRITE(J); SWRITE("\n");
                Js = COMP(J, Js);
            }
        }
    }

    // two-dimensional cells
    SWRITE("2d cells \n");
    while (C2s != NIL) {
        ADV(C2s, &C, &C2s);
        Word Ps = ZeroPols(PF1, r, C);

        // n-2 polynomials Ps will appear in rows of the jacobi, now we generate the lists of indices
        for (int i = 1; i < r; ++i) {
            for (int j = i + 1; j <= r; j++) {
                Word Is = GenerateIndex(r, i,j);
                Word J = JACOBI(r, NIL, 0, Ps, Is);
                if (!IPCONST(r, J)) {
                    LWRITE(Is); SWRITE(" "); LWRITE(J); SWRITE("\n");
                    Js = COMP(J, Js);
                }

            }
        }
    }

    // finally, add the jacobis
    *A_ = PIs;
    *P_ = PFs;
    *J_ = PPs;

    AddJacobis(r, Js, P_, A_);
}

