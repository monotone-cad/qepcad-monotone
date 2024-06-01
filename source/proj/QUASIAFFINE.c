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

void ProcessPolynomials(Word r, Word Ps, Word *Ps1_, Word *PsQ_, Word *Js1_, Word *Fs1_)
{
    // we assume they are all r-variate polynomials, but only add those with nonzero degree in x_r as qepcad polynomials
    Word P, e, P1;
    Word Ps1 = NIL;
    Word PsQ = NIL; // extract factors only.
    Word Js1 = NIL; // zero degree xn x_r, as a polynomial in x_1,...,x_{r-1}
    while (Ps != NIL) {
        ADV(Ps, &P, &Ps);

        FIRST2(P, &e, &P1);
        if (e > 0) {
            Ps1 = COMP(MPOLY(P, NIL, NIL, PO_OTHER, PO_KEEP), Ps1);
        } else {
            Js1 = COMP(P1, Js1);
        }
    }

    // factorise
    Word Fs1 = IPLFAC(r, Ps1);

    // prepare to return Ps and Fs
    *Ps1_ = Ps1;
    *Fs1_ = Fs1;
    *Js1_ = Js1;

    // and prepare Qs
    while (Ps1 != NIL) {
        ADV(Ps1, &P, &Ps1);
        PsQ = COMP(LELTI(P, PO_POLY), Ps);
    }

    *PsQ_ = PsQ;
}

void QepcadCls::QUASIAFFINE(Word r, Word V, Word F, Word* A_, Word* P_, Word* J_)
{
    Word C, Ps, C1s, C2s, Js = NIL, PIs, PFs, PPs;
    SmoothOneTwoDim(r, V, F, &C1s, &C2s, &PIs, &PPs, &PFs);
    Word PF1 = CINV(PFs); // reverse order of proj factors to match cells

    // one-dimensional cells
    // i represents the projection onto one-dimensional coordinate subspace proj(i)
    // i.j represents proj(i,j), but we consider all minors with n-2 polynomials
    while (C1s != NIL) {
        ADV(C1s, &C, &C1s);
        Word Ps = ZeroPols(PF1, r, C);

        // proj(i)
        for (int i = 1; i <= r; ++i) {
            Word Is = GenerateIndex(r, i,0);
            Word J = JACOBI(r, NIL, 0, Ps, Is);

            if (!IPCONST(r, J)) {
                Js = COMP(J, Js);
            }
        }
    }

    // two-dimensional cells
    // i.j represent 2-dimensional coordinate subspace proj(i,j), and all minors for proj(i)
    while (C2s != NIL) {
        ADV(C2s, &C, &C2s);
        Word Ps = ZeroPols(PF1, r, C);

        // proj(i,j)
        for (int i = 1; i < r; ++i) {
            for (int j = i + 1; j <= r; j++) {
                Word Is = GenerateIndex(r, i,j);
                Word J = JACOBI(r, NIL, 0, Ps, Is);
                if (!IPCONST(r, J)) {
                    Js = COMP(J, Js);
                }

            }
        }
    }

    // CAD should be sign-invariant on Js. this means it should be compatible with their projections, too.
    // we add projections here to save recomputing the projection for entire CAD.

    // assign the pointers.
    *A_ = PIs; // input polynomials
    *P_ = PFs; // projection factors
    *J_ = PPs; // projection polynomials

    // add Js...
    // Ps1 contains QEPCAD polynomials, Fs1 contains its factors.
    Word Ps1, PsQ, Js1, Fs1;
    ProcessPolynomials(r, Js, &Ps1, &PsQ, &Js1, &Fs1);
    ADDPOLS(Ps1, r, LFS("Q"), J_);
    ADDPOLS(Fs1, r, LFS("Q"), P_);

    // and its projections
    while (r > 1) {
        LWRITE(PsQ); SWRITE(" ");
        LWRITE(Js1); SWRITE("\n");
        Js = ProjMcxUtil(r, PsQ);
        Js = CONC(Js, Js1); // polynomials with 0 degree in x_r
        --r;

        ProcessPolynomials(r, Js, &Ps1, &PsQ, &Js1, &Fs1);
        ADDPOLS(Ps1, r, LFS("Q"), J_);
        ADDPOLS(Fs1, r, LFS("Q"), P_);
    }
}

