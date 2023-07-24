/*======================================================================
                    F -< MONOTONE(FF)

Adds extra polynomials to projection factor set to ensure monotone cells will be produced
Does not recompute the cad

\Input
  \parm{A} is a list a_1, ..., a_r where a_i is a list of i-level projection factors
  \parm{AA} is a list a_1, ..., a_r where a_ is a list of i-level input polynomials
  \parm{D} is a CAD
  \parm{r} is the space in which D lives

======================================================================*/
#include "qepcad.h"

// for retrieving one-dimensional components of two-dimensional cell indices
// let I = (0,...,0,i_j=1,0,...,0,i_k=1,0,...,0)
// if C is two-dimensional, then I is of this form.
// Then j and k are set accordingly, true is returned.
// Otherwise, false is returned, value of j and k are undefined.
bool TwoDimIndex(Word I, Word *j_, Word *k_)
{
    Word J,d,n,l;

    n = 0, d = 0;
    while (I != NIL) {
        ADV(I, &l, &I);
        ++n;

        if (ODD(l)) {
            ++d;
            J = COMP(n, J);
        }

        // only consider 2-dimensional cells
        if (d > 2) {
            return false;
        }
    }

    if (d < 2) {
        return false;
    }

    // assign j and k, and return
    FIRST2(J, k_, j_);
    return true;
}

// return the cell in list L with (partial) index I.
Word FindByIndex(Word L, Word I, Word j, Word k)
{
    Word J,C,C1;

    while (L != NIL) {
        ADV(L, &C, &L);
        J = LELTI(C, INDX);

        // no partial match
        if (FIRST(I) != LELTI(J, k)) continue;

        // found
        if (j == k) return C;

        // partial match, recursion
        return FindByIndex(LELTI(C, CHILD), RED(I), j, k + 1);
    }

    // not found by recursion
    return NIL;
}

// return a list (f_j,...,f_k), where f_i is a level j polynomial which is 0 on C.
Word ZeroPols(Word A, Word r, Word C, Word i1, Word i2)
{
    Word k, L, A1, P, S, S1, s;
    k = LELTI(C, LEVEL);
    S = REDI(LELTI(C, SIGNPF), k - i2);
    A = REDI(A, r-i2); // throw out polynomials of higher level
    printf("%d %d %d\n", k, LENGTH(S), LENGTH(A));
    L = NIL;

    while (A != NIL && S != NIL) {
        // only look for requested polynomials
        if (i2 < i1) break;
        printf("round %d ", i2);
        ADV(A, &A1, &A);
        ADV(S, &S1, &S);

        // check signs, lvl k.
        while (A1 != NIL && S1 != NIL) {
            ADV(A1, &P, &A1);
            ADV(S1, &s, &S1);
            LWRITE(LELTI(P, PO_POLY)); printf(" [%d] ", s); SWRITE(" ");

            if (s == ZERO) {
                L = COMP(LELTI(P, PO_POLY), L);
            }
        }
        SWRITE("\n");
        // what if none were found?
        --i2;
    }
    printf("%d\n", LENGTH(L));
    return L; // note that l is in ascending order.
}

void QepcadCls::MONOTONE(Word A, Word J, Word D, Word r)
{
    Word AA, Ds, TrueCells, junk, C, I, I1, Ij, Ik, C0, CT, CB, Gs, Fk, Fs, P, P1;

    // loop over true cells
    LISTOFCWTV(D, &TrueCells, &junk);
    Ds = LELTI(D, CHILD); // cells of D, for searching cells
    AA = INV(A); // in reverse order, to match SIGNPF. for later.

    while (TrueCells != NIL) {
        ADV(TrueCells, &C, &TrueCells);
        I = LELTI(C, INDX);

        // consider 2-dimensional cells.
        if (!TwoDimIndex(I, &Ij, &Ik)) continue;

        // TODO
        SWRITE("----------\n");
        LWRITE(I);
        printf(" 2d indx: (%d, %d)\n", Ij, Ik);

        // proj_j(C) is a 0-dimensional cell (c_1,...,c_{j-1})
        if (Ij > 1) {
            C0 = FindByIndex(Ds, I, Ij-1, 1);
        } else {
            C0 = D;
        }

        // top and bottom of proj_k(C) are one-dimensional sections.
        I1 = LCOPY(I);
        // top
        SLELTI(I1, Ik, LELTI(I, Ik) + 1);
        CT = FindByIndex(Ds, I1, Ik, 1);
        // bottom
        SLELTI(I1, Ik, LELTI(I, Ik) - 1);
        CB = FindByIndex(Ds, I1, Ik, 1);

        // TODO
        printf("sub-cad index "); LWRITE(LELTI(C0, INDX)); SWRITE("\n");
        printf("top C index "); LWRITE(LELTI(CT, INDX)); SWRITE("\n");
        printf("bottom C index "); LWRITE(LELTI(CB, INDX)); SWRITE("\n");

        // find polynomials:
        // Gs = g_2,...,g_{k-1} define proj_{k-1}(C).
        Gs = ZeroPols(AA, r, C, Ij + 1, Ik - 1);
        Word GS1 = Gs, P, i = Ij;
        while (GS1 != NIL) {
            ADV(GS1, &P, &GS1);
            LWRITE(P); SWRITE("\n");
            IPDWRITE(i, P, GVVL); SWRITE("\n");
            i++;
        }

        // Fk (f_{k,T}, f_{k,B}) are 0 on CT and CB respectively.

        // Fs = (f_{K+1},...,f_n) define C.

        // semi-monotone

        // monotone
    }
}

