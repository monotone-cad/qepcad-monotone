/*======================================================================
                   F -< MONOTONE(FF)

Adds extra polynomials to projection factor set to ensure monotone cells will be produced
A refinement of the kind { x_i < c }, { x_i = c }, { x_i > c } will be required.

\Input
  \parm{A} is a list a_1, ..., a_r where a_i is a list of i-level projection factors
  \parm{AA} is a list a_1, ..., a_r where a_ is a list of i-level input polynomials
  \parm{D} is a CAD
  \parm{r} is the space in which D lives

======================================================================*/
#include "qepcad.h"

// suppose I = (i_1,...,i_n) is a cell index.
// if i_j = i_k = 1 but all other components are 0, then C is 2-dimensional.
// set j and k to indicate which components are equal to one and return true.
// otherwise, C is not two-dimensional. return false. Values of j and k are not defined.
bool TwoDimIndex(Word I, Word *j_, Word *k_)
{
    Word J,d,n,l;

    // d tracks cell dimension.
    n = 0, d = 0;
    while (d <= 3 && I != NIL) {
        ADV(I, &l, &I);
        ++n;

        if (ODD(l)) {
            ++d;
            J = COMP(n, J);
        }
    }

    // if cell is two-dimensional
    if (d == 2) {
        FIRST2(J, k_, j_);

        return true;
    }

    // otherwise cell is not two-dimensional
    return false;
}

// return the cell in list L with (partial) index I.
// purpose is to find parent cells.
// L : list of cells, level k
// I : partial index, level j
// return : C s.t. C has index I
Word FindByIndex(Word L, Word I, Word j, Word k)
{
    Word J,C,C1;

    while (L != NIL) {
        ADV(L, &C, &L);
        J = LELTI(C, INDX);

        // no partial match
        if (FIRST(I) != LELTI(J, k)) continue;

        // otherwise, partial match. if j == k then we're done
        if (j == k) return C;

        // if j < k, then we have a partial match but need to keep on searching recursively
        return FindByIndex(LELTI(C, CHILD), RED(I), j, k + 1);
    }

    // no partial match.
    return NIL;
}

// return a list (f_j,...,f_k), where f_i is a level j polynomial which is 0 on C.
// A : reversed list of projection factors
// r : ambient dimension of CAD
// C : cell to check
// i1 <= i1 : min and max indices to search
// return : (f_{i_1}, ..., f_{i2})
Word ZeroPols(Word A, Word r, Word C, Word i1, Word i2)
{
    Word k, L, A1, P, Q, S, S1, s;
    k = LELTI(C, LEVEL);
    S = REDI(LELTI(C, SIGNPF), k - i2);
    A = REDI(A, r-i2); // throw out polynomials of higher level
    L = NIL;

    while (A != NIL && S != NIL) {
        // only look for requested polynomials
        if (i2 < i1) break;

        ADV(A, &A1, &A);
        ADV(S, &S1, &S);

        // check signs of level k polynomials
        Q = 0;
        while (A1 != NIL && S1 != NIL) {
            ADV(A1, &P, &A1);
            ADV(S1, &s, &S1);

            if (s == ZERO) {
                Q = LELTI(P, PO_POLY);

                break;
            }
        }

        // append polynomial.
        // note 1: if none found, then 0 is appended.
        // note 2: if more than one polynomial is zero on C, then only the first in the list is added
        L = COMP(Q, L);
        --i2;
    }

    return L; // note: L is in ascending level order.
}



// evaluate r-variate integer polynomial P at sample point S = (c_1,...,c_k)
// return Q in Z[x_{k+1},...,k_r] = P(c_1,...,c_k,x_{k+1},...,k_r)
// r : positive integer
// P : integer polynomial in Z[x_1,...,x_r]
// S : sample point (Q, I, (c_1,...,c_k)), in primitive form
// return : substituted polynomial in Z[x_{k+1},...,x_r]
Word Substitute(Word r, Word P, Word S)
{
    Word LL, P1, Q, J, c;

    // determine if S is rational or algebraic
    FIRST3(S, &Q, &J, &c);

    // sample subset R^0, nothing to do
    if (LENGTH(c) == 0) return P;
    LWRITE(P); SWRITE("\n");

    if (PDEG(Q) == 1) {
        LWRITE(RCFAFC(c)); SWRITE("\n");
        P1 = IPFRP(r - LENGTH(c), IPRNME(r, P, RCFAFC(c)));
    }
    LWRITE(P1); SWRITE("\n");

    return P1;
}

void QepcadCls::MONOTONE(Word A, Word J, Word D, Word r)
{
    Word AA, Ds, TrueCells, junk, C, I, I1, Ij, Ik, C0, S0, CT, CB, Gs, FT, FB, Fs, P, P1;

    // consider true cells
    LISTOFCWTV(D, &TrueCells, &junk);
    Ds = LELTI(D, CHILD); // cells of D, for searching cells
    AA = INV(A); // in reverse order, to match SIGNPF. for later.

    while (TrueCells != NIL) {
        ADV(TrueCells, &C, &TrueCells);
        I = LELTI(C, INDX);

        // only consider two-dimensional cells
        if (!TwoDimIndex(I, &Ij, &Ik)) continue;

        // TODO
        SWRITE("----------\n");
        LWRITE(I);
        printf(" 2d indx: (%d, %d)\n", Ij, Ik);

        // C0 := proj_{j-1}(C) is a 0-dimensional cell (c_1,...,c_{j-1})
        if (Ij > 1) {
            C0 = FindByIndex(Ds, I, Ij-1, 1);
        } else {
            C0 = D;
        }

        // sample point (c_1,...,c_{j-1}) of C0.
        S0 = LELTI(C0, SAMPLE);

        // top and bottom of proj_k(C) are one-dimensional sections by definition.
        I1 = LCOPY(I);
        // top: (i_1,...,i_k + 1)
        SLELTI(I1, Ik, LELTI(I, Ik) + 1);
        CT = FindByIndex(Ds, I1, Ik, 1);
        // bottom: (i_1,...,i_k - 1)
        SLELTI(I1, Ik, LELTI(I, Ik) - 1);
        CB = FindByIndex(Ds, I1, Ik, 1);

        // TODO
        printf("sub-cad level "); IWRITE(LELTI(C, LEVEL)); SWRITE(", ");
        printf("sub-cad index "); LWRITE(LELTI(C0, INDX)); SWRITE("\n");
        printf("top C index "); LWRITE(LELTI(CT, INDX)); SWRITE("\n");
        printf("bottom C index "); LWRITE(LELTI(CB, INDX)); SWRITE("\n");

        // find polynomials:
        // Gs = g_2,...,g_{k-1} define proj_{k-1}(C).
        Gs = ZeroPols(AA, r, C, Ij + 1, Ik - 1);
        Word LL = Gs, P, Q, i = Ij; // TODO
        while (LL != NIL) {
            ADV(LL, &P, &LL);
            IPDWRITE(i, P, GVVL); SWRITE("\n");
            Q = Substitute(i, P, S0);
            IPDWRITE(i+1-Ij, Q, GVVL); SWRITE("\n");

            i++;
        }

        // Fk (f_{k,T}, f_{k,B}) are 0 on CT and CB respectively.
        FT  = ZeroPols(AA, r, CT, Ik, Ik);
        FB  = ZeroPols(AA, r, CB, Ik, Ik);
        SWRITE("f_top := "); IPDWRITE(Ik, FIRST(FT), GVVL); SWRITE("\n");
        Q = Substitute(Ik, FIRST(FT), S0);
        SWRITE("  subst := "); IPDWRITE(Ik + 1 - Ij, Q, GVVL); SWRITE("\n");
        SWRITE("f_bottom := "); IPDWRITE(Ik, FIRST(FB), GVVL); SWRITE("\n");
        Q = Substitute(Ik, FIRST(FB), S0);
        SWRITE("  subst := "); IPDWRITE(Ik + 1 - Ij, Q, GVVL); SWRITE("\n");

        // Fs = (f_{K+1},...,f_n) is a map from proj_{k}(C) to R^{n-k}, of which C is the graph.
        Fs = ZeroPols(AA, r, C, Ik + 1, r);
        LL = Fs, i = Ik + 1;
        while (LL != NIL) {
            ADV(LL, &P, &LL);
            IPDWRITE(i, P, GVVL); SWRITE("\n");
            i++;
        }

        // semi-monotone -- Lagrange on f_T and f_B

        // monotone
    }
}

