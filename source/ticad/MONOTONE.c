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
    if (c == NIL) return P;

    if (PDEG(Q) == 1) {
        P1 = IPFRP(r - LENGTH(c), IPRNME(r, P, RCFAFC(c)));
    }

    return P1;
}

// return a list (f_j,...,f_k), where f_i is a level j polynomial which is 0 on C with sample point S substituted in.
// A : reversed list of projection factors
// r : ambient dimension of CAD
// C : cell to check
// i1 <= i1 : min and max indices to search
// S : sample point of level l < i1 to substitute. give the empty sample point of D to do no substitution.
// return : (f_{i_1}, ..., f_{i2}) in Q[x_{l+1},...,k_n]
Word ZeroPolsSub(Word A, Word r, Word C, Word i1, Word i2, Word S0, Word l, Word n)
{
    Word k, L, A1, P, Q, S, S1, s;
    k = LELTI(C, LEVEL);
    S = REDI(LELTI(C, SIGNPF), k - i2);
    A = REDI(A, r - i2); // throw out polynomials of higher level
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
                // write substituted polynomial Q in Q[x_l+1,...,x_n]
                // P in Q[x_1,...,i_2]
                // Q in Q[x_l+1,...,i_2]
                Q = PADDVS(Substitute(i2, LELTI(P, PO_POLY), S0), l + n - i2);

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

// find the local maxima and minima of function f in Q[x_1,...,x_r] subject to constraints
// Gs in Q[k_1,...,x_i], 2 <= i < r using the method of Lagrange Multipliers.
// return a list of polynomials, each h in Q[x_1]
// such that the roots of h give the x_1-coordinates of the critical points of f.
Word LagrangeRefinement(Word r, Word f, Word i, Word Gs, Word Is)
{
    Word Q = JACOBI(r, f, i, Gs, Is);

    SWRITE("\njacobi det ");
    if (Q == 0) IWRITE(0); else LWRITE(Q); SWRITE("\n");

    return 0;
}

// iterative application of lagrange multipliers then projection to x1.
// r  : number of variables
// Gs : minimal constraints
// P  : polynomial denoting top / bottom of sector cell
// Fs : polynomials constituting map from sector cell to section cell
// lagrange is done on P and each element of Fs, adding more constraints each time.
Word Refinement(Word r, Word Gs, Word P, Word Fs)
{
    // TODO collect, project and return!
    Word Q, i, k, Is;
    // generate sequence Is = (1,...,k-1)
    i = 1;
    k = LENGTH(Gs) + 1;
    Is = NIL;
    while (i < k) {
        Is = COMP(i, Is);
        ++i;
    }

    if (Gs != NIL)
        LagrangeRefinement(r, P, k, Gs, Is);

    while (Fs != NIL) {
        ADV(Fs, &Q, &Fs);
        Gs = COMP(P, Gs);
        Is = COMP(k, Is);
        P = Q;
        ++k;

        LagrangeRefinement(r, Q, k, Gs, Is);
    }

    return NIL;
}

void QepcadCls::MONOTONE(Word A, Word J, Word D, Word r)
{
    Word AA, Ds, TrueCells, junk, C, I, I1, Ij, Ij1, nv, Ik, C0, S0, CT, CB, Gs, FT, FB, Fs, P, P1;

    // consider true cells
    LISTOFCWTV(D, &TrueCells, &junk);
    Ds = LELTI(D, CHILD); // cells of D, for searching cells
    AA = INV(A); // in reverse order, to match SIGNPF. for later.

    while (TrueCells != NIL) {
        ADV(TrueCells, &C, &TrueCells);
        I = LELTI(C, INDX);

        // only consider two-dimensional cells
        if (!TwoDimIndex(I, &Ij, &Ik)) continue;
        Ij1 = Ij - 1;
        nv = r - Ij1;

        // TODO
        SWRITE("----------\n");
        LWRITE(I);
        printf(" 2d indx: (%d, %d)\n", Ij, Ik);

        // C0 := proj_{j-1}(C) is a 0-dimensional cell (c_1,...,c_{j-1})
        if (Ij > 1) {
            C0 = FindByIndex(Ds, I, Ij1, 1);
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
        if (CT != NIL) { printf("top C index "); LWRITE(LELTI(CT, INDX)); SWRITE("\n"); }
        if (CB != NIL) { printf("bottom C index "); LWRITE(LELTI(CB, INDX)); SWRITE("\n"); }

        // find polynomials in sub-cad
        // TODO
        //   it takes a lot of work to find these polynomials, and work is often repeated.
        //   each label is unique, so we can check whether we have already substituted to save on effort.
        //   what other things could i cache? can i use the existing db faciility?
        // (note they will be in Z[x_i,...,x_l] after substitution):
        // Gs = g_2,...,g_{k-1} define proj_{k-1}(C).
        Gs = ZeroPolsSub(AA, r, C, Ij + 1, Ik - 1, S0, Ij1, nv);
        Word LL = Gs, P, Q; // TODO
        SWRITE("Gs:\n");
        while (LL != NIL) {
            ADV(LL, &P, &LL);
            IPDWRITE(nv, P, GVVL); SWRITE("\n");
        }

        // Fk (f_{k,T}, f_{k,B}) are 0 on CT and CB respectively.
        if (CT != NIL) FT = FIRST(ZeroPolsSub(AA, r, CT, Ik, Ik, S0, Ij1, nv));
        if (CB != NIL) FB = FIRST(ZeroPolsSub(AA, r, CB, Ik, Ik, S0, Ij1, nv));
        if (CT != NIL) { SWRITE("f_top := "); IPDWRITE(nv, FT, GVVL); SWRITE("\n"); }
        if (CB != NIL) { SWRITE("f_bottom := "); IPDWRITE(nv, FB, GVVL); SWRITE("\n"); }

        // Fs = (f_{K+1},...,f_n) is a map from proj_{k}(C) to R^{n-k}, of which C is the graph.
        Fs = ZeroPolsSub(AA, r, C, Ik + 1, r, S0, Ij1, nv);
        SWRITE("Fs:\n");
        LL = Fs;
        while (LL != NIL) {
            ADV(LL, &P, &LL);
            IPDWRITE(nv, P, GVVL); SWRITE("\n");
        }

        // semi-monotone -- critical points of top and bottom of sector below, if they exist
        // monotone
        if (CT != NIL) {
            SWRITE("top\n");
            Refinement(nv, Gs, FT, Fs);
        }

        if (CB != NIL) {
            SWRITE("bottom\n");
            Refinement(nv, Gs, FB, Fs);
        }
    }
}

