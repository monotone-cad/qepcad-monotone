/*======================================================================
                   F -< MONOTONE(FF)

Adds extra polynomials to projection factor set to ensure monotone cells will be produced
A refinement of the kind { x_i < c }, { x_i = c }, { x_i > c } will be required.

\Input
  \parm{A} proj factor set
  \parm{J} proj polynomial set
  \parm{r} is the space in which D lives
  \parm{D} is a CAD

Side Effect
  A and J are modified by adding new polynomials

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

// if polynomial deos not depend on variables x1,..,xk, then return a polynomial in x(k+1),...,xr
Word PNotDependVs(Word r, Word P, Word k)
{
    // base case r == 0, P is an integer
    if (r == 0) {
        return P;
    }

    Word P1, e, Q, Q1;

    // k < r: check for dependent variable
    if (r < k) {
        FIRST2(P, &e, &Q);

        if (e > 0) { // depends on x_i
            return NIL;
        } else {
            return PNotDependVs(r - 1, Q, k);
        }
    }

    // k >= r, reconstruct polynomial
    P1 = NIL;
    while (P != NIL) {
        ADV2(P, &e, &Q, &P);

        Q1 = PNotDependVs(r - 1, Q, k);

        // error case
        if (Q1 == NIL) return NIL;

        // otherwise construct the polynomial
        P1 = COMP2(Q1, e, P1);
    }

    return INV(P1);
}

// evaluate r-variate integer polynomial P at sample point S = (c_1,...,c_k)
// return Q in Z[x_{k+1},...,k_r] = P(c_1,...,c_k,x_{k+1},...,k_r)
// r : positive integer
// P : integer polynomial in Z[x_1,...,x_r]
// S : sample point (Q, I, (c_1,...,c_k)), in primitive form
// return : substituted polynomial in Z[x_{k+1},...,x_r]
Word Substitute(Word r, Word P, Word S)
{
    // determine if S is rational or algebraic
    Word Q, J, c;
    FIRST3(S, &Q, &J, &c);

    // sample subset R^0, nothing to do
    if (c == NIL) return P;

    // number of coordinates in sample point
    Word i0 = LENGTH(c);

    // if P does not depend on all x1,...,kk, then no substitution needed
    Word P1 = PNotDependVs(r, P, i0 + 1);
    if (P1 != NIL) return P1;

    // otherwise perform substitution
    Word k = r - i0;
    if (PDEG(Q) == 1) { // rational sample
        P1 = IPRNME(r, P, RCFAFC(c));
    } else { // algebraic sample
        P1 = IPAFME(r, Q, P, c);
        P1 = AFPICR(k, P1);
        printf("algebraic evaluation "); LWRITE(P1); SWRITE(", in Z: ");
        LWRITE(IPFRP(k, P1)); SWRITE("\n");
    }

    // convert to a polynomial with integer coefficients
    return IPFRP(k, P1);
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

// use proj McCallum without leading coeffs, removes constants, see PROJMCx
// resultants and discriminants of input polynomials
// it is assumed that no polynomial in the list is equal to 0. otherwise the function will crash
Word ProjMcx(Word r, Word A)
{
    Word A1,Ap,Ap1,Ap2,App,D,L,Lh,P,R,W,i,t;

    // set of polynomials to project on this round
    A1 = NIL;

    // set of projected polynomials
    P = NIL;

    // construct the list of polynomials to project now and the list to project later
    while (A != NIL) {
        ADV(A, &Ap1, &A);

        // Ap1 = x_r^e Aq1
        Word e, Aq1;
        FIRST2(Ap1, &e, &Aq1);

        // if nonzero degree in x_r
        if (e > 0) {
            // project on this round
            A1 = COMP(Ap1, A1);
        } else {
            // defer to  next round
            P = COMP(Aq1, P);
        }
    }

    // discriminants
    Ap = A1;
    while (Ap != NIL) {
        ADV(Ap,&Ap1,&Ap);

        if (PDEG(Ap1) < 2) continue;

        D = IPDSCRQE(r,Ap1);
        if (!IPCONST(r-1, D)) {
            P = COMP(D,P);
        }
    }

    // resultants
    Ap = A1;
    while (Ap != NIL) {
        ADV(Ap,&Ap1,&Ap);

        App = Ap;
        while (App != NIL) {
            ADV(App,&Ap2,&App);

            R = IPRESQE(r,Ap1,Ap2);
            if (!IPCONST(r-1, P)) {
                P = COMP(R,P);
            }
        }
    }

    return P;
}

// use CAD projection to find x-values of a set (not system) of equations
Word ProjSolve(Word r, Word A)
{
    Word J = A;
    while (r > 1) {
        J = ProjMcx(r, J);
        --r;
    }

    return J;
}

// find the local maxima and minima of function f in Q[x_1,...,x_r] subject to constraints
// Gs in Q[k_1,...,x_i], 2 <= i < r using the method of Lagrange Multipliers.
// return a list of polynomials, each h in Q[x_1]
// such that the roots of h give the x_1-coordinates of the critical points of f.
Word LagrangeRefinement(Word r, Word f, Word i, Word Gs, Word Is)
{
    Word Q = JACOBI(r, f, i, Gs, Is);

    // check for zero polynomial. no solutions
    if (Q == 0) return NIL;

    SWRITE("\njacobi det ");
    LWRITE(Q); SWRITE("\n");

    // find solution in x by projecton
    Word J = ProjSolve(r, COMP(Q, Gs));
    SWRITE("proj solve "); LWRITE(J); SWRITE("\n");

    return J;
}

// iterative application of lagrange multipliers then projection to x1.
// r  : number of variables
// Gs : minimal constraints
// P  : polynomial denoting top / bottom of sector cell
// Fs : polynomials constituting map from sector cell to section cell
// lagrange is done on P and each element of Fs, adding more constraints each time.
Word Refinement(Word r, Word Gs, Word P, Word Fs)
{
    Word Rs, Q, i, k, Is;
    // generate sequence Is = (1,...,k-1)
    i = 1;
    k = LENGTH(Gs) + 1;
    Is = NIL;
    while (i < k) {
        Is = COMP(i, Is);
        ++i;
    }

    Rs = NIL; // refinement polynomials

    // semi-monotone: crptical points of P subject to Gs
    if (Gs != NIL) {
        Rs = CONC(LagrangeRefinement(r, P, k, Gs, Is), Rs);
    }

    // monotone
    while (Fs != NIL) {
        ADV(Fs, &Q, &Fs);
        Gs = COMP(P, Gs);
        Is = COMP(k, Is);
        P = Q;
        ++k;

        // monotone, index k: critical points of Q subject to Gs
        Rs = CONC(LagrangeRefinement(r, Q, k, Gs, Is), Rs);
    }

    // solve for x1
    return Rs;
}

Word PrepHelper(Word k, Word P)
{
    // base case: polynomial in one variable
    if (k == 1) return LIST2(0, P);

    // recursive case
    return LIST2(0, PrepHelper(k - 1, P));
}

// polynomial prepend variables
// k: number of variables to prepend
// P: in I[x_1,...,x_r]
// return: P written in I[y_1,...,y_k,x_1,...,x_r]
Word PPREPVS(Word k, Word P)
{
    if (k == 0) return P;

    Word P1 = NIL;
    Word e, Q;
    while (P != NIL) {
        ADV2(P, &e, &Q, &P);

        P1 = COMP2(PrepHelper(k, Q), e, P1);
    }

    return INV(P1);
}

void QepcadCls::MONOTONE(Word A, Word J, Word D, Word r)
{
    // consider each true cell in D
    Word TrueCells, junk;
    LISTOFCWTV(D, &TrueCells, &junk);
    Word Ds = LELTI(D, CHILD); // cells of D, for searching cells
    Word AA = INV(LCOPY(A)); // in reverse order, to match SIGNPF. for later.

    while (TrueCells != NIL) {
        Word C;
        ADV(TrueCells, &C, &TrueCells);
        Word I = LELTI(C, INDX);

        // C has dimension two, then IJ < Ik are the positions in I where the composent is equal to 1. otherwise skip
        Word Ij, Ik;
        if (!TwoDimIndex(I, &Ij, &Ik)) {
            printf("skip cell\n");
            continue;
        }

        Word Ij1 = Ij - 1;
        Word nv = r - Ij1;

        // TODO debugging
        SWRITE("----------\n");
        LWRITE(I);
        printf(" 2d indx: (%d, %d)\n", Ij, Ik);

        // C0 := proj_{j-1}(C) is a 0-dimensional cell (c_1,...,c_{j-1})
        Word C0 = D;
        if (Ij > 1) {
            C0 = FindByIndex(Ds, I, Ij1, 1);
        }

        // sample point (c_1,...,c_{j-1}) of C0.
        Word S0 = LELTI(C0, SAMPLE);

        // top and bottom of proj_k(C) are one-dimensional sections by definition.
        Word I1 = LCOPY(I);
        // top: (i_1,...,i_k + 1)
        SLELTI(I1, Ik, LELTI(I, Ik) + 1);
        Word CT = FindByIndex(Ds, I1, Ik, 1);
        // bottom: (i_1,...,i_k - 1)
        SLELTI(I1, Ik, LELTI(I, Ik) - 1);
        Word CB = FindByIndex(Ds, I1, Ik, 1);

        // TODO debugging
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
        Word Gs = ZeroPolsSub(AA, r, C, Ij + 1, Ik - 1, S0, Ij1, nv);
        Word LL = Gs, P, Q; // TODO debug
        SWRITE("Gs:\n");
        while (LL != NIL) {
            ADV(LL, &P, &LL);
            IPDWRITE(nv, P, GVVL); SWRITE("\n");
        }

        // Fk (f_{k,T}, f_{k,B}) are 0 on CT and CB respectively.
        Word FT, FB;
        if (CT != NIL) FT = FIRST(ZeroPolsSub(AA, r, CT, Ik, Ik, S0, Ij1, nv));
        if (CB != NIL) FB = FIRST(ZeroPolsSub(AA, r, CB, Ik, Ik, S0, Ij1, nv));

        if (CT != NIL) { SWRITE("f_top := "); IPDWRITE(nv, FT, GVVL); SWRITE("\n"); }
        if (CB != NIL) { SWRITE("f_bottom := "); IPDWRITE(nv, FB, GVVL); SWRITE("\n"); }

        // Fs = (f_{K+1},...,f_n) is a map from proj_{k}(C) to R^{n-k}, of which C is the graph.
        Word Fs = ZeroPolsSub(AA, r, C, Ik + 1, r, S0, Ij1, nv);
        SWRITE("Fs:\n");
        LL = Fs;
        while (LL != NIL) {
            ADV(LL, &P, &LL);
            IPDWRITE(nv, P, GVVL); SWRITE("\n");
        }

        Word Rs = NIL; // list of refinement polynomials, in I[x_Ij]
        // perform refinement
        if (CT != NIL) {
            SWRITE("top\n");
            Rs = CONC(Refinement(nv, Gs, FT, Fs), Rs);
        }

        if (CB != NIL) {
            SWRITE("bottom\n");
            Rs = CONC(Refinement(nv, Gs, FB, Fs), Rs);
        }

        // add refinement polynomials
        Word Rs1 = NIL;
        while (Rs != NIL) {
            Word P;
            ADV(Rs, &P, &Rs);
            LWRITE(P); SWRITE("\n");

            // write in Q[x_1,...,x_Ij]
            Rs1 = COMP(MPOLY(PPREPVS(Ij1, P), NIL, NIL, PO_POLY, PO_KEEP), Rs1);

            LWRITE(PPREPVS(Ij1, P)); SWRITE("\n");
        }

        // add refinement polynomials to J and their factors to A.
        Rs = IPLFAC(Ij, Rs1);
        ADDPOLS(Rs1, Ij, LFS("K"), &J);
        ADDPOLS(Rs, Ij, LFS("M"), &A);
    }
}

