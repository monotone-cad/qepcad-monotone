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
#include "db/CAPolicy.h"
using namespace std;

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
                Q = PADDVS(SUBSTITUTE(i2, LELTI(P, PO_POLY), S0, false), l + n - i2); // with integer coefficients

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

    // factorise
    Word junk, Qs, Q1;
    IPFACDB(r ,Q, &junk, &junk, &Qs);
    while (Qs != NIL) {
        ADV(Qs, &Q1, &Qs);

        Gs = COMP(SECOND(Q1), Gs);
    }

    // simplify Gs by constructing a Groebner basis, if we have Singular
    if (GVCAP->supports("GROEBNER")) {
        Gs = GVCAP->GROEBNER(Gs, NIL, r);
    }

    // find solution in x by projecton
    // Gs now contains factors of the jacobi determinant
    return ProjSolve(r, Gs);
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
    k = LENGTH(Gs) + i;
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

// adds roots of refinement polynomials, as sample points, to a new set. These will be sample points of new 0-cells.
// TODO at the momont we add all roots, but we may only need roots inside specific cells. Fix this.
Word REFINEMENTPOINTS(Word r, Word Rs, Word* A_, Word* J_)
{
    // initialise refinement polynomials set
    Word RPs = NIL;
    for (int i = 0; i < r; ++i) RPs = COMP(NIL, RPs);

    while (Rs != NIL) {
        Word I, S, R1s;
        ADV3(Rs, &I, &S, &R1s, &Rs);

        ADDREFINEMENTPOINTS(I, S, R1s, A_, J_, &RPs);
    }

    return RPs;
}

// Store refinement polynomials Rs with base index I
// TODO also save max/min root of polynomial so we only add the relevant ponts rather than roots not outside the cell
void STOREPOLYNOMIALS(Word Rs, Word I, Word S, Word* A_)
{
    // find cell if exists
    Word A = *A_, Rs1 = NIL;
    while (A != NIL) {
        Word I1, S1, Ps;
        ADV3(A, &I1, &S1, &Ps, &A);

        // found
        if (I == I1) { // avoid proper list comparison because subcads are identical
            Ps = CONC(Ps, Rs);

            return;
        }
    }

    // not found
    *A_ = COMP3(I, S, Rs, *A_);
}

Word QepcadCls::MONOTONE(Word* A_, Word* J_, Word D, Word r)
{
    // to store unfactorised refinement polynomials, in the form (Index, Sample, Ps)
    Word Rs = NIL;

    // consider each true cell in D
    Word TrueCells, junk;
    LISTOFCWTV(D, &TrueCells, &junk);
    Word Ds = LELTI(D, CHILD); // cells of D, for searching cells
    Word AA = INV(LCOPY(*A_)); // in reverse order, to match SIGNPF. for later.

    while (TrueCells != NIL) {
        Word C;
        ADV(TrueCells, &C, &TrueCells);
        Word I = LELTI(C, INDX);

        // C has dimension two, then IJ < Ik are the positions in I where the composent is equal to 1. otherwise skip
        Word Ij, Ik;
        if (!TwoDimIndex(I, &Ij, &Ik)) {
            continue;
        }

        Word Ij1 = Ij - 1;
        Word nv = r - Ij1;

        // C0 := proj_{j-1}(C) is a 0-dimensional cell (c_1,...,c_{j-1})
        Word C0, S0, I0;
        if (Ij == 1) { // base CAD is D
            C0 = D;
            S0 = LIST3(PMON(1,1), LIST2(0,0), NIL);
            I0 = NIL;
        } else { // sub-cad on top of some 0-cell
            C0 = FindByIndex(Ds, I, Ij1, 1);
            S0 = LELTI(C0, SAMPLE);
            I0 = LELTI(C0, INDX);
        }

        // top and bottom of proj_k(C) are one-dimensional sections by definition.
        Word I1 = LCOPY(I);
        // top: (i_1,...,i_k + 1)
        SLELTI(I1, Ik, LELTI(I, Ik) + 1);
        Word CT = FindByIndex(Ds, I1, Ik, 1);
        // bottom: (i_1,...,i_k - 1)
        SLELTI(I1, Ik, LELTI(I, Ik) - 1);
        Word CB = FindByIndex(Ds, I1, Ik, 1);

        // find polynomials in sub-cad
        // TODO
        //   it takes a lot of work to find these polynomials, and work is often repeated.
        //   each label is unique, so we can check whether we have already substituted to save on effort.
        //   what other things could i cache? can i use the existing db faciility?
        // (note they will be in Z[x_i,...,x_l] after substitution):
        // Gs = g_2,...,g_{k-1} define proj_{k-1}(C).
        Word Gs = ZeroPolsSub(AA, r, C, Ij + 1, Ik - 1, S0, Ij1, nv);

        // Fk (f_{k,T}, f_{k,B}) are 0 on CT and CB respectively.
        Word FT, FB;
        if (CT != NIL) FT = FIRST(ZeroPolsSub(AA, r, CT, Ik, Ik, S0, Ij1, nv));
        if (CB != NIL) FB = FIRST(ZeroPolsSub(AA, r, CB, Ik, Ik, S0, Ij1, nv));

        // Fs = (f_{K+1},...,f_n) is a map from proj_{k}(C) to R^{n-k}, of which C is the graph.
        Word Fs = ZeroPolsSub(AA, r, C, Ik + 1, r, S0, Ij1, nv);

        // perform refinement
        if (CT != NIL) {
            STOREPOLYNOMIALS(Refinement(nv, Gs, FT, Fs), I0, S0, &Rs);
        }

        if (CB != NIL) {
            STOREPOLYNOMIALS(Refinement(nv, Gs, FB, Fs), I0, S0, &Rs);
        }
    }

    return REFINEMENTPOINTS(r, Rs, A_, J_);
}

