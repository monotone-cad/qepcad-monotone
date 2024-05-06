/*======================================================================
                   R -< MONOTONE(A, J, D, r)

Refinements of the kind { x_i < c }, { x_i = c }, { x_i > c } above a 0-dimensional cell in R^{i-1} will be performed
such that every cell is monotone. Monotone functions are obtained using Lagrange Multipliers.

\Input
  \parm{A} projection factor set
  \parm{J} projection polynomial set
  \parm{D} base CAD (of R^0)
  \parm{r} positive integer, dimension of ambient space

\Output
  \parm{R} set of *refinement points*: a set R_1, ..., R_r
           such that each R_i is a set of i-dimensional sample points c_j which define the refinement { x_i < c_j }, {
           x_i = c_j }, { x_i > c_j }.

Side Effect
  refinement polynomials (roots of which are the refinement points) re added to J

======================================================================*/
#include "qepcad.h"
#include "db/CAPolicy.h"
using namespace std;

// suppose I = (i_1,...,i_n) is a cell index.
// if i_j = i_k = 1 but all other components are 0, then C is 2-dimensional.
// set j and k to indicate which components are equal to one and return 2.
// otherwise, C is not two-dimensional. return d := dim(C) if dim(C) < 2, otherwise d := 3. Values of j and k are undefined.
int TwoDimIndex(Word I, Word *j_, Word *k_)
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

    // if cell is not 2-dimensional, return false
    if (d != 2) {
        return d;
    }

    // otherwise, cell is two-dimensional
    FIRST2(J, k_, j_);

    return 2;
}

// return a list (f_i1,...,f_i2), where each f_i is a level i polynomial which is 0 on C with sample point S substituted in.
// A : (reversed) projection factor set
// r : positive integer, dimension of CAD
// C : cell to check
// i1 <= i2 : min and max indices to search
// S0 : sample point of level l < i1 to substitute. give the empty sample point of D to do no substitution.
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

    return L; // note: L is in ascending level order, since A is reversed.
}

// use CAD projection to solve a set of multivariate polynomials for x_1
// A : set of polynomials in r variables
// r : positive integer
// return : set of univariate polynomials whose roots are x_1 coordinates of zeros of A
Word ProjSolve(Word r, Word A)
{
    Word J = A;
    while (r > 1) {
        J = ProjMcxUtil(r, J);
        --r;
    }

    return J;
}

// find the local maxima and minima of function f in Q[x_1,...,x_r] subject to constraints
// Gs in Q[k_1,...,x_i], 2 <= i < r using the method of Lagrange Multipliers.
// return a list of polynomials, each h in Q[x_1]
// such that the roots of h give the x_1-coordinates of the critical points of f subject to Gs.
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

    // simplify Gs by constructing a Groebner basis is supported
    // this step makes solving the jacobi determinant for x_1 a lot quicker in practice.
    if (GVCAP->supports("GROEBNER")) {
        Gs = GVCAP->GROEBNER(Gs, NIL, r);
    }

    // find solution in x_1 by projecton
    // Gs now includes factors of the jacobi determinant
    return ProjSolve(r, Gs);
}

// iterative application of lagrange multipliers then projection to x_1.
// r  : number of variables
// Gs : minimal constraints
// P  : polynomial denoting top / bottom of sector cell
// Fs : polynomials constituting map from sector cell to section cell
// lagrange multipliers is done first with P subject to Gs, then on each subsequent round, a polynomial from F is used
// and the optimisation polynomial used in the previous round is appended to the constraints.
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

    // semi-monotone: critical points of P subject to Gs
    if (Gs != NIL) {
        Rs = CONC(LagrangeRefinement(r, P, k, Gs, Is), Rs);
    }

    // monotone: ecah F_i subject to Gs, F_1,...F_{i-1} for each i
    while (Fs != NIL) {
        ADV(Fs, &Q, &Fs);
        Gs = COMP(P, Gs);
        Is = COMP(k, Is);
        P = Q;
        ++k;

        // monotone, index k: critical points of Q subject to Gs
        Rs = CONC(LagrangeRefinement(r, Q, k, Gs, Is), Rs);
    }

    // solve for x_1
    return Rs;
}

// constructs a set of *refinement points* -- roots of refinement polynomials as sample points. Each of these sample
// points will be the sample point of a new 0-cell in the refinement of the CAD.
// TODO at the momont we add all roots, but we may only need roots inside specific cells. Fix this.
// TODO fix to work with total derivatives. i.e., discard roots at which partial f / partial x_r is zero
Word REFINEMENTPOINTS(Word r, Word Rs, Word* A_, Word* J_)
{
    // initialise refinement polynomials set
    Word RPs = NIL;
    for (int i = 0; i < r; ++i) RPs = COMP(NIL, RPs);

    // add refinement points using refinement polnomials
    while (Rs != NIL) {
        Word I, S, R1s, Endpoints;
        ADV4(Rs, &I, &S, &R1s, &Endpoints, &Rs);

        ADDREFINEMENTPOINTS(I, S, R1s, Endpoints, A_, J_, &RPs);
    }

    return RPs;
}

// Store refinement polynomials Rs with base index I
// TODO only save relevant roots, i.e., those inside some interval (1-cell)
void STOREPOLYNOMIALS(Word Rs, Word I, Word S, Word Endpoints, Word* A_)
{
    // find cell index if exists
    Word A = *A_, Rs1 = NIL;
    while (A != NIL) {
        Word I1, S1, E1, Ps;
        ADV4(A, &I1, &S1, &E1, &Ps, &A);

        if (I == I1) { // avoid proper list comparison because subcads are identical
            Ps = CONC(Ps, Rs);

            return;
        }
    }

    // cell index not yet included, make new
    *A_ = COMP4(I, S, Rs, Endpoints, *A_);
}

Word QepcadCls::MONOTONE(Word* A_, Word* J_, Word D, Word r)
{
    // to store unfactorised refinement polynomials, in the form (Index, Sample, Ps)
    Word Rs = NIL;

    // consider each true cell in D
    // note: in practice, it is slower and takes more memory to do a walk of the CAD, saving polynomials along the way.
    Word TrueCells, junk;
    LISTOFCWTV(D, &TrueCells, &junk);
    Word Ds = LELTI(D, CHILD); // cells of D, for searching cells
    Word AA = INV(LCOPY(*A_)); // in reverse order, to match SIGNPF. for later.

    while (TrueCells != NIL) {
        Word C;
        ADV(TrueCells, &C, &TrueCells);

        // C has dimension two, then IJ < Ik are the positions in I where the composent is equal to 1. otherwise skip
        Word I = LELTI(C, INDX);
        Word Ij, Ik;
        Word d = TwoDimIndex(I, &Ij, &Ik);
        if (d > 2) {
            // fail: only implemented for dimension at most 2
            FAIL("source/ticad/MONOTONE", "cell dimension greater than 2 not supported.");
        } else if (d < 2) {
            // cells with dimension 0 and 1 are already monotone.
            continue;
        }

        Word Ij1 = Ij - 1;
        Word nv = r - Ij1;

        // C0 := proj_{j-1}(C) is a 0-dimensional cell (c_1,...,c_{j-1})
        Word C0, S0, I0;
        Word ij = LELTI(I, Ij);
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

        // Top and bottom of 1-sector, level J
        Word C0C = LELTI(C0, CHILD);
        Word I01 = CCONC(I0, LIST1(ij));
        Word S0B = NIL, S0T = NIL, Endpoints;
        if (ij > 1) {
            GETSAMPLEK(
                Ij,
                LELTI(LELTI(C0C, ij - 1), SAMPLE),
                &junk,
                &S0B
            );

            S0B = SECOND(S0B);
        }
        if (ij < LENGTH(C0C)) {
            GETSAMPLEK(
                Ij,
                LELTI(LELTI(C0C, ij + 1), SAMPLE),
                &junk,
                &S0T
            );

            S0T = SECOND(S0T);
        }
        Endpoints = LIST2(S0B, S0T);

        // find polynomials in sub-cad
        // (note they will be in Z[x_i,...,x_l] after substitution):
        // Gs = g_2,...,g_{k-1} define proj_{k-1}(C).
        // TODO can we save the polynomials, maybe using cell position, to save duplicating work?
        Word Gs = ZeroPolsSub(AA, r, C, Ij + 1, Ik - 1, S0, Ij1, nv);

        // Fk (f_{k,T}, f_{k,B}) are 0 on CT and CB respectively.
        Word FT, FB;
        if (CT != NIL) FT = FIRST(ZeroPolsSub(AA, r, CT, Ik, Ik, S0, Ij1, nv));
        if (CB != NIL) FB = FIRST(ZeroPolsSub(AA, r, CB, Ik, Ik, S0, Ij1, nv));

        // Fs = (f_{K+1},...,f_n) is a map from proj_{k}(C) to R^{n-k}, of which C is the graph.
        Word Fs = ZeroPolsSub(AA, r, C, Ik + 1, r, S0, Ij1, nv);

        // perform refinement
        if (CT != NIL) {
            STOREPOLYNOMIALS(Refinement(nv, Gs, FT, Fs), I01, S0, Endpoints, &Rs);
        }

        if (CB != NIL) {
            STOREPOLYNOMIALS(Refinement(nv, Gs, FB, Fs), I01, S0, Endpoints, &Rs);
        }
    }

    return REFINEMENTPOINTS(r, Rs, A_, J_);
}

