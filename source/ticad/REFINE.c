/*======================================================================
                      D <- REFINE(k,D,A,P)

Refine CAD D so that its zero-cells C0 are compatible with relevant new PO_REFINEMENT polynomials in A

\Input
  \parm{k} is the level of D
  \parm{D} is the original cad
  \parm{A} list of refinement polynomials (same structure as projection factor set)
  \parm{P} list of projection factors (list of lists)

Output
  \parm{D} is the refinement of D compatible with new P_MONOTONE polynamials in A

======================================================================*/
#include "qepcad.h"

// compare two algebraic numbers
// if a or b is rational, then this function expects M = PMON(1,1), I = (r,r) where r is the rational number
// if both algebraic and not equal, rational isolating interval refinements until separated are used.
// TODO how many interval refinements needed? can we have infinite refinement?
Word COMPARE(Word M1, Word *I1_, Word M2, Word *I2_)
{
    bool a_rat = PDEG(M1) == 1;
    bool b_rat = PDEG(M2) == 1;

    // easy case: rational comparison
    if (a_rat && b_rat) {
        return RNCOMP(FIRST(*I1_), FIRST(*I2_));
    }

    // otherwise, at least one of a and b is algebraic. compare using rational approximation.
    Word a1, a2, b1, b2;
    while(1) {
        FIRST2(*I1_, &a1, &a2);
        FIRST2(*I2_, &b1, &b2);

        // check whether the intervals are separated, if not, then refine them until they are.
        Word c1 = RNCOMP(a2, b1);

        if (c1 < 0) {
            return -1;
        }

        Word c2 = RNCOMP(b2, a1);
        if (c2 < 0) {
            return 1;
        }

        // now we know that the intervals overlap
        // two possibilities:
        // maybe algebraic and equal?
        if (c1 == 0 && c2 == 0 && EQUAL(M1, M2)) {
            return 0;
        }

        // or we can refine the isolating intervals and do endpoint comparisons again.
        if (!a_rat) *I1_ = IUPIIR(M1, LIST2(a1,a2));
        if (!b_rat) *I2_ = IUPIIR(M2, LIST2(b1,b2));
    }

    return 0;
}

// set cell C index, element k to value a
// also updates all children
void SETINDEXK(Word C, Word k, Word a)
{
    Word I = LELTI(C, INDX);
    SLELTI(I, k, a);

    Word Children = LELTI(C, CHILD);
    while (Children != NIL) {
        Word C1;
        ADV(Children, &C1, &Children);

        SETINDEXK(C1, k, a);
    }
}

// are the first n-1 coordinates of S rational?
bool SampleIsRat(Word SM, Word Sb)
{
    // base polynomial is rational
    if (PDEG(SM) == 1) {
        return true;
    }

    // otherwise, the point is in primitive form. we must check whether all c_1,...,c_{k-1} are all rational.
    while (Sb != NIL) {
        Word b;
        ADV(Sb, &b, &Sb);

        // algebraic? (non-zero and nonconstant polynomial for alpha)
        Word junk;
        if (!AfIsRat(b, &junk)) {
            return false;
        }
    }

    // otherwise first k-1 coordinates are ratinoal
    return true;
}

// helper function: convert sample to primitive, no timekeeping
void ConvertToPrimitive(Word SQ, Word SJ, Word SM, Word SI, Word Sb, Word* M_, Word* I_, Word* b_)
{
    Word G, K, t, u, a, b, junk;
    SIMPLEQE(SM,SI,SQ,SJ,&G,&t,&u,&K,&a,&b,&junk,&junk);

    if (u != 0) {
        MODCRDB(Sb,G,a,b,&Sb);
    } else {
        Sb = CCONC(Sb,LIST1(b));
    }

    *M_ = G;
    *I_ = K;
    *b_ = Sb;
}

// set sample point, recursive helper function
// S is *always a primitive* sample point for a cell C
// Ch is the list of children of C
// M, I is an algebraic or rational value to append to the sample point.
Word SetSampleHelper(Word r, Word S, Word Ch, Word M, Word I, Word PFs)
{
    Word SM, SI, Sb;
    FIRST3(S, &SM, &SI, &Sb);

    // is the new sample rational?
    Word S1;
    if (PDEG(M) == 1) {
        // then append it
        Word c = AFFRN(FIRST(I));
        Word Sb1 = CCONC(Sb, LIST1(c));

        S1 = LIST3(SM, SI, Sb1);
    } else if (SampleIsRat(SM, Sb)) {
        // k-th coordinate is algebraic, but we can use primitive representation.
        Word Sb1 = CCONC(Sb, LIST1(AFGEN()));

        S1 = LIST3(M, I, Sb1);
    } else if (Ch == NIL) {
        // algebraic and no children to update -- extended
        S1 = LIST5(AFPFIP(1,M), I, SM, SI, Sb);
    } else {
        // algebraic but has children to update -- convert to primitive
        Word SM1, SI1, Sb1;
        ConvertToPrimitive(AFPFIP(1,M), I, SM, SI, Sb, &SM1, &SI1, &Sb1);

        S1 = LIST3(SM1, SI1, Sb1);
    }

    // the sample points of all children are now wrong. update those recursively.
    if (Ch == NIL) {
        // if there are 0 children there is nothing to do
        return S1;
    } else if (LENGTH(Ch) == 1) {
        // if there is one child, sample point can be any number. Let it be 0.
        Word C = FIRST(Ch);
        SetSampleHelper(r+1, S1, LELTI(C, CHILD), PMON(1,1), LIST1(RNINT(0)), RED(PFs));

        return S1;
    }

    // recompute sample points using root finding
    // there are k > 1 child cells.
    // since polynomials are delineable before refinement, they are still delineable after refinement,
    // therefore FIRST(PFs) is a system with k roots
    Word Ps;
    ADV(PFs, &Ps, &PFs);

    // evaluate the polynomials at new sample point.
    Word SPs = NIL;
    while (Ps != NIL) {
        Word P;
        ADV(Ps, &P, &Ps);

        Word P1 = SUBSTITUTE(r, LELTI(P, PO_POLY), S1, true); // with rational coefficients
        SPs = COMP(P1, SPs);
    }

    // and find their roots
    Word B = ROOTS(SPs, LIST2(NIL, NIL));
    while (B != NIL && Ch != NIL) {
        Word C1, C2, RM, RI;
        ADV2(Ch, &C1, &C2, &Ch);
        ADV2(B, &RI, &RM, &B);

        // TODO if the first root is not algebraic, then the isolating interval is not sufficient. same is true for the
        // last root. maybe we should do C_B,C instead, and set the first sector to floo or first interval - 1.
        Word SC1 = SetSampleHelper(r+1, S1, LELTI(C1, CHILD), PMON(1,1), LIST1(FIRST(RI)), PFs);

        // if root is rational
        Word SC2;
        if (PDEG(RM) == 1) {
            Word c = IUPRLP(RM);

            SC2 = SetSampleHelper(r+1, S1, LELTI(C2, CHILD), PMON(1,1), LIST1(c), PFs);
        } else { // polynomial and isolating interval define the algebraic number
            SC2 = SetSampleHelper(r+1, S1, LELTI(C2, CHILD), RM, RI, PFs);
        }

        SLELTI(C1, SAMPLE, SC1);
        SLELTI(C2, SAMPLE, SC2);
    }

    return S1;
}

// let C be a level k cell. set index K of sample point to the algebraic number M,I
// if rational, algebraic number should be in the same format as for COMPARE
void SETSAMPLE(Word C, Word M, Word I, Word PFs)
{
    // get the sample point, put it in the correct format for the helper function
    Word k = LELTI(C, LEVEL);
    Word S = LELTI(C, SAMPLE);
    Word Ch = LELTI(C, CHILD);

    // since we will set the k-th coordinate, the last coordinate should be discarded.
    Word junk, SM, SI, Sb;
    if (LENGTH(S) == 5) { // extended (we don't care about the last coordinate)
                          // if in extended form, the first k-1 coordinates are given.
        FIRST5(S, &junk, &junk, &SM, &SI, &Sb);
    } else {
        // otherwise, throw out the last coordinate.
        FIRST3(S, &SM, &SI, &Sb);
        if (LENGTH(Sb) == 1) {
            SM = PMON(1,1);
            SI = LIST2(0,0);
            Sb = NIL;
        } else {
            Sb = INV(RED(INV(Sb))); // clumsy way to delete the last element
        }
    }

    Word S1 = SetSampleHelper(k+1, LIST3(SM, SI, Sb), Ch, M, I, PFs);
    SLELTI(C, SAMPLE, S1);
}

// let C = FIRST(Cs) be a (0,...,0,1)-cell and s be a point in C. refine C into three new cells such that s is a new
// (0,...,0,0)-cell
Word RefineCell(Word k, Word Cs, Word PM, Word PI, Word c, Word PFs, bool* rc)
{
    Word Cs2 = Cs;

    // Let C = (a,b). C becomes (a,s), C2 becomes new cell s and C3 new cell (s,b)
    Word C1;
    ADV(Cs, &C1, &Cs);
    Word C2 = LDCOPY(C1);
    Word C3 = LDCOPY(C1);

    // k-th element of C index.
    Word j = LELTI(LELTI(C1, INDX), k);

    SWRITE("Refine cell "); LWRITE(LELTI(C1, INDX)); SWRITE("\n");

    // update indices
    SETINDEXK(C2, k, ++j);
    SETINDEXK(C3, k, ++j);

    // update sample
    // we will need to update only two of the cells, as the existing sample will be correct for one of them
    Word SQ, SJ;
    GETSAMPLEK(-1, LELTI(C1, SAMPLE), &SQ, &SJ);

    Word sign = COMPARE(SQ, &SJ, PM, &PI);
    // -1: C1 is correct, 0: C2 is correct, +1: C3 is correct.

    if (sign != -1) { // need to update C1 ...
        // to the given rational point c
        SETSAMPLE(C1, PMON(1,1), LIST1(c), PFs);
    }

    if (sign != 0) { // need to update C2 ...
        // to the new "refinement point" given
        SETSAMPLE(C2, PM, PI, PFs);
    }

    // we may will need to update the sample of C3, but this is done later.
    *rc = sign != 1;

    // increment indices of remaining cells.
    Word Cs1 = Cs;
    while (Cs1 != NIL) {
        Word C;
        ADV(Cs1, &C, &Cs1);

        SETINDEXK(C, k, ++j);
    }

    // append new cells
    SRED(Cs2, COMP2(C2, C3, Cs));
    return Cs2;
}

// let P in Q(alpha)[x] such that every coefficient is rational. Convert P into a polynomial in Z[x]
// requires less computation than AFPNORM
Word STRIP(Word P)
{
    Word P1, e, A, a, b;
    P1 = NIL;
    while (P != NIL) {
        ADV2(P, &e, &A, &P);
        FIRST2(A, &a, &b);
        if (PDEG(b) > 0) {
            perror("STRIP: polynomial has algebraic coefficients.\n");
            return 0;
        }

        P1 = COMP2(a, e, P1);
    }

    P1 = INV(P1);
    return IPFRPmod(1, P1);
}

// convenience function for next polynomial
void NextPolynomial(Word Ps, Word* PM_, Word* PI_, Word* J_, Word* Ps_)
{
    if (Ps == NIL) {
        *PM_ = NIL;
        *PI_ = NIL;
        *J_ = NIL;

        return;
    }

    Word P, M, I, b, a;
    ADV(Ps, &P, Ps_);
    Word P1 = FIRST(P);

    // extended form, can assume the polynomial is normalised
    if (LENGTH(P1) == 5) {
        FIRST2(P1, &M, &I);
        M = STRIP(M);
    } else {
        FIRST3(P1, &M, &I, &b);
        b = LAST(b);
        if (AfIsRat(b, &a)) {
            M = PMON(1,1);
            I = LIST2(a,a);
        }

        // otherwise, the algebraic refinement point coordinate is always given by AFGEN(), so just use M and I.
    }

    *PM_ = M;
    *PI_ = I;
    *J_ = SECOND(P);
}

// Refine subcad D to be compatible with level 1 polynomials Ps
Word RefineSubcad(Word k, Word Ch, Word Ps, Word PFs)
{
    Word Ch2 = Ch; // backup of list pointer, to return.
    Word PM, PI, J;

    // Ps is a list of sample points in ascending order, which define new 0-cells we will add to D
    NextPolynomial(Ps, &PM, &PI, &J, &Ps);
            SWRITE("first polynomial ");
            LWRITE(PM); SWRITE(" ");
            LWRITE(PI); SWRITE("\n");


    Word sign = -1;
    Word Ch1 = Ch;
    bool refined = false;
    while (Ch != NIL) {
        Word C, CT, PM1, PI1;
        // next sector cell
        ADV(Ch, &C, &Ch1);

        // get next polynomial if needed.
        if (sign != -1) {
            PM1 = PM, PI1 = PI;
            NextPolynomial(Ps, &PM, &PI, &J, &Ps);
            SWRITE("next polynomial ");
            LWRITE(PM); SWRITE(" ");
            LWRITE(PI); SWRITE("\n");
        }

        // last cell TODO
        //         SETSAMPLE(C, PMON(1,1), RNSUM(RNINT(1), RNCEIL(SECOND(J))), PFs);

        Word SQ, SJ;
        // refine the last cell
        if (Ch1 != NIL) {
            // next section -- top of C
            ADV(Ch1, &CT, &Ch1);

            // get k-th coordinate of the sample point of C
            GETSAMPLEK(-1, LELTI(CT, SAMPLE), &SQ, &SJ);
            SWRITE("cell top ");
            LWRITE(SQ); SWRITE(" ");
            LWRITE(SJ); SWRITE("\n");

            if (PM != NIL) {
                sign = COMPARE(SQ, &SJ, PM, &PI);
                printf("sign = %d\n", sign);
            }
        } else if (PM != NIL) {
            SQ = PMON(1,1);

            Word c1 = RNSUM(SECOND(J), RNINT(2));
            SJ = LIST2(c1, c1);

            sign = 1; // force a refinement.
        }

        // previous cell was refined, need to update
        if (refined) {
            Word c = RNQ(RNSUM(SECOND(PI1), FIRST(SJ)), RNINT(2));
            SWRITE("refine "); RNWRITE(c); SWRITE("\n");
            SETSAMPLE(C, PMON(1,1), LIST1(c), PFs);
            refined = false;
        }

        // no more polynomials
        if (PM == NIL) break;

        if (sign > 0) { // S > P, refine previous sector, FIRST(Ch)
            Word c = PDEG(PM) == 1 ? FIRST(J) : FIRST(PI);
            Ch1 = RED2(RefineCell(k, Ch, PM, PI, c, PFs, &refined));
        }

        Ch = Ch1;
    }

    return Ch2;
}

// are the first k elements in L1 and L2 equal
bool EQUALK(Word k, Word L1, Word L2)
{
    Word i, a, b;

    i = 0;
    while (L1 != NIL && L2 != NIL && i < k) {
        ADV(L1, &a, &L1);
        ADV(L2, &b, &L2);
        ++i;

        if (a != b) return false;
    }

    // all elements checked are equal, also need to check lists have the same length
    return i == k;
}

Word QepcadCls::REFINE(Word k, Word D, Word A, Word PF)
{
    // no children to refine.
    Word Ch = LELTI(D, CHILD);
    if (Ch == NIL) {
        return D;
    }

    PF = RED(PF);
    Word k1 = k-1;
    Word A1;
    ADV(A, &A1, &A); // decomstruct A. A1 is the set of level k+1 polys

    // find the new PO_REFINE polynomials.
    Word Ps = NIL, I = LELTI(D, INDX);
    while (A1 != NIL) {
        Word P;
        ADV(A1, &P, &A1);

        if (EQUALK(k1, LELTI(P, PO_REFINEMENT), I)) { // list equality check not needed, since same cell index pointer is used
            Ps = COMP(LELTI(P, PO_POLY), Ps);
        }
    }

    // do refinement if the list of Ps is non-empty
    if (Ps != NIL) {
        printf("refining subcad of "); LWRITE(LELTI(D, INDX)); SWRITE("\n");
        Ch = RefineSubcad(k, Ch, Ps, PF);
        SLELTI(D, CHILD, Ch);
    }

    // no more refinement polynomials
    if (A == NIL) {
        return D;
    }

    // walk the CAD, sections only.
    Word C, junk;
    ADV(Ch, &junk, &Ch);
    while (Ch != NIL) {
        ADV2(Ch, &C, &junk, &Ch);

        C = REFINE(k+1, C, A, PF);
    }

    return D;
}

