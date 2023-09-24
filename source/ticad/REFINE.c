/*======================================================================
                      D <- REFINE(k,D,A,P)

Refine CAD D so that its zero-cells C0 are compatible with relevant new PO_REFINEMENT polynomials in A

\Input
  \parm{k} is the level of D
  \parm{D} is the original cad
  \parm{A} list of refinement polynomials (list of lists)
  \parm{P} list of projection factors (list of lists)

Output
  \parm{D} is the refinement of D compatible with new P_MONOTONE polynamials in A

======================================================================*/
#include "qepcad.h"

// need a new property of cells, SIGNRP, like SIGNPF. this wiss contain signs of refinement polynomials, level k, above
// the sub-cad as indicated in the polynomial. add in cell write.

// is algebraic field element a rational number? return true if rational and store the rational number in r_
// otherwise, return false, if algebraic then r_ is not modified
inline bool AfIsRat(Word b, Word* r_)
{
    if (b == 0) {
        *r_ = 0;

        return true;
    }

    Word p, Q;
    FIRST2(b, &p, &Q);
    *r_ = p;

    return PDEG(Q) == 0;
}

// compare two algebraic numbers
// if a or b is rational, then this function expects M = PMON(1,1), I = (r,r) where r is the rational number
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
// child recursive update not needed.
void SETINDEXK(Word C, Word k, Word a)
{
    Word I = LELTI(C, INDX);
    SLELTI(I, k, a);
}

// are the first n-1 coordinates of S rational
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
Word SetSampleHelper(Word S, Word Ch, Word M, Word I, Word PFs)
{
    Word SM, SI, Sb;
    FIRST3(S, &SM, &SI, &Sb);

    // is the new sample rational?
    Word S1;
    if (PDEG(M) == 1) {
        // then append it
        Word c = AFFRN(FIRST(I));
        Sb = CONC(Sb, LIST1(c));

        S1 = LIST3(SM, SI, Sb);
    } else if (SampleIsRat(SM, Sb)) {
        // k-th coordinate is algebraic, but we can use primitive representation.
        Sb = CONC(Sb, LIST1(AFGEN()));

        S1 = LIST3(M, I, Sb);
    } else if (Ch == NIL) {
        // algebraic and no children to update -- extended
        S1 = LIST5(AFPFIP(1,M), I, SM, SI, Sb);
    } else {
        // algebraic but has children to update -- convert to primitive
        Word SM1, SI1, Sb1;
        ConvertToPrimitive(AFPFIP(1,M), I, SM, SI, Sb, &SM1, &SI1, &Sb1);

        S1 = LIST3(SM1, SI1, Sb1);
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
        Sb = INV(RED(INV(Sb))); // clumsy way to delete the last element
    }

    Word S1 = SetSampleHelper(S, Ch, M, I, PFs);
    SLELTI(C, SAMPLE, S1);
}

// let C = FIRST(Cs) be a (0,...,0,1)-cell and s be a point in C. refine C into three new cells such that s is a new
// (0,...,0,0)-cell
Word RefineCell(Word k, Word Cs, Word SQ, Word SJ, Word PM, Word PI, Word PFs, bool* rc)
{
    Word Cs2 = Cs;

    // Let C = (a,b). C becomes (a,s), C2 becomes new cell s and C3 new cell (s,b)
    Word C1;
    ADV(Cs, &C1, &Cs);
    Word C2 = LDCOPY(C1);
    Word C3 = LDCOPY(C1);

    // k-th element of C index.
    Word j = LELTI(LELTI(C1, INDX), k);

    printf("refine cell "); LWRITE(LELTI(C1, INDX));

    // update indices
    SETINDEXK(C2, k, ++j);
    SETINDEXK(C3, k, ++j);

    // update sample
    // we will need to update only noe cell, as the existing sample will be correct for one of them
    Word sign = COMPARE(SQ, &SJ, PM, &PI);
    // -1: C1 is correct, 0: C2 is correct, +1: C3 is correct.

    if (sign != -1) { // need to update C1 ...
                      // ... to the left-hand end of the isolating interval
        SETSAMPLE(C1, PMON(1,1), PI, PFs);
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

// convenience function: get k-th coordinate of the sample point as an algebraic number
void GETSAMPLEK(Word S, Word* Q_, Word* J_)
{
    Word SQ, SJ, SM, SI, Sb, junk;
    if (LENGTH(S) == 5) {
        FIRST4(S, &SQ, &SJ, &SM, &SI);
        SQ = AFPNORM(1, SM, SQ);
    } else {
        FIRST3(S, &SM, &SI, &Sb);
        Sb = LAST(Sb);
        Word b;
        if (AfIsRat(Sb, &b)) {
            SQ = PMON(1,1);
            SJ = LIST2(b, b);
        } else {
            ANFAF(SM, SI, Sb, &SQ, &SJ);
        }
    }

    *Q_ = SQ;
    *J_ = SJ;
}

// Refine subcad D to be compatible with level 1 polynomials Ps
Word RefineSubcad(Word k, Word Ch, Word Ps, Word PFs)
{
    Word Ch2 = Ch; // backup of list, to return.
    Word PM, PI, J;

    // Ps is a list of sample points in ascending order, which define new 0-cells we will add to D
    NextPolynomial(Ps, &PM, &PI, &J, &Ps);

    Word sign = -1;
    Word Ch1 = Ch;
    bool refined = false;
    while (Ch != NIL) {
        Word C;
        // next sector cell
        ADV(Ch, &C, &Ch1);

        // sample point may need updating if a refinement of its bottom was just performed.
        if (refined) {
            printf("need to set sample of "); LWRITE(LELTI(C, INDX)); SWRITE("\n");
        }

        // no more cells or no more polynomials
        if (Ch1 == NIL || PM == NIL) {
            break;
        }

        // next section -- top of C
        ADV(Ch1, &C, &Ch1);

        // get k-th coordinate of the sample point of C
        Word SQ, SJ;
        GETSAMPLEK(LELTI(C, SAMPLE), &SQ, &SJ);        sign = COMPARE(SQ, &SJ, PM, &PI);

        if (sign < 0) { // S < P, don't refine, keep looking.
            Ch = Ch1;

            continue;
        } else if (sign > 0) { // S > P, refine previous sector, FIRST(Ch)
            Ch1 = RED2(RefineCell(k, Ch, SQ, SJ, PM, PI, PFs, &refined));
        }

        // and get next polynomial
        NextPolynomial(Ps, &PM, &PI, &J, &Ps);
        Ch = Ch1;
    }

    return Ch2;
}

Word QepcadCls::REFINE(Word k, Word D, Word A, Word PF)
{
    PF = RED(PF);
    Word A1;
    ADV(A, &A1, &A); // decomstruct A. A1 is the set of level k+1 polys

    // find the new PO_REFINE polynomials.
    Word Ps = NIL, I = LELTI(D, INDX);
    while (A1 != NIL) {
        Word P;
        ADV(A1, &P, &A1);

        if (LELTI(P, PO_REFINEMENT) == I) { // list equality check not needed, since same cell index pointer is used
            Ps = COMP(LELTI(P, PO_POLY), Ps);
        }
    }

    // finished?
    Word Ch = LELTI(D, CHILD);
    if (A == NIL || Ch == NIL) {
        return D;
    }

    // do refinement if the list of Ps is non-empty
    if (Ps != NIL) {
        Ch = RefineSubcad(k, Ch, Ps, PF);
        SLELTI(D, CHILD, Ch);
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

