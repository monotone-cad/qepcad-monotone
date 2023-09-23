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

// compare two coordinates of sample points.
Word COMPARE(Word M1, Word I1, Word M2, Word I2)
{
    //LWRITE(M1); SWRITE(" ");
    //LWRITE(I1); SWRITE(" ");
    //LWRITE(M2); SWRITE(" ");
    //LWRITE(I2); SWRITE("\n");

    bool a_rat = PDEG(M1) == 1;
    bool b_rat = PDEG(M2) == 1;

    // easy case: rational comparison
    if (a_rat && b_rat) {
        return RNCOMP(FIRST(I1), FIRST(I2));
    }

    // otherwise, at least one of a and b is algebraic. compare using rational approximation.
    Word a1, a2, b1, b2;
    while(1) {
        FIRST2(I1, &a1, &a2);
        FIRST2(I2, &b1, &b2);

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
        if (!a_rat) I1 = IUPIIR(M1, LIST2(a1,a2));
        if (!b_rat) I2 = IUPIIR(M2, LIST2(b1,b2));
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

void SETSAMPLEALG(Word C, Word Q, Word J, Word M, Word I, Word b, Word PFs)
{
    Word S;
    if (Q == NIL) {
        S = LIST3(M, I, b);
    } else {
        S = LIST5(Q, J, M, I, b);
    }

    SLELTI(C, SAMPLE, S);
}

void SETSAMPLERAT(Word C, Word M, Word I, Word b, Word PFs)
{
    Word S = LIST3(M, I, b);
    SLELTI(C, SAMPLE, S);

    Word Ch = LELTI(C, CHILD);
    if (Ch == NIL) return;

    // head and tail of PFs
    Word PF;
    ADV(PFs, &PF, &PFs);

    // one child, just set it to zero.
    if (LENGTH(Ch) == 1) {
        Word b1 = CONC(LCOPY(b), LIST1(AFFINT(0)));
        SETSAMPLERAT(FIRST(Ch), M, I, b1, PFs);

        return;
    }

    // otherwise, we have roots to find and sample poits to update.
    Word SP = NIL;
    Word r = LELTI(C, LEVEL) + 1;
    while (PF != NIL) {
        Word P1;
        ADV(PF, &P1, &PF);

        Word PP = SUBSTITUTE(r, LELTI(P1, PO_POLY), S, true); // with rational coefficients
        SP = COMP(PP, SP);
    }

    // root finding.
    Word B = ROOTS(SP);

    Word C1, C2, JR, MR;
    while (Ch != NIL && B != NIL) {
        ADV2(Ch, &C1, &C2, &Ch);
        ADV2(B, &JR, &MR, &B);

        // sector cell C1 takes on left hand end of isolating interval
        Word b1 = CONC(LCOPY(b), LIST1(AFFRN(FIRST(JR))));
        SETSAMPLERAT(C1, M, I, b1, PFs);

        // for the sector cell C2, set sample point to the root.
        if (PDEG(MR) == 1) { // rational root
            Word c = AFFRN(IUPRLP(MR));

            Word b2 = CONC(LCOPY(b), LIST1(c));
            SETSAMPLERAT(C2, M, I, b2, PFs);
        } else { // algebraic root, sample point is now algebraic
            Word c = AFGEN();

            Word b2 = CONC(LCOPY(b), LIST1(c));
            SETSAMPLEALG(C2, NIL, NIL, MR, JR, b2, PFs);
        }
    }

    // last cell, integer number, ceiling of right endpoint of last isolating interval + 1
    Word b1 = CONC(LCOPY(b), LIST1(AFFINT(RNCEIL(SECOND(JR)) + 1)));
    SETSAMPLERAT(FIRST(Ch), M, I, b1, PFs);
}

// are the first n-1 coordinates of S rational
bool IsRat(Word SQ, Word SJ, Word SM, Word SI, Word Sb)
{
    // base polynomial is rational
    if (PDEG(SM) == 1) {
        return true;
    }

    // otherwise, SM represents an algebraic extension.
    if (SQ != NIL) {
        // at least one of (c_1,...,c_{k-1}) has to be algebraic
        return false;
    }

    // otherwise, the point is in primitive form. we must check whether all c_1,...,c_{k-1} are all rational.
    Word k = LENGTH(Sb) - 1;
    while (Sb != NIL && k > 0) {
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

// let C be a level k cell. set sample, level k to (M,I,b), updating children as appropriate.
void SETSAMPLE(Word C, Word M, Word I, Word PFs)
{
    Word k = LELTI(C, LEVEL);
    Word S = LELTI(C, SAMPLE);

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

    // is the new sample rational?
    Word S1;
    if (PDEG(M) == 1) {
        // then append it
        Word c = LIST2(FIRST(I), PMON(1,0));
        Sb = CONC(Sb, LIST1(c));

        S1 = LIST3(SM, SI, Sb);
    } else {
        // otherwise (algebraic), store in extended form
        S1 = LIST5(AFPFIP(1,M), I, SM, SI, Sb);
    }

    SLELTI(C, SAMPLE, S1);
}

// let C = FIRST(Cs) be a (0,...,0,1)-cell and s be a point in C. refine C into three new cells such that s is a new
// (0,...,0,0)-cell
Word RefineCell(Word k, Word Cs, Word SQ, Word SJ, Word PM, Word PI, Word c, Word PFs)
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
    LWRITE(PM); SWRITE(" "); LWRITE(PI); SWRITE("\n");

    // update indices
    SETINDEXK(C2, k, ++j);
    SETINDEXK(C3, k, ++j);

    // update sample
    // existing sample point of C1 will be correct for one of the cells. check.
    Word sign = COMPARE(SQ, SJ, PM, PI);
    // -1: C1 is correct, 0: C2 is correct, +1: C3 is correct.

    if (sign != -1) { // need to update C1 ...
        // ... to the left-hand end of the isolating interval
        SETSAMPLE(C1, PMON(1,1), LIST1(c), PFs);
    }

    if (sign != 0) { // need to update C2 ...
        // to the new "refinement point" given
        SETSAMPLE(C2, PM, PI, PFs);
    }

    // update indices of remaining cells.
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

// univariate algebraic number polynomial strip, assuming algebraic coordinate does not appear.
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

// find some rational point between the last coordinates of P1 = (M1, I1, b1) and P2 = (M2, I2, b2)
void MIDPOINT(Word M1, Word I1, Word b1, Word M2, Word I2, Word b2, Word* M_, Word* I_, Word* b_)
{
    Word k = LENGTH(b2);
    Word af2 = LAST(b2);

    // if a is rational...
    Word a1;
    if (AfIsRat(b1, &a1)) {
        // extract rational value
        Word a2;
        Word c;

        // ... and b2 is also rational
        if (AfIsRat(af2, &a2)) {
            // simply extract the rational value and calculate the midpoint between a1 and a2
            c = RNQ(RNSUM(a1, a2), RNINT(2));
        } else { // .. otherwise, b2 is algebraic
            // extract left-hand endpoint of isolating interval.
            a2 = FIRST(I2);

            // if a1 < a2, a2 can be used
            if (RNCOMP(a1, a2) < 0) {
                c = a2;
            } else {
                perror("a1 > a2. need to calculate a new \"midpoint\".\n");
                c = 0;
            }
        }

        // return M, I, (a1 + a2) / 2
        *M_ = LCOPY(M1);
        *I_ = LCOPY(I1);

        Word b = LCOPY(b2);
        SLELTI(b, k, AFFRN(c));
        *b_ = b;

        return;
    }

    // otherwise a1 is algebraic
    // TODO
}

// Refine subcad D to be compatible with level 1 polynomials Ps
Word RefineSubcad(Word k, Word Ch, Word Ps, Word PFs)
{
    Word Ch2 = Ch; // backup of list, to return.
    Word PM, PI, J;

    // Ps is a list of sample points, in ascending order, which define new 0-cells we will add to D
    Word Ch1 = Ch;
    while (Ps != NIL) {
        NextPolynomial(Ps, &PM, &PI, &J, &Ps);

        // find child cell which needs splitting.
        // i.e., the sector cell whose top is aftor the current refinement point, or the last cell.
        Word sign = -1, SQ, SJ;
        while (sign < 0) {
            Word C;
            Ch = Ch1;
            ADV(Ch, &C, &Ch1);

            // last cell?
            if (Ch1 == NIL) {
                break;
            }

            // top of sector is past refinement ponit
            Word S, SM, SI, Sb, junk;
            ADV(Ch1, &C, &Ch1);
            S = LELTI(C, SAMPLE);
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

            sign = COMPARE(SQ, SJ, PM, PI);
        }

        // refinement point coincides with an existing 0-cell
        if (sign == 0) continue; // point P coincides with an existing section cell. no refinement

        // change of sign occurs in cell Cp = FIRST(Ch). refine it
        Ch1 = RefineCell(k, Ch, SQ, SJ, PM, PI, AFFRN(FIRST(J)), PFs);
    }

    // final step: update sample point of last refined cell.
    Ch1 = RED2(Ch1);

    Word Jl, Jr;
    FIRST2(J, &Jl, &Jr);
    Word C1 = FIRST(Ch1);
    Word M1, I1, Sb1;
    FIRST3(LELTI(C1, SAMPLE), &M1, &I1, &Sb1);

    // then it's the last sector
    if (LENGTH(Ch1) == 1) {
        Word b = AFFINT(RNCEIL(Jr) + 1);
        SETSAMPLE(C1, PMON(1,1), LIST1(b), PFs);

        return Ch2;
    }

    // maybe we will get lucky and the sample point of the original cell is already correct.
    //if (COMPARE(PM, PI, Pb, M1, I1, LAST(Sb1)) < 0) {
    //    printf("we got lucky, sample point is already correct.\n");

    //    return  Ch2;
    //}

    // otherwise, the sample point needs updating. it's bounded from the right by a cell C2.
    //Word C2, SM, SI, Sb;
    //C2 = SECOND(Ch1);
    //Word S = LELTI(C2, SAMPLE);
    //FIRST3(S, &SM, &SI, &Sb);

    //Word M, I, b;
    //MIDPOINT(PM, PI, Pb, SM, SI, Sb, &M, &I, &b);
    //LWRITE(LELTI(C1, INDX)); SWRITE(" "); LWRITE(b); SWRITE("\n");
    //SETSAMPLERAT(C1, M, I, b, PFs);

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

