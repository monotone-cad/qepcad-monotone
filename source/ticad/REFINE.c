/*======================================================================
                      D <- REFINE(k,D,A)

Refine CAD D so that its zero-cells C0 are compatible with relevant new PO_REFINEMENT polynomials in A

\Input
  \parm{k} is the level of D
  \parm{D} is the original cad
  \parm{A} is the list~$(A_1,\ldots,A_r)$,
           where $A_i$ is the list of all
           the distinct $i$--level normalized input polynomials.

Output
  \parm{D} is the refinement of D compatible with new P_MONOTONE polynamials in A

======================================================================*/
#include "qepcad.h"

// need a new property of cells, SIGNRP, like SIGNPF. this wiss contain signs of refinement polynomials, level k, above
// the sub-cad as indicated in the polynomial. add in cell write.

// compare two coordinates of sample points.
Word COMPARE(Word Ma, Word Ia, Word Ba, Word Mb, Word Ib, Word Bb)
{
    bool a_rat = !(ISLIST(Ba) && PDEG(SECOND(Ba)) > 0);
    bool b_rat = !(ISLIST(Ba) && PDEG(SECOND(Ba)) > 0);

    // rational vs rational
    if (a_rat && b_rat) {
        Word a = ISLIST(Ba) ? FIRST(Ba) : 0;
        Word b = ISLIST(Bb) ? FIRST(Bb) : 0;

        return RNCOMP(a, b);
    }

    if (b_rat) {
        // a is algebraic while b is rational
        return AFCOMP(Ma, Ia, Ba, Bb);
    }

    if (a_rat) {
        // b is algebraic while b is rational
        return AFCOMP(Mb, Ib, Ba, Bb);
    }

    // otherwise both are algebraic and they might have different generators
    perror("TODO alg vs alg. need to get the same generator.\n");

    return 0;
}

// list deep copy
// like LCOPY but it recursively copies all the sublists.
Word LDCOPY(Word L)
{
    // base: empty
    if (L == NIL) return NIL;

    // base: non-list
    if (!ISLIST(L)) return L;

    // recursive: list
    Word LL = NIL, A;
    while (L != NIL) {
        ADV(L, &A, &L);

        LL = COMP(LDCOPY(A), LL);
    }

    return INV(LL);
}

// set cell C index, element k to value a
// child recursive update not needed.
void SETINDEXK(Word C, Word k, Word a)
{
    Word I = LELTI(C, INDX);
    SLELTI(I, k, a);
}

void SETSAMPLERAT(Word C, Word M, Word I, Word b)
{
    SLELTI(C, SAMPLE, LIST3(M, I, b));

    Word Ch = LELTI(C, CHILD);
    while (Ch != NIL) {
        Word C1;
        ADV(Ch, &C1, &Ch);

        Word b1 = CONC(LCOPY(b), LIST1(AFFINT(42)));
        SETSAMPLERAT(C1, M, I, b1);
    }
}

void SETSAMPLEALG(Word C, Word Q, Word J, Word M, Word I, Word b)
{
    Word S;
    if (Q == NIL) {
        S = LIST3(M, I, b);
    } else {
        S = LIST5(Q, J, M, I, b);
    }

    SLELTI(C, SAMPLE, S);
}

inline bool IsAF(Word b)
{
    return b != 0 && PDEG(SECOND(b)) > 0;
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
        if (IsAF(b)) {
            return false;
        }
    }

    // otherwise first k-1 coordinates are ratinoal
    return true;
}

// let C be a level k cell. set sample, level k to (M,I,b), updating children as appropriate.
void SETSAMPLE(Word C, Word M, Word I, Word b)
{
    Word k = LELTI(C, LEVEL);
    Word S = LELTI(C, SAMPLE);

    // extended or primitive representation?
    Word SQ, SJ, SM, SI, Sb;
    if (LENGTH(S) == 5) { // extended (we don't care about the last coordinate)
        FIRST5(S, &SQ, &SJ, &SM, &SI, &Sb);

        Sb = CONC(Sb, LIST1(b));
    } else { // primitive
        SQ = NIL, SJ = NIL;
        FIRST3(S, &SM, &SI, &Sb);

        SLELTI(Sb, k, b);
    }

    bool is_rat = IsRat(SQ, SJ, SM, SI, Sb);
    bool b_is_af = IsAF(b);

    // completely algebraic
    if (!b_is_af) {
        if (is_rat) { // rational base
            printf("rat rat\n");
            SETSAMPLERAT(C, PMON(1,1), LIST2(0,0), Sb);

            return;
        }

        // algebaric base
        printf("base alg, b rat\n");
        SETSAMPLEALG(C, NIL, NIL, SM, SI, Sb);

        return;
    }

    // otherwise, b is algebraic
    if (is_rat || (EQUAL(SM, M) && EQUAL(SI, I))) {
        // rational base, or extensions are equal, just use M and I
        printf("base rat or same, b alg\n");
        SETSAMPLEALG(C, NIL, NIL, SM, SI, Sb);

        return;
    }

    // otherwise, S and b are both algebraic, and they both have different extensions.
    // write new sample in extended form
    if (SQ != NIL) { // if it was in primitive form, remove the last coordinate.
        SLELTI(Sb, k, NIL);
    }

    printf("base alg, b alg\n");
    SETSAMPLEALG(C, M, I, SM, SI, Sb);
}

// let C = FIRST(Cs) be a (0,...,0,1)-cell and s be a point in C. refine C into three new cells such that s is a new
// (0,...,0,0)-cell
Word RefineCell(Word k, Word Cs, Word M, Word I, Word b, Word c)
{
    Word Cs2 = Cs;

    // Let C = (a,b). C becomes (a,s), C2 becomes new cell s and C3 new cell (s,b)
    Word C1;
    ADV(Cs, &C1, &Cs);
    Word C2 = LDCOPY(C1);
    Word C3 = LDCOPY(C1);

    // k-th element of C index.
    Word j = LELTI(LELTI(C1, INDX), k);

    // update indices
    SETINDEXK(C2, k, ++j);
    SETINDEXK(C3, k, ++j);

    // update sample
    // existing sample point of C1 will be valid for one of the cells. check.
    Word SM, SI, Sb;
    FIRST3(LELTI(C1, SAMPLE), &SM, &SI, &Sb);
    Sb = LAST(Sb);

    Word sign = COMPARE(SM, SI, Sb, M, I, b);
    // -1: C1 is correct, 0: C2 is correct, +1: C3 is correct.

    if (sign != -1) { // need to update C1 ...
        // ... to the left-hand end of the isolating interval
        printf("update C1 ");
        SETSAMPLE(C1, M, I, c);
    }

    if (sign != 0) { // need to update C2 ...
        // to the new "refinement point" given
        printf("update C2 ");
        SETSAMPLE(C2, M, I, b);
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

// convenience function for next polynomial
void NextPolynomial(Word Ps, Word* PM_, Word* PI_, Word* Pb_, Word* J_, Word* Ps_)
{
    Word P;
    ADV(Ps, &P, Ps_);
    FIRST3(FIRST(P), PM_, PI_, Pb_);
    *J_ = SECOND(P);
}

// Refine subcad D to be compatible with level 1 polynomials Ps
Word RefineSubcad(Word k, Word Ch, Word Ps)
{
    Word Ch2 = Ch; // backup of list, to return.
    Word PM, PI, Pb, J;

    // Ps is a list of sample points, in ascending order, which define new 0-cells we will add to D
    Word Ch1 = Ch;
    while (Ps != NIL) {
        NextPolynomial(Ps, &PM, &PI, &Pb, &J, &Ps);
        Pb = LAST(Pb);

        // find child cell which needs splitting.
        // i.e., the sector cell whose top is aftor the current refinement point, or the last cell.
        Word sign = -1;
        while (Ch1 != NIL && sign < 0) {
            Word C;
            Ch = Ch1;
            ADV(Ch, &C, &Ch1);

            // last cell?
            if (Ch1 == NIL) {
                break;
            }

            // top of sector is past refinement ponit
            ADV(Ch1, &C, &Ch1);
            Word SM, SI, Sb, junk;
            FIRST3(LELTI(C,SAMPLE), &SM, &SI, &Sb);
            Sb = LAST(Sb);

            sign = COMPARE(SM, SI, Sb, PM, PI, Pb);
        }

        if (sign == 0) continue; // point P coincides with an existing section cell. no refinement

        // change of sign occurs in cell Cp = FIRST(Ch). refine it
        printf("Refine cell, level %d ", k); LWRITE(LELTI(FIRST(Ch), INDX));
        Ch1 = RefineCell(k, Ch, PM, PI, Pb, AFFRN(FIRST(J)));
    }

    // final step: update sample point of last refined cell.
    Ch = RED2(Ch);
    Word Jl, Jr;
    FIRST2(J, &Jl, &Jr);

    if (LENGTH(Ch) == 1) { // updating last cell
        Word C1 = FIRST(Ch);
        Word M1, I1, Sb1;
        FIRST3(LELTI(C1, SAMPLE), &M1, &I1, &Sb1);
        Word b = AFFINT(RNCEIL(Jr) + 1);

        SETSAMPLE(C1, M1, I1, b);
        return Ch2;
    }

    // update some intermediate cell C1, bounded from the right by C2
    Word C1, C2, SM, SI, Sb;
    FIRST2(Ch, &C1, &C2);
    Word S = LELTI(C2, SAMPLE);
    FIRST3(S, &SM, &SI, &Sb);

    if (PDEG(SM) == 1 && PDEG(PM) == 1) { // rational sample.
        // find the midpoint between P and C2.
        Word a = FIRST(LAST(Sb)); // rational part of coordinate k
        Word b = AFFRN(RNQ(RNSUM(Jr, a), RNINT(2)));

        SETSAMPLE(C1, SM, SI, b);
    } else {
        // otherwise, C2->SAMPLE is algebraic.
        perror("algebraic refinement last cell. what should the sample be?!\n");
    }

    return Ch2;
}

Word QepcadCls::REFINE(Word k, Word D, Word A)
{
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
        // TODO return new children cells
        SWRITE("Refine subcad "); LWRITE(I); SWRITE("\n");
        Ch = RefineSubcad(k, Ch, Ps);
        SLELTI(D, CHILD, Ch);
    }

    // walk the CAD, sections only.
    Word C, junk;
    ADV(Ch, &junk, &Ch);
    while (Ch != NIL) {
        ADV2(Ch, &C, &junk, &Ch);

        C = REFINE(k+1, C, A);
    }

    return D;
}

