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
    bool a_alg = PDEG(Ma) == 1;
    bool b_alg = PDEG(Mb) == 1;

    // rational vs rational
    if (a_alg && b_alg) {
        Word a = ISLIST(Ba) ? FIRST(Ba) : 0;
        Word b = ISLIST(Bb) ? FIRST(Bb) : 0;

        return RNCOMP(a, b);
    }

    if (a_alg) {
        // b is rational
        return AFCOMP(Ma, Ia, Ba, Bb);
    }

    if (b_alg) {
        // a is rational
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

// set cell C sample point, element k to value a
// it overrides the base polynomial and isolating interval with new ones, it is assumed that these are sufficient to
// define the cell.
void SETSAMPLEK(Word C, Word k, Word M, Word I, Word b)
{
    // recursive update children.
    Word Ch = LELTI(C, CHILD), C1;
    while (Ch != NIL) {
        ADV(Ch, &C1, &Ch);

        SETSAMPLEK(C1, k, M, I, b);
    }

    Word S = LELTI(C, SAMPLE);
    Word M1, I1, b1;
    FIRST3(S, &M1, &I1, &b1);

    // k-th coordinate
    SLELTI(b1, k, b);

    // set sample with new minimal polynomial and isolating interval, if reqired.
    if (M == NIL) M = M1;
    if (I == NIL) I = I1;

    SLELTI(C, SAMPLE, LIST3(M, I, b1));
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
    // C1 updates to the left-most end of the isolating interval J (same as in CAD construction)
    SETSAMPLEK(C1, k, NIL, NIL, c);
    // C2 is set to the new sample provided.
    SETSAMPLEK(C2, k, M, I, b);
    // C3 is updated later

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
        Word Sb1 = THIRD(LELTI(C1, SAMPLE));
        Word b = AFFINT(RNCEIL(Jr) + 1);

        SETSAMPLEK(C1, k, NIL, NIL, b);
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

        SETSAMPLEK(C1, k, NIL, NIL, b);
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

