#include "qepcad.h"

// is algebraic field element a rational number? return true if rational and store the rational number in r_
// otherwise, return false, if algebraic then r_ is not modified
bool AfIsRat(Word b, Word* r_)
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

// convenience function: get k-th coordinate of the sample point as an algebraic number
// pass -1 for last
void GETSAMPLEK(Word k, Word S, Word* Q_, Word* J_)
{
    Word SQ, SJ, SM, SI, Sb, junk;
    bool extended = LENGTH(S) == 5;

    if (extended && (k == -1 || k == LENGTH(LELTI(S, 5)))) {
        // extended representation, last coordinate
        FIRST4(S, &SQ, &SJ, &SM, &SI);
        SQ = AFPNORM(1, SM, SQ);
    } else {
        if (extended) {
            FIRST5(S, &junk, &junk, &SM, &SI, &Sb);
        } else {
            FIRST3(S, &SM, &SI, &Sb);
        }

        if (k == -1) k = LENGTH(Sb);
        Sb = LELTI(Sb, k);
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


