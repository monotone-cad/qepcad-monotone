/*======================================================================
                      D <- FRONTIER(C,Q,F,f,P,A)

Frontier condition (lifting above bad points, lazard-style)

\Input
  \parm{C} is the original cad with missing signs of proj factors (may be recursive i.e., level of C may be > 1)
  \parm{Q} is the list of quantifiers
           in the input formula.
  \parm{F} $=F(x_1,\ldots,x_r)$
           is the normalized input quantifier-free formula.
  \parm{f} is the number of free variables in the input formula.
  \parm{P} is the list~$(P_1,\ldots,P_r)$,
           where $P_i$ is the list of
           $i$--level projection factors.
  \parm{A} is the list~$(A_1,\ldots,A_r)$,
           where $A_i$ is the list of all
           the distinct $i$--level normalized input polynomials.

Output
  \parm{D} is a truth--invariant partial CAD of $f$--space
           in which every leaf cell has a determined truth value.
======================================================================*/
#include "qepcad.h"

void ZEROCELLS(Word C, Word level, Word *L_);

bool VANISHES(Word C, Word k);

Word QepcadCls::FRONTIER(Word C, Word Q, Word F, Word f, Word P, Word A)
{
    // if R < 3, frontier condition is obtaoined automatically.
    if (f == 2) return C;

    // gather 0-cells
    Word Cells = NIL;
    ZEROCELLS(C, 2, &Cells);

    // for each base 0-cell B
    Word B, f1, f2, f3, P1, P2, P3;
    FIRST3(P, &P1, &P2, &P3);
    while (Cells != NIL) {
        ADV(Cells, &B, &Cells);

        Word AP3 = P3; Word i = 0;
        while (AP3 != NIL) {
            ADV(AP3, &f3, &AP3);
            i++; // index of f3

            if (VANISHES(B, i)) {
                printf("f3 at %d vanishes identically\n", i);
            }
        }
    }

    return C;
}

void ZEROCELLS(Word C, Word level, Word *L_) {
    Word Children = LELTI(C, CHILD);

    // no children
    if (Children == NIL) return;

    // isolate level 2 cells
    if (level == 0) {
        printf("appendi "); LWRITE(LELTI(C, INDX)); SWRITE("\n");
        *L_ = COMP(C, *L_);
    }

    Word Ch;
    while (Children != NIL) {
        ADV(Children, &Ch, &Children);

        if (EVEN(LAST(LELTI(Ch, INDX)))) {
            ZEROCELLS(Ch, level-1, L_);
        }
    }
}

// does polynomial at level l, index k vanish identically above C?
bool VANISHES(Word C, Word k) {
    Word Children = LELTI(C, CHILD);

    // no children TODO check what if not full CAD
    if (Children == NIL) return false;

    // for each child of C
    Word Ch; Word i = 0;
    while (Children != NIL) {
        ADV(Children, &Ch, &Children);
        i++; // index (in stack above C) of Ch.

        if (EVEN(i)) continue; // can,t "vanish identically" on a section

        Word S = LELTI(FIRST(LELTI(Ch, SIGNPF)), k);
        if (S == ZERO) return true; // f3 vanishes identically on an interval above C.
    }

    // othervise f3 does not vanish on any intervals.
    return false;
}

