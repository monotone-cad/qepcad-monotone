/*======================================================================
                      D <- FRONTIER(r,C,P)

Frontier condition (lifting above bad points, lazard-style)

\Input
  \parm{r} is the number of free variables in the input formula.
  \parm{C} is the original cad with missing signs of proj factors (may be recursive i.e., level of C may be > 1)
  \parm{P} is the list~$(P_1,\ldots,P_r)$,
           where $P_i$ is the list of
           $i$--level projection factors.

Output
  \parm{D} is a truth--invariant partial CAD of $f$--space
           in which every leaf cell has a determined truth value.
======================================================================*/
#include "qepcad.h"
#include <iostream>

// returns the index of first polynomial which is zero, -1 otherwise.
Word IndexOfFirstZero(Word S)
{
    Word i = 1;
    while (S != NIL) {
        Word s;
        ADV(S, &s, &S);

        if (s == 0) {
            return i;
        }

        ++i;
    }

    return -1;
}

// unique comp, with == check, O(n)
Word UCOMP(Word x, Word L)
{
    Word LL = L;
    while (LL != NIL) {
        Word y;
        ADV(LL, &y, &LL);

        if (x == y) return L;
    }

    return COMP(x, L);
}

Word IdentifyBadCells(Word r, Word C, Word Ps)
{
    Word Children = LELTI(C, CHILD);
    if (Children == NIL) return NIL;

    Word P;
    ADV(Ps, &P, &Ps);

    Word Is = NIL;
    Word Js = NIL;
    bool section = true;
    while (Children != NIL) {
        Word C, Is2, i;
        ADV(Children, &C, &Children);
        section = !section;

        if (section && (Is2 = IdentifyBadCells(r + 1, C, Ps)) != NIL) {
            Word j = IndexOfFirstZero(FIRST(LELTI(C, SIGNPF)));
            Word Pj = LELTI(LELTI(P, j), PO_POLY);

            while (Is2 != NIL) {
                Word I, C1;
                ADV2(Is2, &I, &C1, &Is2);

                Is = COMP2(COMP(Pj, I), C1, Is);
            }
        } else if (!section && r > 2 && (i = IndexOfFirstZero(FIRST(LELTI(C, SIGNPF)))) > 0) {
            Js = UCOMP(i, Js);
        }
    }

    // look up the polynomials from indices in Js
    Word JPs = NIL;
    while (Js != NIL) {
        Word i;
        ADV(Js, &i, &Js);

        JPs = COMP2(LIST1(LELTI(LELTI(P, i), PO_POLY)), C, JPs);
    }

    return CONC(Is, JPs);
}

Word QepcadCls::FRONTIER(Word r, Word C, Word P, Word* A_, Word* J_, Word* RPs_)
{
    // if r < 3, frontier condition is obtaoined automatically.
    if (r < 3) return C;

    // find bad cells
    Word L = IdentifyBadCells(1, C, P);
    while (L != NIL) {
        Word Ps, C, As, f;
        ADV2(L, &Ps, &C, &L);
        // need access to the cell for refinement. store two elements, cell and polynomials.

        Word Hs = LazardLifting(Ps, LELTI(C, SAMPLE), P); // TODO pass S

        // add refinement points and refine cell immediately.
        ADDREFINEMENTPOINTS(LELTI(C, INDX), LELTI(C, SAMPLE), Hs, A_, J_, RPs_);
        C = REFINE(1, C, GVREFL, P);
    }

    return C;
}

