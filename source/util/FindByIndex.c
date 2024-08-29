#include "qepcad.h"

// return the cell in list L with (partial) index I.
// purpose is to find parent cells.
// L : list of cells, level k
// I : partial index to match, level j
// return : C s.t. C has index I
Word FindByIndex(Word L, Word I, Word j, Word k)
{
    Word J, C, C1;

    while (L != NIL) {
        ADV(L, &C, &L);
        J = LELTI(C, INDX);

        // no partial match
        if (FIRST(I) != LELTI(J, k)) continue;

        // otherwise, partial match. if j == k then we're done
        if (j == k) return C;

        // if j < k, then we have a partial but incomplete match, continue searching recursively
        return FindByIndex(LELTI(C, CHILD), RED(I), j, k + 1);
    }

    // no partial match.
    return NIL;
}


