/*======================================================================
 * L2 <- LDCOPY(L1)
 * List deep (recursive) copy.
 *
 * Input:
 *     L1 : list
 * Output:
 *     L2 : copy of L1
 *
 *====================================================================*/
#include "qepcad.h"

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

