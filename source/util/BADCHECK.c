/*======================================================================
 * isBad <- BADCHECK(D, P)
 * Returns true if the CAD contains no bad cells of dimension greater than 1.
 *
 * Input:
 *     d : dimension of D.
 *     D : CAD cell
 * Output:
 *     isBad : true if the CAD contains any bad cells of dimension > 1
 *
 *====================================================================*/
#include "qepcad.h"

bool BadCheckHelper(Word dim, Word D)
{
    // base case: dimension 2, check for vanishing polynomial
    if (dim == 2) {
        Word s1 = FIRST(LELTI(D, SIGNPF));

        return MEMBER(0, s1);
    }

    // walk the CAD, proceed by induction
    bool sector = false;
    Word Ch = LELTI(D, CHILD);
    if (Ch == NIL) return false;

    while (Ch != NIL) {
        sector = !sector;
        Word C;
        ADV(Ch, &C, &Ch);

        // recurse
        if (BadCheckHelper(dim + sector, C)) {
            return true;
        }
    }

    return false;
}

Word BADCHECK(Word D)
{
    return BadCheckHelper(0, D);
}

