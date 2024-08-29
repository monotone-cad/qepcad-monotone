#include "qepcad.h"

// suppose I = (i_1,...,i_n) is a cell index.
// if i_j = i_k = 1 but all other components are 0, then C is 2-dimensional.
// set j and k to indicate which components are equal to one and return 2.
// otherwise, C is not two-dimensional. return d := dim(C) if dim(C) < 2, otherwise d := 3. Values of j and k are undefined.
int TwoDimIndex(Word I, Word *j_, Word *k_)
{
    Word J,d,n,l;

    // d tracks cell dimension.
    n = 0, d = 0;
    while (d <= 3 && I != NIL) {
        ADV(I, &l, &I);
        ++n;

        if (ODD(l)) {
            ++d;
            J = COMP(n, J);
        }
    }

    // if cell is not 2-dimensional, return false
    if (d != 2) {
        return d;
    }

    // otherwise, cell is two-dimensional
    FIRST2(J, k_, j_);

    return 2;
}


