/*======================================================================
 * s = SacPolyToMapleString(r, P)
 *
 * converts a saclib polynomial P, in r variables, into a char array.
 *====================================================================*/
#include "qepcad.h"

inline Word AsciiInt(Word a) {
    return a + 48;
}

inline Word AsciiChar(Word a) {
    return a + 96;
}

void IntToString(Word A, Word* S)
{
    Word a, b;
    bool sign = false;
    a = A;
    if (a < 0) {
        sign = true;
        a = -a;
    }

    do {
        b = a % 10;
        a = a / 10;

        *S = COMP(AsciiInt(b), *S);
    } while (a > 10);

    if (sign) {
        *S = COMP('-', *S);
    }
}

void SacPolyToStringHelper(Word r, Word P, Word* S)
{
    if (r == 0) {
        IntToString(P, S);

        return;
    }

    Word e, A;
    ADV2(P, &e, &A, &P);

    if (e == 0) {
        return SacPolyToStringHelper(r - 1, A, S);
    }

    *S = COMP(')', *S);
    while (true) {
        // exponent
        if (e > 1) {
            *S = COMP2('^', AsciiInt(e), *S);
        }

        // variable
        if (e > 0) {
            *S = COMP2('*', AsciiChar(r), *S);
        }

        // print A, polynomial in r-1 variables
        SacPolyToStringHelper(r - 1, A, S);

        // done?
        if (P == NIL) {
            break;
        }

        *S = COMP('+', *S);
        ADV2(P, &e, &A, &P);
    }

    *S = COMP('(', *S);
}

const char* SacPolyToMapleString(Word r, Word P)
{
    // convert to distributed form
    LWRITE(P); SWRITE("\n");
    Word S = NIL;
    SacPolyToStringHelper(r, P, &S);
    CLOUT(S); SWRITE("\n");
    return "x";
}

