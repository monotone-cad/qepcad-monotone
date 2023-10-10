/*======================================================================
 * s = SacPolyToMapleString(r, P)
 *
 * converts a saclib polynomial P, in r variables, into a char array.
 *====================================================================*/
#include "qepcad.h"
#include <iostream>
#include <string>
using namespace std;

string WriteMonomial(Word r, Word D, Word V)
{
    bool write_mul = false;
    string term = "";

    while (D != NIL) {
        Word e, v;
        ADV(D, &e, &D);
        ADV(V, &v, &V);

        // skip variables of 0 degree
        if (e == 0) {
            continue;
        }

        // write multiplication sign
        if (write_mul) {
            term += " * ";
        }

        write_mul = true;

        // write variable...
        term += string(1, v);

        // ... and the exponent if needed
        if (e > 1) {
            term += "^" + to_string(e);
        }
    }

    return term;
}

string SacPolyToMapleString(Word r, Word P, Word V)
{
    // put polynomial in distributed form.
    // list of pairs (a, E) where a in Z is the coefficient and E is a list (e1,...,er) of degrees in each variable
    Word P1 = DIPFP(r, P);

    Word a, D;
    ADV2(P1, &a, &D, &P1);

    string out = "";
    while (true) {
        if (a != 1) {
            out += to_string(a) + " ";
        }

        out += WriteMonomial(r, D, V);

        if (P1 == NIL) {
            break;
        }

        ADV2(P1, &a, &D, &P1);
        // write sign
        if (a < 0) {
            out += " - ";
            a = -a;
        } else {
            out += " + ";
        }
    }

    return out;
}

