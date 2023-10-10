/*======================================================================
 * V = CreateMapleVariableList(r)
 *
 * generate a variable list for use with Maple functions.
 * Variables will be uppercase letters [A, B, C, ...]
 *====================================================================*/
#include "qepcad.h"
#include <iostream>
#include <string>
using namespace std;

Word CreateMapleVariableList(Word r)
{
    Word i = 0;
    Word L = NIL;
    while (i < r) {
        L = COMP('A' + i, L);

        ++i;
    }

    return L;
}

