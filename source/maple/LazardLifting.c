/*======================================================================
 * L = LazardLifting(Ps, As)
 *
 * Perform Lazard Lifting. Ps := [g_1,...,g_k,f],
 * where p = { g_1 = 0, ..., g_k = 0 } defines a 0-dimensional cell
 * and f is a function containing a blow-up above p.
 *====================================================================*/
#include "qepcad.h"
#include <iostream>
#include <string>
#include "db/CAPolicy.h"
using namespace std;

// write variable list into maple, variable name VL
void WriteVariableList(MKernelVector kv, Word V)
{
    string vl = "];";
    while (true) {
        Word v;
        ADV(V, &v, &V);

        vl.insert(0, 1, v);

        if (V == NIL) {
            break;
        }

        vl.insert(0, 1, ',');
    }

    EvalMapleStatement(kv, ("VL := [" + vl).c_str());
}

string SatIdealToString(Word rf, Word f, Word rg, Word g, Word rh, Word V)
{
    string out = "I1 := [" + SacPolyToMapleString(rf, f, V) + ", " + SacPolyToMapleString(rg, g, REDI(V, rf - rg))
        + ", 1 - z * f" + to_string(rh)
        + "];";
    cout << out;
    return out;
}

Word SolveSaturatedIdeal(MKernelVector kv, Word rf, Word f, Word rg, Word g, Word rh, Word V)
{
    // Write saturated ideal I1 into Maple
    EvalMapleStatement(kv, SatIdealToString(rf, f, rg, g, rh, V).c_str());

    // compute lex Grobner basis B for I, eliminating new variable z.
    ALGEB B = EvalMapleStatement(
        kv,
        "B1 := remove(depends, Basis(I1, lexdeg([z], VL)), z);"
    );

    if (!IsMapleList(kv, B)) {
        printf("Groebner basis is not a list\n");

        return NIL;
    }

    return NIL;
}

// from the list of polynomials (Grobner basis) in r+1 variables, find polynomials defining refinement points.
// TODO discard variable r+1 and any polynomials depending on it.
// keep if it has nonzero degree in variable r
// eval at p = f1,f2 and find roots.
Word SectionCells(Word r, Word Ps)
{
    Word Ps1 = NIL;
    while (Ps != NIL) {
        Word P, e, A;
        ADV(Ps, &P, &Ps);
        FIRST2(P, &e, &A);

        // if P depends on variable r, discard it
        if (e > 0) continue;

        // refactor evaluation code from MONOTONE, eval P at sample x and append.
        Ps1 = COMP(A, Ps1);
        LWRITE(A); SWRITE("\n");
    }

    return Ps1;
}

Word LazardLifting(MKernelVector kv, Word Ps, Word As)
{
    Word r = LENGTH(Ps);

    if (r > 3) {
        printf("TODO: lift with bad points for r > 3");

        return NIL;
    }

    Word g, f1, f2, f3, h, A1, A2;
    FIRST3(Ps, &f1, &f2, &f3);
    FIRST2(As, &A1, &A2);

    // f3 is blow-up polynomial, p = { f1 = 0, f2 = 0 }
    // construct a 0-dimensinoal ideal (f, g, h); compute saturation w.r.t h, for h in { f1, f2 }

    // first save the variable list, f1 and f2 in maple.
    Word V = CreateMapleVariableList(3);
    WriteVariableList(kv, V);
    EvalMapleStatement(
        kv,
        ("f1 := " + SacPolyToMapleString(1, f1, RED2(V)) + ";").c_str()
    );
    EvalMapleStatement(
        kv,
        ("f2 := " + SacPolyToMapleString(2, f2, RED(V)) + ";").c_str()
    );

    // fix h = f1 and consider g in A2
    // h = 1 - z*f1
    h = IPSUM(
        4,
        PPREPVS(PMON(1,0), 3), // 1
        LIST2(1, PPREPVS(IPNEG(1, f1), 2)) // z * - h
    );
    while (A2 != NIL) {
        ADV(A2, &g, &A2);

        // TODO
        SolveSaturatedIdeal(kv, 3, f3, 2, LELTI(g, PO_POLY), 1, V);
        SectionCells(
            3,
            GVCAP->GROEBNER(LIST3(f3, LELTI(g, PO_POLY), h), LIST3(3,2,4), 4)
        );
    }

    // now fix h = f2 and consider g in A1
    h = IPSUM(
        4,
        PPREPVS(PMON(1,0), 3), // 1
        LIST2(1, LIST2(0, IPNEG(1, f2))) // z * - h
    );
    while (A1 != NIL) {
        ADV(A1, &g, &A1);

        // TODO
        SolveSaturatedIdeal(kv, 3, f3, 1, LELTI(g, PO_POLY), 2, V);
        SectionCells(
            r,
            GVCAP->GROEBNER(LIST3(f3, LELTI(g, PO_POLY), h), LIST3(3,1,4), 4)
        );
    }

    return NIL;
}

