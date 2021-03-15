/*======================================================================
                      D <- TICAD(Q,F,f,P,A)

Truth Invariant CAD Construction Phase. Recompute based on new (linear) proj factors.

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

// given a list of level k cells, and a list of level k polynomials, ensure the list L is properly sign invariant
Word SPLIT(Word L, Word P);

// moves n items along L
Word ADVN(Word L, Word n);

// -1 if < 0, 1 if > 0, 0 is there is a root within the cell
Word INTERVALSIGN(Word L, Word idx, Word x);

// convert k-th coordinate of sample point to rational
// TODO remove this
// Word SAMPLETORATIONAL(Word C, Word k);
// get the k-the coordinate of sample point of C
Word SAMPLEK(Word C, Word k);


// find the solution of a linear polynomial p
Word LINEARROOT(Word p);

// Append sign of projection factor to C at (C, LEVEL)
void APPENDSIGNPF(Word C, Word S);

// increment index at level k by t
void INCINDEXL(Word C, Word k, Word t);
void INCINDEX(Word C, Word k, Word t);

// Deep cell copy
Word CELLCOPY(Word C);

// Given cell C which is not sign invariant and rational number p which is the point at which new polynomial changes
// sign, return a list of 3 cells having the same truth values as C but which are sign invarient
Word SPLITCELL(Word C, Word p);

// insert new cells into list L
Word INSERTCELLS(Word L, Word NewCells);

Word QepcadCls::RECOMPUTE(Word C, Word Q, Word F, Word f, Word P, Word A)
{
    Word D, k, L, Children, Child, P1, S;

Step1: /* Initialise. */
    D = C;

Step2: /* Recursive cad walk, looking for proj factors that need fixing. (Setup) */
    k = LELTI(C, LEVEL);
    Children = LELTI(C, CHILD);
    if (Children == NIL) goto Return;

    L = Children;
    P1 = LELTI(P, k+1);

Step3: /* Recurse on children. */
    while (Children != NIL) {
        ADV(Children, &Child, &Children);

        S = FIRST(LELTI(Child, SIGNPF)); // Note SIGNPF comes in reverse order
        if (LENGTH(S) < LENGTH(P1)) {
            Children = SPLIT(L, P1);

            // goto Step3; // go back to the beginning to ensure all cells are OK
            goto Return;
        }

        Child = RECOMPUTE(Child,Q,F,f,P,A);
    }

Return: /* return. */
    return D;
}

Word SPLIT(Word L, Word Ps)
{
    Word P, p, x, s, C, k;

    if (LENGTH(L) < 2) goto Return;

    C = FIRST(L);
    k = LELTI(C, LEVEL);
    Ps = ADVN(Ps, LENGTH(LELTI(LELTI(C, SIGNPF), k)));

Step1: /* Iterate through the new polynomials */
    while (Ps != NIL) {
    ADV(Ps, &P, &Ps);
    p = LELTI(P, PO_POLY);
    x = LINEARROOT(p);

    if (x == NIL) {
        return L;
    }

    // convert to rationar representation, same as sample points
    x = AFFRN(x);
    AFWRITE(x, 'x');
    Word NewCells = NIL; // list of lists of new cells.
    for (int i = 1; i <= LENGTH(L); i++) {
        s = INTERVALSIGN(L, i, x);

        // not sign invariant
        if (s == 0) {
            NewCells = SPLITCELL(LELTI(L, i), x);

            continue;
        }

        APPENDSIGNPF(LELTI(L, i), s);
    }

    // add new cells (Ds) to list, deleting C
    L = INSERTCELLS(L, NewCells);

    }
Return:
    return L;
}

Word ADVN(Word L, Word n)
{
    Word v;

    for (int i = 0; i < n; i++) {
        if (L == NIL) break;
        ADV(L, &v, &L);
    }

    return L;
}

Word LINEARROOT(Word p)
{
    Word n = PDEG(p);
    if (n > 1) return NIL; // TODO what to do if it's not linear

    // directly solve P := den X - num = 0
    Word den = PLDCF(p);
    Word num = PCOEFF(p,0);

    // x = 0
    if (num == 0) {
        return 0;
    }

    Word lm, ln;
    IFCL2(den,&lm,&ln);
    if (ICOMP(lm,ln) != 0) return NIL;

    Word x = RNRED(INEG(num),den);
    return x;
}

// convert k-th coordinate of sample point to rational
Word __SAMPLETORATIONAL(Word C, Word k)
{
    Word S, I, p, B;

    S = LELTI(C, SAMPLE);

    if (LENGTH(S) != 3) return 0;

    I = LELTI(S, 3);
    p = LELTI(I, k);
    if (p == 0) return 0;

    // compute B* (b)
    // TODO should be B*b(0)
    B = IUPBREV(SECOND(p), FIRST(p));
    return B;
}

Word SAMPLEK(Word C, Word k)
{
    Word S, P, M, I;

    S = LELTI(C, SAMPLE);

    if (LENGTH(S) != 3) return 0;

    FIRST3(S, &M, &I, &P);
    return LIST3(LELTI(P, k), M, I);
}

Word INTERVALSIGN(Word L, Word idx, Word x)
{
    Word C, k, Cl, Cr, Sl, Sr, l, r, Ml, Il, Mr, Ir;

Step1: /* Initialise */
    C = LELTI(L, idx);
    k = LELTI(C, LEVEL);

Step2: /* special case for even cell */
    // next and previous cells
    Sl = NIL; Sr = NIL;

    if (EVEN(idx)) {
        Sl = SAMPLEK(C, k);
        Sr = Sl;

        goto Step3;
    }

    if (idx > 1) {
        Cl = LELTI(L, idx - 1);

        Sl = SAMPLEK(Cl, k);
    }

    if (idx < LENGTH(L)) {
        Cr = LELTI(L, idx + 1);

        Sr = SAMPLEK(Cr, k);
    }

    // special case for leftmost and rightmost cells (whele l == NIL and r == NIL resp).
    if (Sl == NIL) Sl = Sr;
    if (Sr == NIL) Sr = Sl;

Step3: /* endpoint comparison */
    // get point, M, I for each point
    FIRST3(Sl, &l, &Ml, &Il);
    FIRST3(Sr, &r, &Mr, &Ir);

    // endpoint comparisons
    Word xl = AFCOMP(Ml, Il, x, l); // x is: left < equal < right of l
    Word xr = AFCOMP(Mr, Ir, x, r); // x is: less < equal < right of r

    if (xl > 0 && xr < 0) return 0;
    if (xr <= 0) return 1;
    if (xl >= 0) return -1;
    return NIL; // if this happened something went badly wrong
}

void APPENDSIGNPF(Word C, Word sign)
{
    Word S, M, D, l;

    S = FIRST(LELTI(C, SIGNPF));
    l = LENGTH(S); // how many SIGNPF
    M = LELTI(C, MULSUB);
    D = LELTI(C, DEGSUB);

    // append sign
    CONC(S, LIST1(sign));

    // append new degsub (degsub = [d_1, ..., d_n], degrees of k-level projection factors
    // we know it's a linear polynomial
    if (D != NIL) { CONC(D, LIST1(1)); }

    // append multiplicity - ((i_1,e_1),...,(i_n,e_n)) where i_j is the index at level k and e_j is the multiplicity
    if (sign == 0) { // only applies to section cells
        // TODO does order matter?
        if (M != NIL) {
            CONC(M, LIST1(LIST2(l+1, 1)));
        } else {
            SLELTI(C, MULSUB, LIST1(LIST2(l+1, 1)));
        }
    }
}

void INCINDEX(Word C, Word k, Word t)
{
    Word I = LELTI(C, INDX);
    SLELTI(I, k, LELTI(I, k) + t);
}

void INCINDEXL(Word L, Word k, Word t)
{
    Word C;

    while (L != NIL) {
        ADV(L, &C, &L);

        INCINDEX(C, k, t);
    }
}

Word CELLCOPY(Word C)
{
    Word Children, Ch, NewChildren;

    // TODO figure out which ones should be LLCOPied
    Children = LELTI(C, CHILD);
    NewChildren = NIL;
    while (Children != NIL) {
        ADV(Children, &Ch, &Children);

        NewChildren = COMP(
            CELLCOPY(Ch),
            NewChildren
        );
    }

    if (NewChildren != NIL) NewChildren = INV(NewChildren);

    return MCELL(
        LELTI(C, LEVEL),
        NewChildren,
        NIL,
        LELTI(C, TRUTH),
        LCOPY(LELTI(C, SAMPLE)),
        LCOPY(LELTI(C, INDX)),
        LLCOPY(LELTI(C, SIGNPF)),
        LELTI(C, HOWTV),
        LCOPY(LELTI(C, DEGSUB)),
        LCOPY(LELTI(C, MULSUB))
    );
}

Word SPLITCELL(Word C, Word p)
{
    Word Cl, Cr, k;

    k = LELTI(C, LEVEL);
    Cl = CELLCOPY(C);
    Cr = CELLCOPY(C);

    // Update indices
    INCINDEX(C, k, 1);
    INCINDEX(Cr, k, 2);

    // add signs of new projection factor
    APPENDSIGNPF(Cl, -1);
    APPENDSIGNPF(C, 0);
    APPENDSIGNPF(Cr, 1);

    return LIST3(Cl, C, Cr);
}

Word INSERTCELLS(Word LL, Word NewCells)
{
    Word C, k, i, t, L, LRed, NewRed;

    C = FIRST(NewCells);
    NewRed = RED(NewCells);
    k = LELTI(C, LEVEL);
    i = LELTI(LELTI(C, INDX), k);

    L = ADVN(LL, i-1); // L = [C_i, ...rest]
    LRed = RED(L);

    // increment LRed indices
    t = LENGTH(NewCells) - 1;
    INCINDEXL(LRed, k, t);

    // concotenate new cells onto remaining cells, and replace first by first new cell
    SRED(L, CONC(NewRed, LRed));
    SFIRST(L, C);

    return LL;
}

