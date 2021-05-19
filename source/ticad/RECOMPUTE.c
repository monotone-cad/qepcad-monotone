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

// -1 if < 0, 1 if > 0, 0 if there is a root within the cell
Word INTERVALSIGN(Word L, Word idx, Word x);

// get the k-th coordinate of sample point of C
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
Word SPLITCELL(Word C, Word p, Word Sl, Word Sr);

// insert new cells into list L
Word INSERTCELLS(Word L, Word NewCells);

// algebraic number (b-a)/2
Word MIDPOINT(Word M, Word l, Word r);

// set sample coordinate, level k, to m
void SETSAMPLEK(Word C, Word k, Word m);

Word QepcadCls::RECOMPUTE(Word C, Word Q, Word F, Word f, Word P, Word A)
{
    Word D, k, L, Children, Child, P1, S;

    /* Initialise. */
    D = C;

    /* Recursive cad walk, looking for proj factors that need fixing. (Setup) */
    k = LELTI(C, LEVEL);
    Children = LELTI(C, CHILD);
    if (Children == NIL) return D;

    L = Children;
    P1 = LELTI(P, k+1);

    /* Recurse on children and figure out if any splitting is needed */
    bool needs_splitting = false;
    int n_p1 = LENGTH(P1);
    while (Children != NIL) {
        ADV(Children, &Child, &Children);

        // we short circuit to prevent doing the length comparison if we already know we need to recompute
        // Note SIGNPF comes in reverse order
        needs_splitting = needs_splitting ||
            LENGTH(FIRST(LELTI(Child, SIGNPF))) < n_p1;

        Child = RECOMPUTE(Child,Q,F,f,P,A);
    }

    /* if we need to split, do it now */
    if (needs_splitting) Children = SPLIT(L, P1);

    /* return. */
    return D;
}

Word SPLIT(Word L, Word Ps)
{
    Word P, p, x, s, C, k;

    // TODO whyyy???
    if (LENGTH(L) < 2) goto Return;

    // initialiseng variables. we assume all cells in the list have the same lever, we also assume that the new proj
    // factors with missing signs are at the end of the list, and that all cells are missing the same number of proj
    // factors.
    C = FIRST(L);
    k = LELTI(C, LEVEL);
    Ps = ADVN(Ps, LENGTH(FIRST(LELTI(C, SIGNPF))));

Step1: /* Iterate through the new polynomials */
    while (Ps != NIL) {
        ADV(Ps, &P, &Ps);
        p = LELTI(P, PO_POLY);
        x = LINEARROOT(p);

        if (x == NIL) {
            SWRITE("ERROR: unable to find a root of new polynomial.\n");

            return L;
        }

        // convert to algebraic extension representation, to be consistent with sample points
        x = AFFRN(x);
        Word NewCells = NIL; // list of lists of new cells.
        // TODO it's linked lists. can we do it more efficiently with ADV rather than LELTI?
        for (int i = 1; i <= LENGTH(L); i++) {
            s = INTERVALSIGN(L, i, x);

            // current cell is not sign invarient - split according to point x
            if (s == 0) {
                NewCells = SPLITCELL(
                    LELTI(L, i),
                    x,
                    SAMPLEK(LELTI(L, MAX(i-1, 1)), k),
                    SAMPLEK(LELTI(L, MIN(i+1, LENGTH(L))), k)
                );
            } else { // it's sign invarient already, simply add the new sign
                APPENDSIGNPF(LELTI(L, i), s);
            }
        }

        // add new cells to list, replacing any cells that were split by a new polynomial
        L = INSERTCELLS(L, NewCells);
    }

Return:
    return L;
}

Word ADVN(Word L, Word n)
{
    Word v;

    for (int i = 0; i < n && L != NIL; i++) {
        ADV(L, &v, &L);
    }

    return L;
}

Word LINEARROOT(Word p)
{
    Word n = PDEG(p);

    // not linear
    if (n > 1) return NIL; // TODO what to do if it's not linear

    // directly solve P := den X + num = 0
    Word den = PLDCF(p);
    Word num = PCOEFF(p,0);

    // x = 0
    if (num == 0) {
        return 0;
    }

    // div by 0
    if (den == 0) return NIL;

    // directly solve x = -num/den
    return RNRED(INEG(num),den);
}

Word SAMPLEK(Word C, Word k)
{
    Word S, P, M, I;

    S = LELTI(C, SAMPLE);

    // fail if the sample point is not in the standard representation (which it should be by now)
    if (LENGTH(S) != 3) return 0;

    // return [ coordinate k, algebraic extension polynomial M, interval I ]
    FIRST3(S, &M, &I, &P);
    return LIST3(LELTI(P, k), M, I);
}

Word INTERVALSIGN(Word L, Word idx, Word x)
{
    Word C, k, Cl, Cr, Sl, Sr, l, r, Ml, Il, Mr, Ir;

Step1: /* Initialise */
    C = LELTI(L, idx);
    k = LELTI(C, LEVEL);

    // next and previous cells
    Sl = NIL; Sr = NIL;

    // special case for even (section) cell
    // we set cells on either side to the value of C
    if (EVEN(idx)) {
        Sl = SAMPLEK(C, k);
        Sr = Sl;
    } else {
        // standard case - sample point of cell to the left of C
        if (idx > 1) {
            Cl = LELTI(L, idx - 1);

            Sl = SAMPLEK(Cl, k);
        }

        // standard case - sample point of cell to the right of C
        if (idx < LENGTH(L)) {
            Cr = LELTI(L, idx + 1);

            Sr = SAMPLEK(Cr, k);
        }

        // special case for leftmost and rightmost cells
        // set Sl and Sr equal
        if (Sl == NIL) Sl = Sr;
        if (Sr == NIL) Sr = Sl;
    }

    /* endpoint comparison */
    // get point, M, I for each point
    FIRST3(Sl, &l, &Ml, &Il);
    FIRST3(Sr, &r, &Mr, &Ir);

    // comparisons
    Word xl = AFCOMP(Ml, Il, x, l); // x is: left < equal < right of l
    Word xr = AFCOMP(Mr, Ir, x, r); // x is: less < equal < right of r

    // decide whether our x is inside, left of, or right of our interval
    if (xl > 0 && xr < 0) return 0;
    if (xr <= 0) return 1;
    if (xl >= 0) return -1;
    return NIL; // if this happened something went badly wrong
}

void APPENDSIGNPF(Word C, Word sign)
{
    Word S, M, D, l;

    // pulling values we need to update off the cell
    S = FIRST(LELTI(C, SIGNPF));
    l = LENGTH(S); // how many SIGNPF
    M = LELTI(C, MULSUB);
    D = LELTI(C, DEGSUB);

    // append new sign
    CONC(S, LIST1(sign));

    // append new degsub (degsub = [d_1, ..., d_n], degrees of k-level projection factors
    // we know it's a linear polynomial
    // TODO do we?? what if it isn't?
    if (D != NIL) CONC(D, LIST1(1));

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

    // TODO check if all the deep copies were done right

    // deep copy children
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

    // deep copy C, adding new deep copy of children
    return MCELL(
        LELTI(C, LEVEL),
        NewChildren,
        NIL,
        LELTI(C, TRUTH),
        LLCOPY(LELTI(C, SAMPLE)),
        LCOPY(LELTI(C, INDX)),
        LLCOPY(LELTI(C, SIGNPF)),
        LELTI(C, HOWTV),
        LCOPY(LELTI(C, DEGSUB)),
        LCOPY(LELTI(C, MULSUB))
    );
}

void SETSAMPLEK(Word C, Word k, Word m)
{
    Word S = LELTI(LELTI(C, SAMPLE), 3);
    SLELTI(S, k, m);
}

Word SPLITCELL(Word C, Word x, Word Sl, Word Sr)
{
    Word Cl, Cr, k, M, l, r, m;

    k = LELTI(C, LEVEL);
    Cl = CELLCOPY(C);
    Cr = CELLCOPY(C);

    // Update indices - Cl takes C's current index, C and Cr are slipped in after it
    INCINDEX(C, k, 1);
    INCINDEX(Cr, k, 2);

    // set sample points
    // Sl : less than x
    FIRST2(Sl, &l, &M);
    m = MIDPOINT(M, l, x);
    SETSAMPLEK(Cl, k, m);

    // C : equal to x
    SETSAMPLEK(C, k, x);

    // Cr : greater than x
    FIRST2(Sr, &r, &M);
    m = MIDPOINT(M, x, r);
    SETSAMPLEK(Cr, k, m);

    // add signs of new projection factor
    APPENDSIGNPF(Cl, -1);
    APPENDSIGNPF(C, 0);
    APPENDSIGNPF(Cr, 1);

    return LIST3(Cl, C, Cr);
}

Word INSERTCELLS(Word LL, Word NewCells)
{
    Word C, k, i, t, L, LRed, NewRed;

    // initialising variables
    C = FIRST(NewCells);
    NewRed = RED(NewCells);
    k = LELTI(C, LEVEL);
    i = LELTI(LELTI(C, INDX), k);

    L = ADVN(LL, i-1); // L = [C_i, ...rest]
    LRed = RED(L);

    // increment LRed indices
    t = LENGTH(NewCells) - 1;
    INCINDEXL(LRed, k, t);

    // concatenate new cells onto remaining cells, and replace first by first new cell
    SRED(L, CONC(NewRed, LRed));
    SFIRST(L, C);

    return LL;
}

Word MIDPOINT(Word M, Word l, Word r)
{
    Word two = AFFINT(2);
    return AFSUM(l, AFQ(
        M,
        AFDIF(r, l),
        two
    ));
}

