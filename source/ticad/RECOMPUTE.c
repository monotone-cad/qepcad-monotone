/*======================================================================
                      D <- RECOMPUTE(C,Q,F,f,P,A)

Recompute cad taking into account new projection factors. Takes as input an existing CAD along with a projection factor
structure which may contain projection factors not included in the CAD. Cells are split, and their data (signature and
sample points) are updated according to these new polynomials.

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
Word INTERVALSIGN(Word Cl, Word C, Word Cr, Word x);

// -1 if < 0, 1 if > 0. same as intervalsign, but for section cells
Word CELLSIGN(Word C, Word x);

// get the k-th coordinate of sample point of C
Word SAMPLEK(Word C, Word k);

// solution of a polynomial as algebraic number
// if 0 then 0,
// if rational then algebraic representation of rational number
// otherwise (b, M, I) where M = P, I is the isolating interval, and b = f(x) = x
Word ROOTS(Word p);

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

// "mid" point of two pseudo sample points.
// if either l or r is not rational, then we take the inside end of isolating interval,
// otherwise, both are rational, and we take the fraction (a-b)/2
// returns a pseudo sample point, which is rational
Word MIDPOINT(Word l, Word Ml, Word Il, Word r, Word Mr, Word Ir);

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
    Word P, p, xs, xs1, x, s, C, Next, Prev, Lp, Lp1, Lc, Lc1, Ln, Ln1, NewCells, k;

    // initialising variables. we assume all cells in the list have the same level, we also assume that the new proj
    // factors with missing signs are at the end of the list, and that all cells have the same number of proj factors so
    // far
    C = FIRST(L);
    k = LELTI(C, LEVEL);
    // the missing proj factors
    Ps = ADVN(Ps, LENGTH(FIRST(LELTI(C, SIGNPF))));

    // construct list pointers for looping over cells later
    // walk along the list L of cells, keeping track of Current, Next and Previous.
    // this is done by having 3 list pointers, for Prev, Current and Next, and advancing them in sync
    Lc = LCOPY(L), Lp = Lc, Ln = INV(Lc);
    Lp = COMP(FIRST(Lp), Lp); // = [a1, a1, ..., an], Note we won't need an
    Ln = RED(INV(COMP(FIRST(Ln), Ln))); // = [a2, ..., an, an]

    /* find roots of the new polynomials represented as algebraic numbers */
    xs = NIL;

    while (Ps != NIL) {
        ADV(Ps, &P, &Ps);
        p = LELTI(P, PO_POLY);
        xs1 = ROOTS(p);

        if (xs1 == NIL) {
            SWRITE("ERROR: unable to find a root of new polynomial.\n");

            return L;
        }

        xs = CONC(xs, xs1);
    }

    printf("number of roots: %d\n", LENGTH(xs));
    /* iterate through the new roots, add new signatures and split if needed */
    while (xs != NIL) {
        ADV(xs, &x, &xs);
        SWRITE("### splitting point: ");
        Word Ix = THIRD(x);
        SWRITE("I=("); RNWRITE(FIRST(Ix)); SWRITE(", "); RNWRITE(SECOND(Ix)); SWRITE(") "); LWRITE(SECOND(x));
        SWRITE(" "); RNWRITE(FIRST(FIRST(x))); SWRITE("\n");

        NewCells = NIL, Lp1 = Lp, Lc1 = Lc, Ln1 = Ln, Prev = NIL, C = NIL, Next = NIL;
        bool section = true; // currently looking at a section cell? (initially we are looking at -\infty)
        // LENGTH(Ln) == LENGTH(L) - Ln is the only list having the correct length
        // since Lc has concatenated last element (from Ln) and Lp has concatenated first and last elements
        // this is the price of pass by reference I guess :(
        while (Ln1 != NIL) {
            // advance the 3 list pointers in step
            ADV(Lp1, &Prev, &Lp1);
            ADV(Lc1, &C, &Lc1);
            ADV(Ln1, &Next, &Ln1);
            section = !section;

            // section cell, append new sign
            if (section) {
                APPENDSIGNPF(C, CELLSIGN(C, x));

                continue;
            }

            // sector cell. may contain splitting point, otherwise just add the sign
            s = INTERVALSIGN(Prev, C, Next, x);
            if (s == 0) {
                // we have a change of sign within C
                printf("splitting cell "); LWRITE(LELTI(C, INDX)); printf("\n");
                NewCells = SPLITCELL(
                    C,
                    x,
                    SAMPLEK(Prev, k),
                    SAMPLEK(Next, k)
                );

                // add new cells to list, replacing any cells that were split by a new polynomial
                // TODO hopefully editing L won't screw things up
                // maybe we can append newCells to the end and add all at once??
                L = INSERTCELLS(L, NewCells);
            } else {
                APPENDSIGNPF(C, s);
            }
        }
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
    Word den = PLDCF(p);
    Word num = PCOEFF(p,0);

    // x = 0
    if (num == 0) {
        return 0;
    }

    // div by 0
    if (den == 0) return NIL;

    // directly solve x = -num/den
    return AFFRN(
        RNRED(INEG(num),den)
    );
}

Word ROOTS(Word p)
{
    Word n = PDEG(p);

    // don't think this can happen, but just in case
    if (n == 0) {
        return NIL;
    }

    // P is linear: directly solve P := den X + num = 0
    // return in pseudo-sample point representation
    if (n == 1) {
        return LIST1(
            LIST3(
                LINEARROOT(p),
                PMON(1,1),
                LIST2(0,0)
            )
        );
    }

    // not linear: find isolating intervals for the roots and return in rational number form
    // we will use P as the minimal polynomial. beta = P(alpha)
    Word B = LIST2(LIST2(1,1), PMON(1,1)); // univariate polynomial: (1/1)(1*x^1)
    Word Is = IPRIM(p);
    Word M = NIL, I = NIL, J = NIL;
    while (Is != NIL) {
        ADV(Is, &I, &Is);
        // refine the isolating interval I // TODO what should p be?
        J = LBRIBRI(ANR(1, p, BRILBRI(I))); // TODO is there a way of avoiding all the logarithmic conversion and back?

        M = COMP(
            LIST3(B, p, J),
            M
        );
    }

    return INV(M);
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

inline Word RATIONALPART(Word x)
{
    return x == 0 ? 0 : FIRST(x);
}

Word CONVERTPOLYNOMIAL(Word P)
{
    Word e, a;
    Word P1 = NIL;
    while (P != NIL) {
        ADV2(P, &e, &a, &P);
        P1 = COMP2(AFFINT(a), e, P1);
    };

    return INV(P1);
}

// algebraic number (not rational) comparison
// since r is splitting point, we may assume that beta == alpha.
Word ANCOMPARE(Word Mx, Word Ix, Word y, Word My, Word Iy)
{
    // TODO something with the isolating interval for efficiency, before messing about with algebraic numbers

    // compute the common extension Gamma, and comvert y into that
    Word G, J, t, u, a, b, junk, L, a1, a2;
    // convert My into Q(alpha)[y] representation.
    My = CONVERTPOLYNOMIAL(My);
    // printf("Mx = "); LWRITE(Mx); printf(" ");
    // printf("My = "); LWRITE(My); printf("\n");

    SIMPLEQE(Mx, Ix, My, Iy, &G, &t, &u, &J, &a, &b, &junk, &junk);
    MODCRDB(LIST1(y), G, a, b, &L);
    FIRST2(L, &a1, &a2);
    // printf("t = %d, a = ", t); (a == 0 ? IWRITE(0) : LWRITE(a)); printf(" ");
    // printf("u = %d; b = ", u); (b == 0 ? IWRITE(0) : LWRITE(b)); printf("\n");
    // printf("GAMMA = "); LWRITE(G); printf(" ");
    // printf("J = "); LWRITE(J); printf("\n");
    // printf("AFCOMP of ");
    // (a1 == 0 ? IWRITE(0) : LWRITE(a1));
    // printf(" ");
    // (a2 == 0 ? IWRITE(0) : LWRITE(a2));
    // printf("\n");

    // -1 because x and y are swapped // TODO maybe swap it back so it is less confusing?
    return -AFCOMP(G, J, a1, a2);
}

// helper function for comparing two pseudo-sample points. R is our "splitting point" and S is the sample point of the
// cell we are examining
Word COMPARE(Word R, Word S)
{
    const Word Px = PMON(1,1); // P(x) = x - this is the convention used for an algebraic representation of a rational number
    Word x, Mx, Ix, y, My, Iy;

    // get coordinate, minimal polynomial, and isolating interval
    FIRST3(R, &x, &Mx, &Ix);
    FIRST3(S, &y, &My, &Iy);

    // we determine the type of point (algebraic extension or rational) to figure out which comparison we can do
    bool x_rat = PDEG(Mx) == 1;
    bool y_rat = PDEG(My) == 1;

    // both are rational: easy! just RNCOMP the rational bit
    if (x_rat && y_rat) {
        return RNCOMP(
            RATIONALPART(y),
            RATIONALPART(x)
        );
    }

    // both algebraic (non rational)
    if (!x_rat && !y_rat) {
        return ANCOMPARE(
            Mx, Ix,
            y, My, Iy
        );
    }

    // one is rational, we can use the minimal polynomial and rational point to figure out the sign
    Word ratPart = RATIONALPART(x_rat ? x : y);
    Word minPol = x_rat ? My : Mx;

    // the sign of M at x
    Word s = IUPBES(minPol, ratPart);
    // may need to flip the sign if M and x are the other way around
    return s * (x_rat ? -1 : 1);
}

Word CELLSIGN(Word C, Word R)
{
    Word k = LELTI(C, LEVEL);
    Word S = SAMPLEK(C, k);

    return COMPARE(R, S);
}

Word INTERVALSIGN(Word Cl, Word C, Word Cr, Word x)
{
    Word k = LELTI(C, LEVEL);
    Word Sl = SAMPLEK(Cl, k);
    Word Sr = SAMPLEK(Cr, k);

    // TODO remove, some day.
    int comp_l = COMPARE(x, Sl);
    int comp_r = COMPARE(x, Sr);

    // TODO is there a less iffy way of doing that?
    if (COMPARE(x, Sl) > 0) { // right of left endpoint, + on cell
        return 1;
    } else if (COMPARE(x, Sr) < 0) { // left of right endpoint, - on cell
        return -1;
    } else { // middle
        return 0;
    }
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

// helper for cell copy
// copies a child cell, handling the structure sharing some pointers with it's parent
Word CELLCOPYCHILD(Word C, Word Parent)
{
    // TODO everywhere where it's LLCOPY, should take list from parent and append first of list from child to be copied.

    Word signpf = LELTI(Parent, SIGNPF);
    signpf = COMP(
        LCOPY(FIRST(LELTI(C, SIGNPF))),
        signpf
    );

    Word NewC = MCELL(
        LELTI(C, LEVEL),
        NIL, // children, set later
        NIL,
        LELTI(C, TRUTH),
        LLCOPY(LELTI(C, SAMPLE)),
        LCOPY(LELTI(C, INDX)),
        signpf,
        LELTI(C, HOWTV),
        LCOPY(LELTI(C, DEGSUB)),
        LCOPY(LELTI(C, MULSUB))
    );

    Word Children = LELTI(C, CHILD);
    if (Children == NIL) return NewC;

    Word NewChildren = NIL, Child = NIL;
    while (Children != NIL) {
        ADV(Children, &Child, &Children);

        NewChildren = COMP(
            CELLCOPYCHILD(Child, NewC),
            NewChildren
        );
    }

    SLELTI(NewC, CHILD, INV(NewChildren));
    return NewC;
}

Word CELLCOPY(Word C)
{
    Word Children, Ch, NewChildren, NewC;

    // TODO check if all the deep copies were done right

    // deep copy of C
    // initially set children to NULL, as they depend on the copy of C and will be updated later
    NewC = MCELL(
        LELTI(C, LEVEL),
        NIL, // children, set later
        NIL,
        LELTI(C, TRUTH),
        LLCOPY(LELTI(C, SAMPLE)),
        LCOPY(LELTI(C, INDX)),
        LLCOPY(LELTI(C, SIGNPF)),
        LELTI(C, HOWTV),
        LCOPY(LELTI(C, DEGSUB)),
        LCOPY(LELTI(C, MULSUB))
    );

    // deep copy children.
    Children = LELTI(C, CHILD);

    if (Children == NIL) { // we're done, no children to copy
        return NewC;
    }

    NewChildren = NIL;
    while (Children != NIL) {
        ADV(Children, &Ch, &Children);

        NewChildren = COMP(
            CELLCOPYCHILD(Ch, NewC),
            NewChildren
        );
    }

    SLELTI(NewC, CHILD, INV(NewChildren));

    return NewC;
}

void SETSAMPLEK(Word C, Word k, Word m)
{
    // set sample point for all of C's children
    Word Children, Child;
    Children = LELTI(C, CHILD);
    while (Children != NIL) {
        ADV(Children, &Child, &Children);
        SETSAMPLEK(Child, k, m);
    }

    // set sample point of cell C
    Word S, Ms, Is, B, My, Iy, y;
    S = LELTI(C, SAMPLE);
    FIRST3(S, &Ms, &Is, &B); // qepcad sample point
    FIRST3(m, &y, &My, &Iy); // pseudo-sample point (yes it's confusing, should probably change it to be consistent)

    // set coordinate k
    SLELTI(B, k, y);

    // set the new algebraic extension

    if (PDEG(My) == 1) return; // we're good, new point is rational

    if (PDEG(Ms) == 1) { // all other points are rational, we can overwrite M and I.
        SLELTI(S, 1, My);
        SLELTI(S, 2, Iy);
    } else { // we will have to find a new algebraic extension encompassing both M and M1.
        // TODO
        printf("TODO: comput a new algebraic extension\n");
    }
}

Word SPLITCELL(Word C, Word x, Word Sl, Word Sr)
{
    Word Cl, Cr, k, M, I, b, m, xx, Mx, Ix;

    k = LELTI(C, LEVEL);
    Cl = CELLCOPY(C);
    Cr = CELLCOPY(C);

    // Update indices - Cl takes C's current index, C and Cr are slipped in after it
    INCINDEX(C, k, 1);
    INCINDEX(Cr, k, 2);

    // set sample points
    FIRST3(x, &xx, &Mx, &Ix);
    // Sl : less than x
    FIRST3(Sl, &b, &M, &I);
    m = MIDPOINT(b, M, I, xx, Mx, Ix);
    SETSAMPLEK(Cl, k, m);

    // C : equal to x
    SETSAMPLEK(C, k, x);

    // Cr : greater than x
    FIRST3(Sr, &b, &M, &I);
    m = MIDPOINT(xx, Mx, Ix, b, M, I);
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

Word MIDPOINT(Word l, Word Ml, Word Il, Word r, Word Mr, Word Ir)
{
    Word Px = PMON(1,1);
    Word m;

    // if endpoint is rational, use that, if not, use the endpoint of isolating interval
    // convert to binary rational representation for later.
    if (PDEG(Ml) > 1) {
        l = SECOND(Il);
    } else {
        l = RATIONALPART(l);
    }

    if (PDEG(Mr) > 1) {
        r = FIRST(Ir);
    } else {
        r = RATIONALPART(r);
    }

    // return midpoint of two rational endpoints (approximate, returning a simple fraction)
    m = LBRNRN(SLBRNIM(RNLBRN(l), RNLBRN(r)));

    return LIST3(
        AFFRN(m),
        Px,
        LIST2(0,0)
    );
}

