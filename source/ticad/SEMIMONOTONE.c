/*======================================================================
                    F -< SEMIMONOTONE(FF)

Adds extra polynomials to projection factor set to ensure semi-monotone cells will be produced

\Input
  \parm{A} is a list a_1, ..., a_r where a_ is a list of i-level projection factors
  \parm{AA} is a list a_1, ..., a_r where a_ is a list of i-level input polynomials
  \parm{D} is a CAD
  \parm{F} is the normalised input QFF
  \parm{r} is the space in which D lives

  Output
  \parm{D} new cad with semi-monotone cells

======================================================================*/
#include "qepcad.h"

// list of levels at which the cell is one-dimensional
Word LEVELIDX(Word C);

// factorise and null parents of a lest of level k polynomials
Word FAC(Word L, Word k);

// get the parent (in D) of cell C, at level k
Word PARENT(Word C, Word k, Word D);

// split cell C based on level k proj factors
void SPLITCELL(Word A, Word C, Word k, Word D);

// increment level k index in L (list of cells) by *by*
Word INCINDEX(Word L, Word k, Word by);

// recursive (deep, incl. children) cell duplication
Word DEEPCELLDUP(Word C);

// does ADV n times and returns what is left of the list. if n > LENGTH(L), then NIL
Word ADVN(Word L, Word n);

// append a new value to C's signpf at level k
void SETSIGNPF(Word C, Word k, Word i, Word sign);


Word QepcadCls::SEMIMONOTONE(Word A, Word AA, Word DD, Word F, Word r)
{
    Word D, Ct, Cf, C, I, i ,j, A1, P, p, P1, p1, L, L1, Cs, k, S, Cc;

Step1: /* Initialise */
    D = DD;
    Cs = NIL; /* list of cells that need re-lifting */

Step2: /* calculate: looping through true cells */
    LISTOFCWTV(D, &Ct, &Cf);
    while (Ct != NIL) {
        ADV(Ct, &C, &Ct);
        I = LEVELIDX(C);
        if (LENGTH(I) > 2) {
            SWRITE("*** ERROR doesn't handle cells of dimension > 2");

            return D;
        } else if (LENGTH(I) < 2) {
            printf("# cell dim < 2 - continue\n");

            continue;
        }

        /* We have a 2-dimensional cell */
        FIRST2(I, &i, &j);

        printf("# two dimensional index: %d %d\n", i, j);
        /* adding partials of level j polynomials with respect to variable i */
        L = NIL;
        A1 = LELTI(A, j);
        while (A1 != NIL) {
            ADV(A1, &P, &A1);
            p = LELTI(P, PO_POLY);
            p1 = IPDER(j, p, i);

            if (IPCONST(j, p1)) continue;

            P1 = MPOLY(p1, NIL, NIL, PO_OTHER, PO_KEEP);
            L = COMP(P1, L);
        }

        printf("# %d ", LENGTH(LELTI(A, i)));
        L = FAC(L, j);
        A = APPEND(A, j, L);
        printf("%d\n", LENGTH(LELTI(A, i)));

        /* Append cell that needs re-lifting */
        Cs = COMP(PARENT(C, i, D), Cs);
    }

Step3: /* splitting cells */
    printf("# %d cells to split\n", LENGTH(Cs));fflush(0);
    while (Cs != NIL) {
        ADV(Cs, &C, &Cs);
        k = LELTI(C, LEVEL);
        A1 = LELTI(A, k);
        SPLITCELL(A1, C, k, D);

        // INSERT CELL
        // - get siblings and insert index
        // - iterate to insert index (see LINS)
        // - insert into list and concat the rest
        // - iterate over the rest (including all children) and increment level i index
        // update index function, C, k, n) recursively set level k index on all cells to i
    }

Return: /* returning */
    return D;
}

Word LEVELIDX(Word C)
{
    Word I, k, j, L;

    L = NIL;
    k = 0;
    I = LELTI(C, INDX);
    while (I != NIL) {
        ADV(I, &j, &I);
        k++;

        if (ODD(j)) {
            L = COMP(k, L);
        }
    }

    return INV(L);
}

Word FAC(Word L, Word k)
{
    Word LL, P;

    L = IPLFAC(k, L);

    LL = L;
    while (L != NIL) {
        ADV(L, &P, &L);
        SLELTI(P, PO_PARENT, NIL);
    }

    return LL;
}

Word PARENT(Word C, Word k, Word D)
{
    Word I, i, j, Cp, Cc;

    Cp = D;
    I = LELTI(C, INDX);

    for (i = 1; i <= k; i++) {
        j = LELTI(I, i);
        Cc = LELTI(Cp, CHILD);
        Cp = LELTI(Cc, j);
    }

    return Cp;
}


void SPLITCELL(Word A, Word C, Word k, Word D)
{
    Word L, Siblings, Rest, j, i, Cs, Cnew, Parent, NewPolynomials, newIndx, idx;

    Parent = PARENT(C, k-1, D);
    Siblings = LELTI(Parent, CHILD);
    L = Siblings;
    j = LELTI(LELTI(C, INDX), k);

    /* shift along siblings until we reach cell C */
    /* Siblings = [C, ...Rest] */
    Siblings = ADVN(Siblings, j-1);

    /* figure out how many cells are needed */
    // TODO something something CSA something something
    NewPolynomials = ADVN(A, LENGTH(LELTI(LELTI(C, SIGNPF), k)));
    if (LENGTH(NewPolynomials) > 1) {
        SWRITE("TODO: I can only handle one polynomial - need a proper csa algorithm please.\na");

        return;
    }

    idx = LENGTH(A) - LENGTH(NewPolynomials) + 1;
    printf("# new idx = %d\n", idx);
    /* there will be 2 * LENGTH(NewPolynomials) + 1 cells, re-use cell C as first */

    /* we know in advance that we are adding 2 cells, so we increment future indices by 2 */
    Rest = RED(Siblings);
    Rest = INCINDEX(Rest, k, 2);

    /* duplicate C and insert cell for ZERO and POSITIVE */
    Cnew = DEEPCELLDUP(C);
    /* update index */
    SLELTI(LELTI(Cnew, INDX), k, j+2);
    SLELTI(Cnew, CHILD, INCINDEX(LELTI(Cnew, CHILD), k, 2));
    SETSIGNPF(Cnew, k, idx, POSITIVE);

    Rest = COMP(Cnew, Rest);
    // TODO repeaty code is bad!!! fix me please
    Cnew = DEEPCELLDUP(C);
    /* update index */
    SLELTI(LELTI(Cnew, INDX), k, j+1);
    SLELTI(Cnew, CHILD, INCINDEX(LELTI(Cnew, CHILD), k, 1));

    // SETSIGNPF(Cnew, k, idx, ZERO);

    Rest = COMP(Cnew, Rest);

    // SETSIGNPF(C, k, idx, NEGATIVE);

    /* finally assign the new Rest (w/ new cells) te Siblings */
    SRED(Siblings, Rest);

    printf("# %d\n", LENGTH(L));

    /* re-assign the siblings */
    SLELTI(Parent, CHILD, L);
}

Word INCINDEX(Word L, Word k, Word by)
{
    Word C, I, Children, j, LL;
    LL = L;

    while (L != NIL) {
        ADV(L, &C, &L);
        I = LELTI(C, INDX);
        Children = LELTI(C, CHILD);
        if (Children != NIL) {
            SLELTI(C, CHILD, INCINDEX(Children, k, by));
        }

        j = LELTI(I, k);
        SLELTI(I, k, j + by);
    }

    return LL;
}

Word DEEPCELLDUP(Word C)
{
    Word Children, Ch, NewChildren;

    // TODO figure out which ones should be LLCOPied
    Children = LELTI(C, CHILD);
    NewChildren = NIL;
    while (Children != NIL) {
        ADV(Children, &Ch, &Children);

        NewChildren = COMP(
            DEEPCELLDUP(Ch),
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

Word ADVN(Word L, Word n)
{
    Word v;

    for (int i = 0; i < n; i++) {
        if (L == NIL) break;
        ADV(L, &v, &L);
        printf("~ advancing\n");
    }

    printf("# %d (in ADVN)\n", LENGTH(L));
    return L;
}

void SETSIGNPF(Word C, Word k, Word i, Word sign)
{
    Word S, Children, Ch, level;
    level = LELTI(C, LEVEL);
    S = LELTI(LELTI(C, SIGNPF), level - k + 1);

    if (i > LENGTH(S)) {
        S + CONC(S, LIST1(sign));
    } else {
        SLELTI(S, i, sign);
    }

    SLELTI(LELTI(C, SIGNPF), k, S);

    // recursively append sign to children
    Children = LELTI(C, CHILD);
    while (Children != NIL) {
        ADV(Children, &Ch, &Children);

        SETSIGNPF(Ch, k, i, sign);
    }
}

