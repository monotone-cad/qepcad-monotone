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

// find the solution of a linear polynomial p
Word LINEARROOT(Word p);

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
            printf("Missing proj factors\n");
            Children = SPLIT(L, P1);

            // goto Step3; // go back to the beginning to ensure all cells are OK
            goto Return;
        }

        Child = RECOMPUTE(Child,Q,F,f,P,A);
    }

Return: /* return. */
    return D;
}

Word SPLIT(Word LL, Word Ps)
{
    Word L, P, p, x, s, C, k;

    if (LENGTH(LL) < 2) goto Return;

    C = FIRST(LL);
    k = LELTI(C, LEVEL);
    Ps = ADVN(Ps, LENGTH(LELTI(LELTI(C, SIGNPF), k)));

Step1: /* Iterate through the new polynomials */
    if (Ps == NIL) goto Return;
    ADV(Ps, &P, &Ps);
    p = LELTI(P, PO_POLY);
    x = LINEARROOT(p);

    if (x == NIL) {
        printf("ERROR: root is NIL!\n");

        return L;
    }

    L = LL;

    for (int i = 1; i <= LENGTH(L); i++) {
        s = INTERVALSIGN(L, i, x);
        printf("%d ", s);
    }
    printf("\n");

    goto Step1; // next polynomial in the list Ps
    // get new polynomials via advN
    // for each cell C in L and each new polynomial
    //     check if any roots of P in C.
    //     if no: get sign via evaliating C sample point
    //     if yes: we need to split C

Return:
    return LL;
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

    printf("%d/%d\n", num, den);
    Word x = RNLBRN(RNRED(INEG(num),den));
    return x;
}

// convert k-th coordinate of sample point to rational
Word SAMPLETORATIONAL(Word C, Word k)
{
    Word S, I, p;

    S = LELTI(C, SAMPLE);

    if (LENGTH(S) != 3) return 0;

    I = LELTI(S, 3);
    p = LELTI(I, k);
    if (p == 0) return 0;

    // compute B* (b)
    return IUPBREV(SECOND(p), FIRST(p));
}

Word INTERVALSIGN(Word L, Word idx, Word x)
{
    Word C, k, Cl, Cr, l, r;

Step1: /* Initialise */
    C = LELTI(L, idx);
    k = LELTI(C, LEVEL);

Step2: /* special case for even cell */
    // next and previous cells
    l = NIL; r = NIL;

    if (EVEN(idx)) {
        l = SAMPLETORATIONAL(C, k);
        r = l;

        goto Step3;
    }

    if (idx > 1) {
        Cl = LELTI(L, idx - 1);

        l = SAMPLETORATIONAL(Cl, k);
    }

    if (idx < LENGTH(L)) {
        Cr = LELTI(L, idx + 1);
        r = SAMPLETORATIONAL(Cr, k);
    }

    // special case for leftmost and rightmost cells (whele l == NIL and r == NIL resp).
    if (l == NIL) l = r;
    if (r == NIL) r = l;

Step3: /* endpoint comparison */
    // endpoint comparisons
    Word xl = LBRNCOMP(x, l); // x is: left < equal < right of l
    Word xr = LBRNCOMP(x, r); // x is: less < equal < right of r

    if (xl > 0 && xr < 0) return 0;
    if (xr <= 0) return 1;
    if (xl >= 0) return -1;
    return NIL; // if this happened something went badly wrong
}
