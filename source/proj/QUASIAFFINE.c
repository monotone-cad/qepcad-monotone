/*======================================================================
                    A -< QUASIAFFINE(A, J, r; P, J)

Adding first derivatives of projections onto coordinate axes, to result in 2d cells which are either increasing,
decreasing or constant along coordinate axes.

\Input
  \parm{AA} set of input polynomials (I_1,...,i_r), each A_i is a list of polynomials in Z[x_1,...,x_i]
  \parm{r} number of variables

Output
  \parm{A}: modified projection factors

SideEffect
  \parm{AA} is modified.

======================================================================*/
#include "qepcad.h"

// return the index immediately after I in {1,...,m1} times ... times {1,...,ik} with respect to lex order.
// if no greater index exists, returns NIL
Word LexNext(Word I, Word n, Word *j_)
{
    Word m, I1, J1;
    ADV(I, &m, &I1);

    // increment m, if we can
    if (m < n) {
        return COMP(m + 1, I1);
    }

    // maximal index
    if (I1 == NIL) {
        return NIL;
    }

    // try to roll over
    J1 = LexNext(I1, n, j_);

    if (J1 == NIL) { // fail, maximal index
        return NIL;
    }

    // roll over
    ++(*j_);
    return COMP(0, J1);
}

// return list 0^k
Word ZEROS(Word k)
{
    Word L = NIL;
    while (k > 0) {
        L = COMP(0, L);
        --k;
    }

    return L;
}

// vector component-wise max
Word VCMAX(Word V)
{
    Word m, a, V1;
    ADV(V, &a, &V1);

    // base case
    if (V1 == NIL) return a;

    // recursive case
    return MAX(a, VCMAX(V1));
}

Word LPROD(Word L)
{
    Word a, L1;
    ADV(L, &a, &L1);

    // base
    if (L1 == NIL) return a;

    return a * LPROD(L1);
}

// hash index
// Word

// generate derivatives required for smooth stratification
Word STRAT(Word r, Word Fs, Word Is, Word Hs)
{
    Word j, k, k1, l, r1, j1, i1, I, P;
    k = LENGTH(Fs), k1 = k - 1, l = LENGTH(Is), r1 = r - l, j = k1;

    int Ms[k]; Word* Ds[k];
    // initialise
    while (Fs != NIL) {
        ADV(Fs, &P, &Fs);

        // degree of P_j
        Ms[j] = VCMAX(PDEGV(r, P));

        // Ds: (..., (j,0,...,0), P_j, ... )
        // TODO is this bad for using loads of memory?
        Ds[j] = (Word*) malloc(sizeof(Word) * IEXP(Ms[j], r1));

        --j;
    }

    // generate derivatives
    j1 = 2, i1 = j1; I = COMP2(0, 1, ZEROS(r1 - 1)); // initial index of derivative
    while (I != NIL) {
        j = FIRST(I);
        printf("j1 = %d ", j1); LWRITE(I); SWRITE("\n");

        int inc = 0;
        I = LexNext(I, Ms[j], &inc);

        // update j1, i1
        if (inc > 0) j1 = inc;
        if (i1 < j1) i1 = j1;
    }

    // free the Dps array.
    for (int i = 0; i < k; ++i) free(Ds[i]);

    return NIL;
}

// list of all partials, sufficient for smooth stratification
Word PARTIALS(Word r, Word L)
{
    return STRAT(r, L, NIL, NIL);
}

void QepcadCls::QUASIAFFINE(Word A, Word r, Word *A_)
{
    Word AA, A1, A11, P, k, r1, L;

Step1: /* decide based on dimension */
    if (r <= 1) {
        // nothing to do

        goto Return;
    }

Step3: /* dim >= 2: all partials of input polynomials */
    /* levels. */
    AA = A, L = NIL, k = 0, r1 = r;
    while (AA != NIL) {
        ADV(AA, &A1, &AA);
        ++k; --r1;

        /* polynomials. */
        while (A1 != NIL) {
            ADV(A1, &A11, &A1);
            P = LELTI(A11, PO_POLY);

            L = COMP(PADDVS(LELTI(A11, PO_POLY), r1), L);
        } /* END polynomials. */
    } /* END level. */

    PARTIALS(r, L);
    // factorise and append -- same function as used on input formula
    //ADDPOLS(IPLFAC(k, L),k,LFS("D"), &A);

Return: /* prepare for return */
    // put proj fac in correct order
    *A_ = A;
}

