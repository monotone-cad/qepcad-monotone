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

// convenience function for multi-degree of polynomial, returns (d_r,...,d_1) where d_i is degree of P in variable i
inline Word DEG(Word r, Word P)
{
    return PDEGV(r, P);
}

// recover index from count
// count = (c1 + c2 * m1 + ... + cn * m(n-1))
Word IndexHelper(Word *m, Word* count, Word M)
{
    Word a, M1;
    ADV(M, &a, &M1);

    // update m = (m1 * ... * mk)
    // special case for a = 0
    if (a > 0)
        *m = (a * *m);

    // base case: last element
    if (M1 == NIL) {
        return NIL;
    }

    Word I = IndexHelper(m, count, M1);

    // special case for a == 0
    if (a == 0)
        return COMP(0, I);

    Word j = *count / *m;

    *count = *count % *m;
    return COMP(j, I);
}

inline Word INDEX(Word count, Word M)
{
    Word J, m = 1;
    J = IndexHelper(&m, &count, M);

    return COMP(count, J);
}

// recursive smooth stratification of polynomials Fs
// np: number of polynomials
// r: number of variables
// Fs: list of polynomials and degrees (..., deg(P), Ps, ...)
// return: TODO NIL
Word STRAT(Word np, Word r, Word** Fs)
{
    Word* Ps[np]; // pointer to next polynomial to take derivative of, index (m_1,...,m_i1,0,...,0)
    Word* Qs[np]; // chaser list of polynomials (m_1,...,m_i1 - 1,0,...,0) // TODO don't need, i used Ps
    Word Ds[np]; // degrees of input polynomials, to calculate differentiation index
    Word Ms[np]; // how many steps until we update v
    Word v_index[np];

    Word p_index = 0; // initial index, variable to be differentiated in h1 to obtain s1

    // initialise
    printf("initialising...\n");
    while (p_index < np) {
        // for each input polynomial P_i
        Word* F1s = Fs[p_index];
        Word D = SECOND(F1s[0]);

        // store backup and chaser list
        Ps[p_index] = F1s;
        Qs[p_index] = F1s;

        // save initial degree
        Ds[p_index] = D;

        // set Ms: initially d1
        Ms[p_index] = FIRST(D);

        // initial v_index
        v_index[p_index] = 1;

        // increment
        ++p_index;
    }

    // take derivatives
    Word nv = LENGTH(SECOND(Fs[0][0])); // length of index == number of variables
    Word n_rollover = 0; // when we rolled over nv * k times, we're done.
    Word count = 0; // counts how many derivatives we have computed so far
    while (n_rollover < nv * np) {
        // loop polynomial list
        if (p_index == np) {
            p_index = 0;
            ++count;
        }

        Word v = v_index[p_index];

        // no more derivatives
        if (v > nv) {
            ++p_index; // move to next polynomial

            continue;
        }

        // get data for polynomial P_i
        Word D = Ds[p_index];
        Word P = FIRST(*Ps[p_index]);

        // s1 = partial h_1 / x_i
        Word Q = IPDER(r, P, v);

        // TODO recursion
        // we know P == h1, i1 == v

        // append new derivative Q
        IWRITE(p_index); SWRITE(", "); LWRITE(INDEX(count, D)); SWRITE(" -- "); Q == 0 ? SWRITE("00") : LWRITE(Q); SWRITE("\n");
        Fs[p_index][count] = LIST2(Q, DEG(r, Q));

        // update variable v and chase list
        if (count < Ms[p_index]) {
            // same variable, higher order
            Ps[p_index] = Ps[p_index] + 1; // move the "chase" pointer along by one
        } else {
            // next variable, order 0 (index "roll-over"
            v_index[p_index] = v + 1;
            Ps[p_index] = Fs[p_index]; // reset chase ist
            Ms[p_index] = Ms[p_index] * LELTI(Ds[p_index], v);
            printf("update i1 = %d, m = %d\n", v_index[p_index], Ms[p_index]);

            // keep track of how many times this was done. if we rollower n+1 times, that polynomial is note.
            ++n_rollover;
        }

        // increment
        ++p_index;
    }
    return NIL;
}

inline Word LPROD(Word L)
{
    Word a, m;

    m = 1;
    while (L != NIL) {
        ADV(L, &a, &L);
        m *= a;
    }

    return m;
}

// list of all partials, sufficient for smooth stratification
Word PARTIALS(Word r, Word L)
{
    Word D, P, P1, k, j;

    k = LENGTH(L);
    Word** Fs = (Word**) malloc(k * sizeof(Word*));

    j = 0;
    while (L != NIL) {
        ADV(L, &P, &L);

        // polynomials
        // each F_1 is of the form (..., P, deg(P), )
        D = DEG(r, P);
        P1 = LIST2(P, D);
        Fs[j] = (Word*) malloc(LPROD(D) * sizeof(Word));
        Fs[j][0] = P1;

        ++j;
    }

    STRAT(k, r, Fs);

    // free 2d Fs
    for (j = 0; j < k; j++) free(Fs[j]);
    free(Fs);

    return NIL;
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

