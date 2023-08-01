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
    return INV(PDEGV(r, P));
}

// recover index from count
// TODO this is only for debugging, can be deleted when done.
// count = (c1 + c2 * m1 + ... + cn * m(n-1))
Word IndexHelper(Word m, Word I, Word* count_)
{
    Word a, I1;
    ADV(I, &a, &I1);

    Word J = NIL;
    if (I1 != NIL) {
        J = IndexHelper(a + 1 * m, I1, count_);
    }

    Word count = *count_;
    Word j = count / m;
    *count_ = count % m;

    return COMP(j, J);
}

inline Word INDEX(Word count, Word M)
{
    if (count < 0) return NIL;

    Word J = IndexHelper(1, M, &count);

    return J;
}

// allocates room for, and sets up, a list of derivatives
Word* Allocate(Word P, Word D)
{
    // polynomials
    // each F_1 is of the form (..., P, deg(P), )
    Word P1 = LIST2(P, D);

    // calculate max number of derivatives possible in list
    Word len = 1, d = 0;
    while (D != NIL) {
        ADV(D, &d, &D);
        len = len * (d + 1);
    }

    len += 2; // includes P itself and a null at the end

    Word* F1s = (Word*) malloc(len * sizeof(Word));
    if (F1s == NULL) {
        FAIL("QUASIAFFINE", "malloc F1s failed.");
    }

    // first element is P1, rest are NIL (empty lists)
    F1s[0] = P1;

    // rest are null
    Word j = 1;
    while (j < len) {
        F1s[j] = NIL;

        ++j;
    }

    return F1s;
}

// list sum
inline Word LSUM(Word L)
{
    Word sum, a;

    sum = 0;
    while (L != NIL) {
        ADV(L, &a, &L);

        sum += a;
    }

    return sum;
}

// recursive smooth stratification of polynomials Fs
// np: number of polynomials
// r: number of variables
// Fs: list of polynomials and degrees (..., P_i, deg(P_i), ...)
// Is: indices (i1,...,ik)
// Hs: list of polynomials h_i1,...,h_ik
// Minor: (LENGTH(I1) - 1) * (LENGTH(I1) - 1) - matrix |a_ij| = partial h_i / partial x_j for h in Hs and i in Is
// return list of all differentials computed by alrogithm
Word STRAT(Word np, Word r, Word** Fs, Word Is, Word Hs, Word Minor)
{
    printf("Round %d\n", LENGTH(Is)); // TODO debugging

    // end of recursion
    Word i0 = FIRST(Is); // number of variables considered so far
    if (np == 0 || i0 == r) {
        return NIL;
    }

    // TODO debugging
    SWRITE("STRAT: Is = "); LWRITE(Is); SWRITE("\n");
    SWRITE("       Hs = "); LWRITE(Hs); SWRITE("\n");

    // set up return value
    Word Gs1 = NIL; // list of all differentials computed in this round, to return

    // set up working array
    Word g_len = np * 2; // length of memory allocated for Gs
    Word g_count = 0; // how many differentials computed so far, index in Gs
    Word** Gs = (Word**) malloc(g_len * sizeof(Word*)); // list of all differentials computed in this round, working set

    // metadata
    Word* Ps[np]; // "chase list", first element is h1
    Word Degrees[np]; // degree of Fs[0][0], gives max index
    Word Ms[np]; // number of steps before maximum differentiation variable (i1) should be incremented
    Word Dvs[np]; // list of current differentiation variable (i1)
    Word Chase_index[np]; // chase indices TODO debugging

    Word init_v = i0 + 1; // initial differentiation variable (i1)
    Word p_index = 0; // initial polynomial index, ranges over 0 <= p_index < np

    // initialise metadata for each input polynomial
    while (p_index < np) {
        Word* F1s = Fs[p_index];
        Word D = REDI(SECOND(F1s[0]), i0);

        Ps[p_index] = F1s;
        Degrees[p_index] = D;
        Dvs[p_index] = init_v;
        Ms[p_index] = FIRST(D);

        Chase_index[p_index] = 0; // TODO debugging

        // increment
        ++p_index;
    }

    // main loop, consider each index (j, m_{i0 + 1}, ..., m_r) in lex order
    Word n_finished = 0; // each polynomial has a different max index, keep track of how many indices are maxed out
    Word count = 0; // number of derivatives computed so far, scalar value of (m_{i0 + 1}, ..., m_r)

    while (n_finished < np) { // stop once differential index for every polynomial is maxed
        // reached the end of polynomial list, cycle back to beginning and consider next differential index I
        if (p_index == np) {
            p_index = 0;
            n_finished = 0;

            ++count;
            printf("increment count %d. cycle polynomial list\n", count);
        }

        // compute s_k = partial_{(h_1,...,h_{k-1}),(i_1,...,i_{k-1}),v} h_k
        Word v = Dvs[p_index]; // differentiation variable
        printf("v = %d\n", v);
        Word ch_index = Chase_index[p_index]; // TODO debugging

        // maxed out index, no more differentials possible -- skip
        if (v > r) {
            printf("skipping %d\n", p_index);
            ++n_finished;
            ++p_index;

            continue;
        }

        // get h_k and its degree
        Word D = Degrees[p_index];
        Word P = FIRST(*Ps[p_index]);

        // TODO partial differential operator with h1,...,hk
        // throw out constant polynomials
        Word Q = IPDER(r, P, v);
        Word Qdeg = DEG(r,Q);

        if (LSUM(Qdeg) == 0) { // s_k is constant, may as well be zero
            Q = 0;
        }

        // TODO debugging
        IWRITE(count); SWRITE(": compute derivative, variable "); IWRITE(v); SWRITE(" polynomial "); IWRITE(p_index);
        SWRITE(", degree "); LWRITE(Degrees[p_index]); SWRITE("\n  ");
        LWRITE(INDEX(ch_index, D)); SWRITE(" ");
        P == 0 ? SWRITE("0") : LWRITE(P); SWRITE("\n  ");
        LWRITE(INDEX(count, D)); SWRITE(" "); Q == 0 ? SWRITE("0") : LWRITE(Q); SWRITE("\n\n");

        if (Q != 0) {
            printf("recurse...\n");
            Word Gs2 = STRAT(g_count, r, Gs, COMP(v, Is), COMP(P, Hs), Minor);
            Gs1 = CONC(Gs1, Gs2);

            // append new derivative to working list Gs and return list
            printf("add derivative to Gs %d\n", g_count);
            if (g_count >= g_len) { // enlarge list
                printf("realloc() %d\n", g_len);
                g_len += np;
                Gs = (Word**) realloc(Gs, g_len * sizeof(Word*));
            }

            Gs[g_count] = Allocate(Q, Qdeg);
            ++g_count;
        }

        // append derivative to Fs, preserving zeroes
        Fs[p_index][count] = LIST2(Q, Qdeg);

        // update variable v and chaser list
        if (count < Ms[p_index]) { // same variable, higher order -- get next element in "chase" list
            Ps[p_index] = Ps[p_index] + 1;
            Chase_index[p_index] = Chase_index[p_index] + 1; // TODO debugging
        } else { // rollover - increment differentiation variable
            Dvs[p_index] = v + 1; // next variable
            Ps[p_index] = Fs[p_index]; // back to beginning of chase list
            Ms[p_index] = Ms[p_index] * LELTI(Degrees[p_index], v - i0); // how many derivatives until we need to roll over again

            Chase_index[p_index] = 0; // TODO debugging
        }

        // consider next polynomial
        ++p_index;
    }

    // construct list Gs1 of all functions in Gs
    Word i = g_count;
    while (i > 0) { // comp backwards
        --i;
        Gs1 = COMP(FIRST(Gs[i][0]), Gs1);
    }

    // free Gs
    for (int i = 0; i < g_count; ++i) {
        free(Gs[i]);
    }

    free(Gs);

    // returns list of derivatives
    return Gs1;
}

// list of all partials, sufficient for smooth stratification
Word PARTIALS(Word r, Word L)
{
    Word D, P, P1, k, j;

    k = LENGTH(L);
    Word** Fs = (Word**) malloc(k * sizeof(Word*));
    if (Fs == NULL) {
        FAIL("QUASIAFFINE", "malloc Fs failed.");
    }

    j = 0;
    while (L != NIL) {
        ADV(L, &P, &L);

        D = DEG(r, P);
        Fs[j] = Allocate(P, D);

        ++j;
    }

    // initial i0 = FIRST(I1) = 0. h0 = FIRST(Hs) = 0, Minor is the empty matrix
    STRAT(k, r, Fs, LIST1(0), LIST1(0), NIL);

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

