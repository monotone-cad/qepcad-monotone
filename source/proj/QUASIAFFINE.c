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
// return: TODO NIL
Word STRAT(Word np, Word r, Word** Fs, Word Is, Word Hs)
{
    printf("Round %d\n", LENGTH(Is));

    // end of recursion
    Word i0 = FIRST(Is); // number of variables already differentiated
    if (i0 == r || np == 0) {
        printf("base of recorsion\n");

        return NIL;
    }

    SWRITE("STRAT: Is = "); LWRITE(Is); SWRITE("\n");
    SWRITE("       Hs = "); LWRITE(Hs); SWRITE("\n");

    Word Gs1 = NIL; // the list of all derivatives constructed.

    // derivatives are stored is Gs
    Word g_len = np * 2;
    Word g_count = 0;
    Word** Gs = (Word**) malloc(g_len * sizeof(Word*));

    Word* Ps[np]; // array of pointers to elements of Fs, first element is equal to h1
    Word Degrees[np]; // degrees of input polynomials, gives the max index
    Word Ms[np]; // how many steps until we update v
    Word V_index[np];
    Word Chase_index[np]; // chase indices TODO delete

    Word init_v = i0 + 1; // first variable to differentiate on this round
    Word p_index = 0; // initial index, variable to be differentiated in h1 to obtain s1

    // initialise Ps, Degrees, Ms and variable indices
    while (p_index < np) {
        // for each input polynomial P_i
        Word* F1s = Fs[p_index];
        Word D = REDI(SECOND(F1s[0]), i0);

        // initialise "chaser" lists, allows quickly finding function h1 such that s1 = partial_i h1
        Ps[p_index] = F1s;

        // Degree of initial polynomial, gives max index
        Degrees[p_index] = D;

        // initial differentiation variable
        V_index[p_index] = init_v;
        Chase_index[p_index] = 0;

        // M: maximum order of derivative at variable 1
        Ms[p_index] = FIRST(D);

        // increment
        ++p_index;
    }

    // take derivatives
    Word n_finished = 0; // how many are finished
    Word count = 0; // number of derivatives computed so far, can be viewed as a scalar representation of the index

    while (n_finished < np) {
        // cycle back to the beginnig of polynomial list when we reached the end
        if (p_index == np) {
            p_index = 0;
            ++count;
            printf("increment count %d. cycle polynomial list\n", count);

            n_finished = 0; // start counting finished from the start
        }

        // differentiation variable
        Word v = V_index[p_index];
        printf("v = %d\n", v);
        Word ch_index = Chase_index[p_index];

        // no more derivatives to take -- skip.
        if (v > r) {
            printf("skipping %d\n", p_index);
            ++n_finished;
            ++p_index;

            continue;
        }

        // polynomial h_1 and its degree
        Word D = Degrees[p_index];
        Word P = FIRST(*Ps[p_index]);

        // these values are used if P == 0
        Word Q = 0;
        Word Qdeg = D;

        if (P != 0) { // otherwise we have to take derivatives
            // s1 = partial h_1 / x_i
            // TODO partial differential operator with h1,...,hk
            // throw out constant polynomials
            Q = IPDER(r, P, v);
            Qdeg = DEG(r,Q);
            if (LSUM(Qdeg) == 0) { // constant -- cheaper than ipconst
                Q = 0;
            }

            // TODO debugging
            IWRITE(count); SWRITE(": compute derivative, variable "); IWRITE(v); SWRITE(" polynomial "); IWRITE(p_index);
            SWRITE(", degree "); LWRITE(Degrees[p_index]); SWRITE("\n  ");
            LWRITE(INDEX(ch_index, D)); SWRITE(" ");
            P == 0 ? SWRITE("0") : LWRITE(P); SWRITE("\n  ");
            LWRITE(INDEX(count, D)); SWRITE(" "); Q == 0 ? SWRITE("0") : LWRITE(Q); SWRITE("\n\n");

            // TODO recursion
            if (Q != 0) {
                printf("recurse.\n");
                // v = i1, P = h1
                Word Gs2 = STRAT(g_count, r, Gs, COMP(v, Is), COMP(P, Hs));
                Gs1 = CONC(Gs1, Gs2);

                // add derivative Q to Gs
                if (g_count >= g_len) { // enlarge list
                    printf("realloc() %d\n", g_len);
                    g_len += np;
                    Gs = (Word**) realloc(Gs, g_len * sizeof(Word*));
                }

                printf("add derivative to Gs %d\n", g_count);
                Gs[g_count] = Allocate(Q, Qdeg);
                ++g_count;
            }
        }

        // append new derivative Q
        Fs[p_index][count] = LIST2(Q, Qdeg);

        // update variable v and chaser list
        if (count < Ms[p_index]) { // don't need to rollover
            // same variable, higher order -- get next element in "chase" list
            Ps[p_index] = Ps[p_index] + 1;
            Chase_index[p_index] = Chase_index[p_index] + 1;
        } else { // need to rollover
            V_index[p_index] = v + 1; // next variable
            printf("increment differentiation variable %d\n", v);
            Ps[p_index] = Fs[p_index]; // back to beginning of chase list
            Ms[p_index] = Ms[p_index] * LELTI(Degrees[p_index], v - i0); // how many derivatives until we need to roll over again
            Chase_index[p_index] = 0;
        }

        // next polynomial
        ++p_index;
    }
    printf("finished loop\n");

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

    // initial i0 = FIRST(I1) = 0. h0 = FIRST(Hs) = 0
    STRAT(k, r, Fs, LIST1(0), LIST1(0));

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

