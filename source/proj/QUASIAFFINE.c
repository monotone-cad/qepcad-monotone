/*======================================================================
  A -< QUASIAFFINE(A, r; P)

Add partial differentials of (factorised) input polynomials to ensure cells are quasi-affine. Projection factors are
sufficient for smoothness and such that each cell will either be increasing, decreasing or constant along each
coordinate axis.

Input
  \parm{A} projection factor set
  \parm{r} number of variables

Output
  \parm{A*}: projection factor set with added derivatives

SideEffect
  \parm{A} is modified.

  ======================================================================*/
#include "qepcad.h"

// convenience function for multi-degree of polynomial, returns (d_r,...,d_1) where d_i is degree of P in variable i
// returns degree + 1
inline Word DEG(Word r, Word P)
{
    Word D = PDEGV(r, P);
    Word d, D1 = NIL;

    while (D != NIL) {
        ADV(D, &d, &D);

        D1 = COMP(d+1, D1);
    }

    return D1;
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

// list product
inline Word LPROD(Word L)
{
    Word m = 1, a = 0;
    while (L != NIL) {
        ADV(L, &a, &L);

        // multiply by zero short-circuits
        if (a == 0) return 0;

        // otherwise perform the multiplication
        m *= a;
    }

    return m;
}


// recover index from count
// TODO this is only for debugging, can be deleted when done.
// count = (c1 + c2 * m1 + ... + cn * m(n-1))
Word IndexHelper(Word m, Word M, Word* count_)
{
    Word count = *count_;
    Word a, M1;
    ADV(M, &a, &M1);

    if (M1 == NIL) {
        *count_ = count % m;

        return LIST1(count / m);
    }

    // skip round
    if (a == 1) {
        return COMP(0, IndexHelper(m, M1, count_));
    }

    Word I1 = IndexHelper(m * a, M1, count_);
    count = *count_;

    Word j = count / m;
    *count_ = count % m;

    return COMP(j, I1);
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
    Word len = LPROD(D) + 2; // includes P itself and a null at the end
    Word* F1s = (Word*) calloc(len, sizeof(Word));
    if (F1s == NULL) {
        FAIL("QUASIAFFINE", "calloc F1s failed.");
    }

    // first element is P1, rest are 0 (by calloc)
    F1s[0] = P1;

    return F1s;
}


// construct a jacobi (k + 1) * (k + 1)-matrix from a minor (k*k-matrix) by appending one column
// (\H_k / \x_j,...,\H_1) and one row (\P / \x_j, \P / \x_ik,..., \P / \x_i1)
Word JacobiFromMinor(Word r, Word P, Word j, Word Hs, Word Is, Word Minor)
{
    // Step 1: add column
    Word Jacobi = NIL; // list of k rows
    Word k = 0; // dimension of Minor

    Word H, Row;
    while (Minor != NIL) {
        ADV(Minor, &Row, &Minor);
        ADV(Hs, &H, &Hs);
        ++k;

        // compute \H / \x_j
        Word H1 = IPDER(r, H, j);

        // append H1 to row
        Row = COMP(H1, Row);
        Jacobi = COMP(Row, Jacobi);
    }

    // Jacobi was constructed backwards
    Jacobi = INV(Jacobi);

    Word i;
    Row = NIL;
    while (k > 0) { // Is is assumed to have k+1 elements
        ADV(Is, &i, &Is);
        --k;

        // compute \P / \x_i and append
        Word P1 = IPDER(r, P, i);
        Row = COMP(P1, Row);
    }

    // row constructed in order (i1,...,ik)
    Row = INV(Row);

    // compute \P / \x_j and append
    Word P1 = IPDER(r, P, j);
    Row = COMP(P1, Row);

    // add row to jacobi matrix and return
    return COMP(Row, Jacobi);
}

// recursive smooth stratification of polynomials Fs
// np: number of polynomials
// r: number of variables
// Fs: list of polynomials and degrees (..., P_i, deg(P_i), ...)
// Is: indices (i1,...,ik)
// Hs: list of polynomials h_i1,...,h_ik
// Minor: (LENGTH(I1) - 1) * (LENGTH(I1) - 1) - matrix |a_ij| = partial h_i / partial x_j for h in Hs and i in Is
// return list of all differentials computed by alrogithm
Word STRAT(Word np, Word r, Word*** Fs_, Word Is, Word Hs, Word Minor)
{
    // end of recursion
    Word i0 = FIRST(Is); // number of variables considered so far
    if (np == 0 || i0 == r) {
        return NIL;
    }

    // TODO debugging
    SWRITE("STRAT: Is = "); LWRITE(Is); SWRITE("\n");
    SWRITE("       Hs = "); LWRITE(Hs); SWRITE("\n");

    // Dereference Fs_
    Word **Fs = *Fs_;

    // set up return value
    Word Gs1 = NIL; // list of all differentials computed in this round, to return

    // set up working array
    Word g_len = np * 2; // length of memory allocated for Gs
    Word g_count = 0; // how many differentials computed so far, index in Gs
    Word** Gs = (Word**) calloc(g_len, sizeof(Word*)); // list of all differentials computed in this round, working set
    if (Gs == NULL) {
        FAIL("QUASIAFFINE", "calloc Gs failed.");
    }

    // metadata
    Word* Ps[np]; // "chase list", first element is h1
    Word Degrees[np]; // degree of Fs[0][0], gives max index
    Word Ms[np]; // number of steps before maximum differentiation variable (i1) should be incremented
    Word Dvs[np]; // list of current differentiation variable (i1)
    Word Chase[np]; // TODO debugging

    Word p_index = 0; // initial polynomial index, ranges over 0 <= p_index < np

    // initialise metadata for each input polynomial
    while (p_index < np) {
        Word* F1s = Fs[p_index];
        Word D = REDI(SECOND(F1s[0]), i0);

        Ps[p_index] = F1s;
        Degrees[p_index] = D;
        Ms[p_index] = COMP(1, LCOPY(D)); // neet do copy as we modify later
        Dvs[p_index] = i0;
        Chase[p_index] = 0;

        // increment index
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
        }

        Word v = Dvs[p_index]; // differentiation variable
        Word m = FIRST(Ms[p_index]);

        // update variable v and chaser list
        if (count >= m && v == r) { // rollover, but we're finished
            ++n_finished; // this polynomial is done.
            ++p_index; // next polynomial

            continue; // skip it
        } else if (count >= m) { // rollover - increment differentiation variable
            ++v; // next variable ...
            Dvs[p_index] = v; // ... and store
            Ps[p_index] = Fs[p_index]; // back to beginning of chase list
            Chase[p_index] = 0;

            // calculate next m
            Word M1 = RED(Ms[p_index]);
            Word d = FIRST(M1);
            SFIRST(M1, d * m);
            Ms[p_index] = M1;

            // degree zero - no derivatives taken for this variable. next iteration will increment the variable.
            if (d == 1) continue;
        }
        printf("p_index %d, variable %d, count = %d, chase_index = %d\n", p_index, v, count, Chase[p_index]);

        // compute s_k = partial_{(h_1,...,h_{k-1}),(i_1,...,i_{k-1}),v} h_k
        // get h_k and its degree
        Word D = Degrees[p_index];
        Word P = FIRST(*Ps[p_index]);

        // construct jacobi matrix using h1 = P and i1 = v
        Word Jacobi = JacobiFromMinor(r, P, v, Hs, Is, Minor);

        // compute partial differential, determinant of jacobi matrix
        Word Q = MAIPDE(r, Jacobi); // next derivative is the jacobi determinant
        Word Qdeg = DEG(r,Q);

        if (LSUM(Qdeg) == r) { // s_k is constant, may as well be zero
            Q = 0;
        }

        // TODO debugging
        IWRITE(count); SWRITE(", variable: "); IWRITE(v); SWRITE(" polynomial: "); IWRITE(p_index);
        SWRITE(", degree: "); LWRITE(Degrees[p_index]); SWRITE("\n  ");
        LWRITE(INDEX(count, D)); SWRITE(" ");
        P == 0 ? SWRITE("0") : LWRITE(P); SWRITE("\n  ");
        LWRITE(INDEX(Chase[p_index], D)); SWRITE(" ");
        Q == 0 ? SWRITE("0") : LWRITE(Q); SWRITE("\n\n");

        if (Q != 0) {
            // Gs2 contains derivatives computed during recursion
            Word Gs2 = STRAT(g_count, r, &Gs, COMP(v, Is), COMP(P, Hs), Jacobi);
            Gs1 = CONC(Gs1, Gs2);

            // append new derivative to working list Gs and return list
            if (g_count >= g_len) { // enlarge list
                g_len += np;
                Gs = (Word**) reallocarray(Gs, g_len, sizeof(Word*));
            }

            Gs[g_count] = Allocate(Q, Qdeg);
            ++g_count;
        }

        // append derivative to Fs, preserving zeroes
        Fs[p_index][count] = LIST2(Q, Qdeg);

        // next polynomial please.
        Ps[p_index] = Ps[p_index] + 1;
        Chase[p_index] = Chase[p_index] + 1;
        ++p_index;
    }

    // construct list Gs1 of all functions in Gs
    Word i = g_count;
    while (i > 0) { // comp backwards
        --i;
        // future functions expect QEPCAD polynomials
        Word P1 = MPOLY(FIRST(Gs[i][0]), NIL, NIL, PO_POLY, PO_KEEP);
        Gs1 = COMP(P1, Gs1);
    }

    // free Gs
    for (int i = 0; i < g_count; ++i) {
        free(Gs[i]);
    }

    free(Gs);

    // returns list of derivatives
    *Fs_ = Fs;
    return Gs1;
}

// list of all partials, sufficient for smooth stratification
Word PARTIALS(Word r, Word L)
{
    Word D, P, P1, k, j;

    k = LENGTH(L);
    Word** Fs = (Word**) calloc(k, sizeof(Word*));
    if (Fs == NULL) {
        FAIL("QUASIAFFINE", "calloc Fs failed.");
    }

    j = 0;
    while (L != NIL) {
        ADV(L, &P, &L);

        D = DEG(r, P);
        Fs[j] = Allocate(P, D);

        ++j;
    }

    // initial i0 = FIRST(I1) = 0. h0 = FIRST(Hs) = 0, Minor is the empty matrix
    Word Gs = STRAT(k, r, &Fs, LIST1(0), LIST1(0), NIL);

    // free 2d Fs
    for (j = 0; j < k; j++) free(Fs[j]);
    free(Fs);

    return Gs;
}

void QepcadCls::QUASIAFFINE(Word A, Word r, Word *A_)
{
    Word AA, A1, A11, P, k, r1, L;

Step1: /* decide based on dimension */
    if (r <= 1) {
        // nothing to do
        *A_ = A;
        return;
    }

    // collect every polynomial in the input set, every polynomial will be in r variables regardless of lever
    AA = A, L = NIL, k = 0, r1 = r;
    while (AA != NIL) {
        ADV(AA, &A1, &AA);
        ++k; --r1;

        while (A1 != NIL) {
            ADV(A1, &A11, &A1);
            P = LELTI(A11, PO_POLY);

            L = COMP(PADDVS(LELTI(A11, PO_POLY), r1), L);
        }
    }

    Word LL = PARTIALS(r, L);
    // factorise and append -- same function as used on input formula
    ADDPOLS(IPLFAC(k, LL),k,LFS("D"), &A);

    // put proj fac in correct order
    *A_ = A;
}

