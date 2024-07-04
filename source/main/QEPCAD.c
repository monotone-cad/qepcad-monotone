/*======================================================================
  QEPCAD(Fs; t,F_e,F_n,F_s)

  Quantifier Elimination by Partial Cylindrical Algebraic Decomposition.

  \Input
  \parm{F*} $=(Q_{f+1} x_{f+1})\cdots(Q_r x_r)\hat{F}(x_1,\ldots,x_r)$,
  $0\leq f < r$, is a quantified formula.

  \Output
  \parm{t}  is either \c{EQU} or \c{INEQU}.
  \parm{Fe} is a quantifier-free formula equivalent to~\v{F*}
  if \v{t} is \c{EQU}, otherwise \v{Fe} is undefined.
  \parm{Fn} is a quantifier-free formula necssary for~\v{F*}
  if \v{t} is \c{INEQU}, otherwise \v{Fn} is undefined.
  \parm{Fs} is a quantifier-free formula sufficient for~\v{F*}
  if \v{t} is \c{INEQU}, otherwise \v{Fs} is undefined.
  ======================================================================*/
#include "qepcad.h"

void QepcadCls::QEPCAD(Word Fs, Word *t_, Word *F_e_, Word *F_n_, Word *F_s_)
{
    Word A,D,F,F_e,F_n,F_s,Fh,J,P,Q,Ths,f,i,r,t, T;
    /* hide Ths,i,t; */
    Word cL,**cC,cr,ce,ci,*cT,cj,cs,cl,ct; /* Chris variables. */
    Word Cs,Ps,Qs,Pps,Cps,Qps,SF; /* Chris variables. */
    char c1,c2; /* Chris variables. */

Step1: /* Normalize. */
    FIRST4(Fs,&r,&f,&Q,&Fh);
    /*Int*/ PCNSTEP = 1;
    /*Int*/ if (INTERACT()) USERINT(LFS("Before Normalization"),'a');
    /*Int*/ if (PCCONTINUE == TRUE) { goto Return; }
    /*Int*/ Ths = ACLOCK();
    F = NORMQFF(Fh);
    if (GVUA != NIL) GVNA = NORMQFF(GVUA);
    /*Int*/ Ths = ACLOCK() - Ths;
    /*Int*/ TMNORMQFF = Ths;
    /*Int*/ GVNQFF = F;
    //       if (TYPEQFF(F) != UNDET) { t = EQU; F_e = F; goto Return; }
    /*Int*/ GVREFL = NIL;
    /*Int*/ GVTD = NIL;

Step2: /* Projection. */
    if (GVUA != NIL) F = LIST3(ANDOP,GVNA,F);
    A = EXTRACT(r,F);
    if (GVUA != NIL) {
        GVNA = SECOND(F);
        F = THIRD(F);
    }

    /*Int*/ for (i = 1; i <= r; i++) NMNIP[i] = LENGTH(LELTI(A,i));
    /*Int*/ GVPF = LLCOPY(A);
    /*Int*/ GVNIP = A;
    /*Int*/ GVLV = r;
    /*Int*/ PCNSTEP = 1;

    /*Int if (INTERACT()) USERINT(LFS("After Normalization"),'A'); */
    /*Int PCNSTEP = 1; */

    // project andadd jacobi determinants for quasi-affine cells
    if (PCMCT == 'y') {
        // note that quasi-affine does projeection, too.
        /*Int*/     USERINT(LFS("Before Projection (quasi-affine)"),'b');
        QUASIAFFINE(r, GVVL, F, &A, &P, &J);
        GVNIP = A;
        GVPF = P;
        GVPJ = J;
    } else { // standard case
        PROJECT(r,A,&P,&J);
    }

    /*Int*/ if (PCCONTINUE == TRUE) { goto Return; }

Step3: /* Truth-invariant CAD. */
    /*Int*/ NMFPF = 0;
    /*Int*/ for (i=1; i<=f; i++) NMFPF=NMFPF+LENGTH(LELTI(P,i));
    /*Int*/ PCNSTEP = 1;
    D = TICAD(Q,F,f,P,A);

Step5: /* Monotone cells, if needed */
    if (PCMCT == 'y') {
        /*Int*/ GVPC = D;
        /*Int*/ PCNSTEP = 1;
        /*Int*/ if (INTERACT()) USERINT(LFS("Before Refinement For Monotone Cells"),'m');
        /*Int*/ if (PCCONTINUE == TRUE) { goto Return; }
        /*Int*/ Ths = ACLOCK();

        // add extra polynomials for [semi]-monotone cells and recompute the cad if needed
        Word RPs = MONOTONE(&P, &J, D, r);
        GVREFL = RPs;

        // refine the CAD
        D = REFINE(1, D, GVREFL, P);

        // if we cached an ESPCAD, it will no longer be valid. If CAD was projection definable before, it will be still
        // be projection definable now.
        if (GVTD != NIL && FIRST(GVTD) == 0) {
            SWRITE("*** Invalidating cached ESPCAD. ***\n");
        }

        D = FRONTIER(f, 1, D, P, &P, &J, &RPs);
        GVREFL = RPs;
        // TODO refine, perhaps make a new refinement list?
        D = REFINE(1, D, GVREFL, P);

        /*Int*/ if (PCCONTINUE == TRUE) { goto Return; }
    }

Step6: /* Solution. */
    /*Int*/ GVPC = D;
    /*Int*/ PCNSTEP = 1;
    /*Int*/ if (INTERACT()) USERINT(LFS("Before Solution"),'e');
    /*Int*/ if (PCCONTINUE == TRUE) { goto Return; }
    T = ACLOCK();
    if (!PCMZERROR)
        SFC3(GVPC,GVPF,GVPJ,GVNFV,CCONC(LIST10(0,0,0,1,0,3,2,4,1,5),LIST1(-1)));
    else
        SFCFULLD(GVPC,GVPF,GVPJ,GVNFV);
    T = ACLOCK() - T;
    TMSFCONST = COMP(T,TMSFCONST);

Return: /* Prepare for return. */
    *t_ = t;
    *F_e_ = F_e;
    *F_n_ = F_n;
    *F_s_ = F_s;
    return;
}
