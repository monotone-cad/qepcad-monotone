#include "qepcad.h"

// helper function: convert sample to primitive, no timekeeping
void ConvertToPrimitive(Word SQ, Word SJ, Word SM, Word SI, Word Sb, Word* M_, Word* I_, Word* b_)
{
    Word G, K, t, u, a, b, junk;
    SIMPLEQE(SM,SI,SQ,SJ,&G,&t,&u,&K,&a,&b,&junk,&junk);

    if (u != 0) {
        MODCRDB(Sb,G,a,b,&Sb);
    } else {
        Sb = CCONC(Sb,LIST1(b));
    }

    *M_ = G;
    *I_ = K;
    *b_ = Sb;
}


