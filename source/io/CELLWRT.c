/*======================================================================
                      CELLWRT(c)

Cell Write (Tarski).

Write out the cell C as a tarski formula
======================================================================*/
#include "qepcad.h"

void QepcadCls::CELLWRT(Word c)
{
       Word S,S1,k,t,i,D,M;
       /* hide t; */

Step1: /* Heading. */
       k = LELTI(c,LEVEL);
       SWRITE("---------- Information about the cell "); LWRITE(LELTI(c,INDX));
       SWRITE(" ----------\n\n");

Step2: /* Signs of Projection Factors. */
       S = LELTI(c,SIGNPF);
       CELLIPLLDWR(GVVL, GVPF, S, k); SWRITE("\n"); goto Return;

Step3: /* Finish. */
       SWRITE("\n----------------------------------------------------\n");
       goto Return;

Return: /* Prepare for return. */
       return;
}

