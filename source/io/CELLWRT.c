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
       SWRITE("Signs of Projection Factors\n");
       S = LELTI(c,SIGNPF);
       for (i = 1; i <= k; i++)
         {
         S1 = LELTI(S,k-i+1);
         SWRITE("Level "); GWRITE(i); SWRITE("  : ");
         if (S1 == 0)
           SWRITE("Not determined");
         else
           SIGNLWR(S1);
         SWRITE("\n");
         }

Step3: /* Finish. */
       SWRITE("\n----------------------------------------------------\n");
       goto Return;

Return: /* Prepare for return. */
       return;
}

