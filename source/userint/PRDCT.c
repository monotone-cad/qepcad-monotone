/*======================================================================
                      PRDCT()

Process "display cell as tarski formula" command.
======================================================================*/
#include "qepcad.h"

void QepcadCls::PRDCT()
{
       Word c,t;
       /* hide t; */

Step1: /* Read in an argument. */
       CELLRDR(GVPC,&c,&t); if (t == 0) goto Return;

Step2: /* Display the cell info. */
       printf("WIP hello!");

Return: /* Prepare for return. */
       return;
}

