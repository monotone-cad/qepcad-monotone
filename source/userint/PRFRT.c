/*======================================================================
                      PRFRT()

Process FRT command.
======================================================================*/
#include "qepcad.h"

void QepcadCls::PRFRT()
{
    Word C;
    /* hide C; */

Step1: /* Get the argument. */
    C = CREADB();

Step2: /* Check for the validity. */
    if (C != 'y' && C != 'n') {
        SWRITE("Error PRFRT 'y' or 'n' was expected.\n");
        DIELOC();
        goto Return;
    } else if (LENGTH(GVVL) > 3) { // only implemented in up to R^3
        SWRITE("Error PRFRT only implemented for up to R^3\n");
        DIELOC();
        goto Return;
    }

Step3: /* Set. */
    PCFRT = C;
    PCFULL = C;
    goto Return;

Return: /* Prepare for return. */
    return;
}
