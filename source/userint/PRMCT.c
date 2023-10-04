/*======================================================================
                      PRMCT()

Process MCT command.
======================================================================*/
#include "qepcad.h"

void QepcadCls::PRMCT()
{
    Word C;
    /* hide C; */

Step1: /* Get the argument. */
    C = CREADB();

Step2: /* Check for the validity. */
    if (C != 'y' && C != 'n') {
        SWRITE("Error PRMCT 'y' or 'n' was expected.\n");
        DIELOC();
        goto Return;
    }

Step3: /* Set. */
    PCMCT = C;
    PCFULL = C;
    goto Return;

Return: /* Prepare for return. */
    return;
}
