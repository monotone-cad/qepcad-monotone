/*======================================================================
 * LoadMaplePackage(packageName)
 *
 * loads the package called <packageName> into maple.
 *====================================================================*/
#include "qepcad.h"

void LoadMaplePackage(MKernelVector kv, string name)
{
    EvalMapleStatement(kv, ("with(" + name + ");").c_str());
}

