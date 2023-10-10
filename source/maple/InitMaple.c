/*======================================================================
 * kern = InitMaple()
 *
 * initialise, returning a maple kernel
 *====================================================================*/
#include "qepcad.h"

static void M_DECL textCallBack( void *data, int tag, const char *output )
{
    printf("%s\n",output);
}

MKernelVector InitMaple()
{
    // init code from "./samples/OpenMaple/simple/simple.c"
    char err[2048];  /* command input and error string buffers */
    MKernelVector kv;  /* Maple kernel handle */
    MCallBackVectorDesc cb = {
        textCallBack,
        0,   /* errorCallBack not used */
        0,   /* statusCallBack not used */
        0,   /* readLineCallBack not used */
        0,   /* redirectCallBack not used */
        0,   /* streamCallBack not used */
        0,   /* queryInterrupt not used */
        0    /* callBackCallBack not used */
    };

    /* initialise Maple */
    if((kv=StartMaple(0,NULL,&cb,NULL,NULL,err)) == NULL) {
        FAIL("InitMaple", "Maple initialisation failed.");
    }

    // load packages
    LoadMaplePackage(kv, "RootFinding");
    LoadMaplePackage(kv, "Groebner");

    return kv;
}

