/***************************************************************
 * Define a CA system based on Maple
 ***************************************************************/
#ifndef _MAPLE_
#define _MAPLE_

#include "qepcad.h"
#include "db/unnamedpipe.h"
#include <sstream>
#include <signal.h>

class MapleServer
{
public:
  UnnamedPipe pin, pout;
  pid_t childpid;

  MapleServer(string Base);
  ~MapleServer();
  void kill() { ::kill(childpid,SIGKILL); }

  Word test(Word a, Word b);

  const string name() { return "Maple"; }
};


#endif
