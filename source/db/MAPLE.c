#include "MAPLE.h"
#include <iostream>
#include <string>
using namespace std;

MapleServer::MapleServer(string Base)
{
    string call = Base + "/maple";

    // Fork
    childpid = fork();

    // check fork was successful
    if (childpid == -1) { perror("Failed to fork."); exit(1); }

    // initialise only once.
    if (childpid != 0) return;

    pin.setStdinToPipe();
    pout.setStdoutToPipe();

    execlp(
        call.c_str(),
        call.c_str(),
        "-q",
        NULL
    );

    perror("MapleServer(): Maple initialisation failed.");
    pout.closeOut();
    exit(0);
}

MapleServer::~MapleServer()
{
    pin.out() << endl << ":quit" << endl << flush;
}

string RatToString(Word n)
{
    // a hack to allow us to use SACLIB's writing functions to create strings
    string out;
    {
        ostringstream sout;
        PushOutputContext(sout);
        RNWRITE(n);
        PopOutputContext();
        out = sout.str();
    }

    return out;
}

// testing method to add two rational numbers.
Word MapleServer::test(Word a, Word b)
{
    // write the expression into maple
    std::cout << RatToString(a) << " + " << RatToString(b) << ";";
    pin.out() << RatToString(a) << " + " << RatToString(b) << ";" << endl;
    std::cout << "successful write to maple.";

    // read maple output
    string out;
    pout.in() >> out;
    std::cout << "successful read from maple.";
    std::cout << out;
    return 0;
}


