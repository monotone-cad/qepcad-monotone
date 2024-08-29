#ifndef _SINGSAC_POLICY_
#define _SINGSAC_POLICY_
#include "db/OriginalPolicy.h"
#include "db/SingularPolicy.h"

class SingSacPolicy : public CAPolicy
{
  CAPolicy *sing, *sac;
public:
  SingSacPolicy()
  {
    sac = new OriginalPolicy;
    sing = new SingularPolicy;
  }
  ~SingSacPolicy()
  {
    delete sing;
    delete sac;
  }
  void IPFAC(Word r, Word P, Word *s_, Word *c_, Word *L_)
  {
    if (r > 1)
      sing->IPFAC(r,P,s_,c_,L_);
    else
      sac ->IPFAC(r,P,s_,c_,L_);
  }

  Word IPRES(Word r, Word A, Word B)
  {
    if (r > 2)
      return sing->IPRES(r,A,B);
    else
      return sac->IPRES(r,A,B);
  }

  Word IPDSCR(Word r, Word A)
  {
    if (r > 2)
      return sing->IPDSCR(r,A);
    else
      return sac ->IPDSCR(r,A);
  }

  Word IPFACTGB(Word r, Word I, Word N)
  {
    return sing->IPFACTGB(r,I,N);
  }

  Word CONSTORDTEST(Word r, Word A, Word L)
  {
    return sing->CONSTORDTEST(r,A,L);
  }
  Word GROEBNER(Word Ps, Word Rs, Word r) {
    return sing->GROEBNER(Ps, Rs, r);
  }


  const string name() { return "SingSacPolicy"; }
  bool supports(const string &s)
  {
      return s == "IPFAC"
          || s == "IPRES"
          || s == "IPDSCR"
          || s == "IPFACTGB"
          || s == "CONSTORDTEST"
          || s == "GROEBNER";
  }

};

#endif
