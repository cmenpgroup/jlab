#ifndef filter_hh
#define filter_hh 1

#include "cpp_headers.hh"

class Filter
{
 int KScut;

 public:
  Filter();
  ~Filter();

  void init();
  void SetKScut(int ks);
  int GetKScut() {return KScut;};
  bool Cut();
};

#endif
