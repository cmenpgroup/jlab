#ifndef filter_hh
#define filter_hh 1

#include "cpp_headers.hh"

class Filter
{
 public:
  Filter();
  ~Filter();

  void init();
  bool Cut();
};

#endif
