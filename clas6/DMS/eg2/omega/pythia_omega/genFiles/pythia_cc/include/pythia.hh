#ifndef pythia_hh
#define pythia_hh 1

#include "cpp_headers.hh"

class Pythia
{
 public:
  Pythia();
  ~Pythia();

  void init();
  void event(int);
};

#endif
