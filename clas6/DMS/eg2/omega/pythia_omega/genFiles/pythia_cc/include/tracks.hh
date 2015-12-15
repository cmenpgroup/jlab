#ifndef tracks_hh
#define tracks_hh 1

#include "cpp_headers.hh"
#include "root.hh"

#define MAX_T 200

class TRACKS
{
 public:
  int    Ntracks;        //
  double Eb;             //
  int    tarA;           //
  int    tarZ;           //
  int    process;        //
  double nu;
  double y;

  int    type[MAX_T];    //
  int    parent[MAX_T];  //
  double px[MAX_T];      //
  double py[MAX_T];      //
  double pz[MAX_T];      //
  double p[MAX_T];       //
  double E[MAX_T];       //
};

extern TRACKS trk;
extern TTree *tr;

#endif
  