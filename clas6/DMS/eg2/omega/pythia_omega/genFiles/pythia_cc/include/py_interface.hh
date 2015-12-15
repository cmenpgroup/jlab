#ifndef py_interface_hh 
#define py_interface_hh 1

#include "cpp_headers.hh"

struct pyjets
{
  int n;
  int npad;
  int k[5][4000];
  double p[5][4000];
  double v[5][4000];    
};


struct pydat1
{
  int    mstu[200];
  double paru[200];
  int    mstj[200];
  double parj[200];
};


struct pydat2
{
  int    kchg[4][500];
  double pmas[4][500];
  double parf[2000];
  double vckm[4][4];
};


struct pydat3
{
  int    mdcy[3][500];
  int    mdme[2][8000];
  double brat[8000];
  int    kfdp[5][8000];
};


struct pydat4
{
  string chaf[2][500];
};


struct pydatr
{
  int    mrpy[6];
  double rrpy[100];
};


struct pysubs
{
  int    msel;
  int    mselpd;
  int    msub[500];
  int    kfin[81][2];
  double ckin[200];
};


struct pypars
{
  int    mstp[200];
  double parp[200];
  int    msti[200];
  double pari[200];
};


struct pyint1
{
  int    mint[400];
  double vint[400];
};


struct pyint2
{
  int    iset[500];
  int    kfpr[2][500];
  double coef[20][500];
  int    icol[2][4][40];
};


struct pyint3
{
  double xsfx[81][2];
  int    isig[3][1000];
  double sigh[1000];
};


struct pyint4
{
  int    mwid[500];
  double wids[5][500];
};


struct pyint5
{
  int   ngenpd;
  int   ngen[3][501];
  float xsec[3][501];
};


struct pyint6
{
  string proc[501];
};


struct pyint7
{
  double sigt[7][7][6];
};


struct pyint8
{
  double xpvmd[13];
  double xpanl[13];
  double xpanh[13];
  double xpbeh[13];
  double xpdir[13];
};


struct pyint9
{
  double vxpvmd[13];
  double vxpanl[13];
  double vxpanh[13];
  double vxpdgm[13];
};


extern "C" pyjets pyjets_;
extern "C" pydat1 pydat1_;
extern "C" pydat2 pydat2_;
extern "C" pydat3 pydat3_;
extern "C" pydat4 pydat4_;
extern "C" pydatr pydatr_;
extern "C" pysubs pysubs_;
extern "C" pypars pypars_;
extern "C" pyint1 pyint1_;
extern "C" pyint2 pyint2_;
extern "C" pyint3 pyint3_;
extern "C" pyint4 pyint4_;
extern "C" pyint5 pyint5_;
extern "C" pyint6 pyint6_;
extern "C" pyint7 pyint7_;
extern "C" pyint8 pyint8_;
extern "C" pyint9 pyint9_;


// pythia functions 
extern "C" double  pygive_(const char*, int);
//extern "C" double  pyname_(int&, const char*, int);
extern "C" void    pylist_(int&);
extern "C" void    pystat_(int&);
//extern "C" double  pymass_(int&);
//extern "C" void    rndmq_(int&, int&, int&, const char*, int);
extern "C" void    pyinit_(const char*, const char*, const char*, double&, int, int, int);
extern "C" void    pyevnt_();
//extern "C" void    pyrget_(int&, int&);
//extern "C" void    pyrset_(int&, int&);
//extern "C" void    pyrobo_(int&, int&, double&, double&, double&, double&, double&);

#endif
