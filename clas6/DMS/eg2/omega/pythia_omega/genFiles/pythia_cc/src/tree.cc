#include "cpp_headers.hh"
#include "tracks.hh"
#include "root.hh"

void initTree()
{
  tr->Branch("Ntracks",      &trk.Ntracks,      "Ntracks/I");

  tr->Branch("Eb",           &trk.Eb,           "Eb/D");
  tr->Branch("tarA",         &trk.tarA,         "tarA/I");
  tr->Branch("tarZ",         &trk.tarZ,         "tarZ/I");
  tr->Branch("process",      &trk.process,      "process/I");
  tr->Branch("nu",           &trk.nu,           "nu/D");

  tr->Branch("type",        trk.type,         "type[Ntracks]/I");
  tr->Branch("parent",      trk.parent,       "parent[Ntracks]/I");
  tr->Branch("px",          trk.px,           "px[Ntracks]/D");
  tr->Branch("py",          trk.py,           "py[Ntracks]/D");
  tr->Branch("pz",          trk.pz,           "pz[Ntracks]/D");
  tr->Branch("p",           trk.p,            "p[Ntracks]/D");
  tr->Branch("E",           trk.E,            "E[Ntracks]/D");
} 


void init_trk()
{
  trk.process  = -10000;
  trk.nu       = -10000;

  for(int ii=0;ii<MAX_T;ii++)
  {
    trk.type[ii]   = -10000;
    trk.parent[ii] = -10000;
    trk.px[ii]     = -10000;
    trk.py[ii]     = -10000;
    trk.pz[ii]     = -10000;
    trk.p[ii]      = -10000;
    trk.E[ii]      = -10000;
  }
}
