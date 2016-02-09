#include "pythia.hh"
#include "py_interface.hh"
#include "funcs.hh"
#include "tracks.hh"
#include "read_optFile.hh"
#include "root.hh"

Pythia::Pythia()
{
  init();
}

Pythia::~Pythia()
{;}

void Pythia::init()
{
  double fEp = m_prot();
  double S=sqrt(4.*fEp*rof->fEe);
  cout<<"*********************************************"<<endl;
  cout<<"proton beam energy: "<< fEp <<" GeV"<<endl;
  cout<<"lepton beam energy: "<< rof->fEe <<" GeV"<<endl;
  cout<<"resulting sqrt(s):  "<<S<<" GeV"<<endl;
  cout<<"*********************************************"<<endl;

  //lepton beam
  pyjets_.p[0][0] =  0.0;
  pyjets_.p[1][0] =  0.0;
  pyjets_.p[2][0] = rof->fEe;

  //proton target
  pyjets_.p[0][1] = 0.0;
  pyjets_.p[1][1] = 0.0;
  pyjets_.p[2][1] = 0.0;

  //fix target experiment
  string fixt = "FIXT";

  //electron beam
  string ge  = "gamma/e-";

  string p       = "p+";

  pyinit_(fixt.c_str(), ge.c_str(), p.c_str(), rof->fEe, fixt.size(), ge.size(), p.size());
}


void Pythia::event(int ie)
{
  init_trk();

  //generate event
  pyevnt_();

  while(pypars_.msti[60] == 1)
  {
    cout<<"call pyevt again"<<endl;
    pyevnt_();  
    break;    
  }

  TLorentzVector Beam;
  TLorentzVector Electron;
  TLorentzVector BeamMinusElectron;

  double elecM = m_elec();
  Beam.SetXYZM(0.0,0.0,rof->fEe,elecM);
    
  trk.Ntracks = pyjets_.n;
  trk.process = pypars_.msti[0];
  trk.nu      = rof->fEe - pyjets_.p[3][2];

  //loop over tracks
  for(int i=0; i<trk.Ntracks; i++)
  {
    trk.ks[i]     = pyjets_.k[0][i];
    trk.type[i]   = pyjets_.k[1][i];
    trk.parent[i] = pyjets_.k[2][i];
    trk.px[i]     = pyjets_.p[0][i];
    trk.py[i]     = pyjets_.p[1][i];
    trk.pz[i]     = pyjets_.p[2][i];

    trk.p[i]      = sqrt(sqr(trk.px[i])+sqr(trk.py[i])+sqr(trk.pz[i]));
    trk.E[i]      = pyjets_.p[3][i];
      
    if(trk.type[i]==pdg_elec()){
        Electron.SetPxPyPzE(trk.px[i],trk.py[i],trk.pz[i],trk.E[i]);
    }
  }
  BeamMinusElectron = Beam - Electron;
  trk.Q2 = -BeamMinusElectron.Mag2();
}
