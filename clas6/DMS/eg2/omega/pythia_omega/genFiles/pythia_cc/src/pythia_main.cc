#include "cpp_headers.hh"
#include "root.hh"
#include "pythia.hh"
#include "py_interface.hh"
#include "read_optFile.hh"
#include "funcs.hh"
#include "tracks.hh"
#include "filter.hh"

ReadOptFile *rof;
TRACKS trk;
TTree *tr;

int main(int argc, char *argv[])
{
  ios::sync_with_stdio(true);

  if(argc != 2) {help_mess(); return 0;}

  //parameters file
  string OptFile(argv[1]);

  // read user parameters
  rof = new ReadOptFile();
  rof->optfile = OptFile;
  rof->read_optfile();

  TFile *out_file = TFile::Open("pythia.root","RECREATE");
  if(!out_file->IsOpen())
  {
    cout << "Failed to open ntuple file " << endl;
    exit(1);
  }
  out_file->SetFileBytesWritten(0); // reset written bytes counter
  TTree::SetMaxTreeSize(2000000000);

  tr = new TTree("MC","hadr_tree");
  initTree();

  trk.Eb      = rof->fEe;
  trk.tarA    = rof->fAt;
  trk.tarZ    = rof->fZt;

  //Pythia initialization
  Pythia *py = new Pythia();

  // event generation
  cout<<"start event generation."<<endl;

  // filter initialization
  Filter *filt = new Filter();

  for(int ie=1; ie<=rof->fNevts; ie++)
  {
    if((ie%1000) == 0) cout << "generating event : " << ie << "/" << rof->fNevts << endl;

    py->event(ie);

    if(filt->Cut()) tr->Fill();

  }//end of loop over events

  //...Print cross sections.
  int i1=1;
  int i4=4;
  pystat_(i1);
  pystat_(i4);

  cout<<"============================================"<<endl;
  cout<<"Pythia total x sec normalisation: "<<pypars_.pari[0]*1000<<" microbarn" << endl;
  cout<<"Total Number of generated events: "<<pypars_.msti[4]<<endl;
  cout<<"Total Number of trials: "<<pyint5_.ngen[2][0]<<endl;
  cout<<"============================================"<<endl;

  out_file->Write();
  out_file->Close();

  delete out_file;
  delete rof;

  return 0;
}
