#include "cpp_headers.hh"
#include "options.hh"
#include "read_optFile.hh"
#include "py_interface.hh"
#include "funcs.hh"

ReadOptFile::ReadOptFile()
{;}


ReadOptFile::~ReadOptFile()
{;}


int ReadOptFile::read_optfile()
{
  int ir = 0;
  vector<int>    itmp1;
  vector<int>    itmp2;
  vector<double> dtmp1;
  vector<double> dtmp2;
  vector<string> stmp1;

  parse_file(optfile.c_str());

  //Lepton Beam Type
  get_opt("SET", "BTYPE", itmp1);
  if(itmp1.size() != 1)
  {
    cerr<<"BTYPE: missing data "<<itmp1.size() <<endl; ir = -1; exit(0);
  }
  fLtype = itmp1[0];
  itmp1.clear();

  //Lepton Beam Energy
  get_opt("SET", "BENERGY", dtmp1);
  if(dtmp1.size() != 1)
  {
    cerr<<"BENERGY: missing data "<<dtmp1.size() <<endl; ir = -1; exit(0);
  }
  fEe = dtmp1[0];
  dtmp1.clear();

  //TarA, TarZ
  get_opt("SET", "TARGET", dtmp2);
  if(dtmp2.size() != 2)
  {
    cerr<<"TARGET: missing data "<< dtmp2.size() << endl; ir = -1; exit(0);
  }
  fAt = dtmp2[0];
  fZt = dtmp2[1];
  dtmp2.clear();
  if(!(fAt == 1 && (fZt == 0 || fZt == 1)))
  {
    cerr<<"Target not supported"<<endl;ir = -1; exit(0);
  }

  //Number of events to generate
  get_opt("SET", "EVENTS", itmp1);
  if(itmp1.size() != 1)
  {
    cerr<<"EVENTS: missing data "<<itmp1.size() <<endl; ir = -1; exit(0);
  }
  fNevts = itmp1[0];
  itmp1.clear();


  //Number of events per ntuple
  get_opt("SET",   "EVTSPERNTP",  itmp1);
  if(itmp1.size() != 1)
  {
    cerr<<"EVTSPERNTP: missing data "<<itmp1.size() <<endl; ir = -1; exit(0);
  }
  fNevtsPerNtp = itmp1[0];
  itmp1.clear();

  dump_optfile();

  // passing all pythia parameters
  string setpar;
  while(get_all_opt("PYTHIA", "SET", setpar))
  {
    pygive_(setpar.c_str(),setpar.size());
  }

  return ir;
}


void ReadOptFile::dump_optfile()
{
  cout << "Reading option file.............................  " << optfile << endl;

  cout << "Beam Type set to................................  " << fLtype  << endl;
  cout << "Beam Energy set to..............................  " << fEe     << " GeV" << endl;
  cout << "Events to generate set to.......................  " << fNevts << endl;
  cout << "Target A and Z set to...........................  " << fAt <<", "<< fZt << endl;
  cout<<"Number of events per root ntuple set to...........  " << fNevtsPerNtp << endl;

  cout << "Reading input file  DONE WITH SUCCESS " << endl;
}


int ReadOptFile::getStatus()
{
  return status;
}