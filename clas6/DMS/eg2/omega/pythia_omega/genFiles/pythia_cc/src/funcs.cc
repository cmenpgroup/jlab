#include "cpp_headers.hh"
#include "py_interface.hh"

void help_mess()
{
  cout<<"Usage:"<<endl;
  cout<<"./pythia.exe <Option_file>"<<endl;
}

double m_prot() {return 0.93827;}
double m_elec() {return 0.00051;}
double pi() {return 3.14159265359;}
double sqr(double x) {return x*x;}
double todeg(double a) {return a*180./pi();}
double torad(double a) {return a*pi()/180.;}



