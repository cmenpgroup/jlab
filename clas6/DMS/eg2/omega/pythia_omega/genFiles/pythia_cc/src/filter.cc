#include "filter.hh"
#include "funcs.hh"
#include "tracks.hh"

Filter::Filter()
{
  init();
}

Filter::~Filter()
{;}

void Filter::init()
{

    cout<<"*********************************************"<<endl;
    cout<<"Initialize the event filter"<<endl;
    cout<<"Particle list: 11, 211, -211, 22"<<endl;
    cout<<"*********************************************"<<endl;

}


bool Filter::Cut()
{
    bool ret = false;
    
    int nElectron = 0;
    int nPip = 0;
    int nPim = 0;
    int nGamma = 0;

    //loop over tracks
    for(int i=0; i<trk.Ntracks; i++){
        if(trk.type[i]==1){
            switch(trk.parent[i]){
                case 11: nElectron++; break;
                case 211: nPip++; break;
                case -211: nPim++; break;
                case 22: nGamma++; break;
                default: break;
            }
        }
    }
    ret = (nElectron>=1 && nPip>=1 && nPim>=1 && nGamma>=2);
    return ret;
}

