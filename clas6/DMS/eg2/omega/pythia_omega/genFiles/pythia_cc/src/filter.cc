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

    this->SetKScut(1);
    
    partType.pushback(11);
    partType.pushback(211);
    partType.pushback(-211);
    partType.pushback(22);

    partQty.pushback(1);
    partQty.pushback(1);
    partQty.pushback(1);
    partQty.pushback(2);
    
    cout<<"*********************************************"<<endl;
    cout<<"Initialize the event filter"<<endl;
    cout<<"Particle list:"<<endl;
    cout<<"Type\t Quantity"<<endl;
    for(int i=0; partType.size(); i++){
        cout<<partType[i]<<"\t"<<partQty[i]<<endl;
    }
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
        if(trk.ks[i]==this->GetKScut()){
            switch(trk.type[i]){
                case this->GetPartType(0): nElectron++; break;
                case this->GetPartType(1): nPip++; break;
                case this->GetPartType(2): nPim++; break;
                case this->GetPartType(3): nGamma++; break;
                default: break;
            }
        }
    }
    ret = (nElectron>=this->GetPartQty(0) && nPip>=this->GetPartQty(1) && nPim>=this->GetPartQty(2) && nGamma>=this->GetPartQty(3));
    return ret;
}

void Filter::SetKScut(int ks)
{
    this->KScut = ks;
}

