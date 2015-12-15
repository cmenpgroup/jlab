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
    
    partType.push_back(11);
    partType.push_back(211);
    partType.push_back(-211);
    partType.push_back(22);

    partQty.push_back(1);
    partQty.push_back(1);
    partQty.push_back(1);
    partQty.push_back(2);
    
    cout<<"Size "<<partType.size()<<"\t"<<partQty.size()<<endl;
    
    cout<<"*********************************************"<<endl;
    cout<<"Initialize the event filter"<<endl;
    cout<<"Particle list:"<<endl;
    cout<<"Type\t Quantity"<<endl;
    for(int i=0; i<partType.size(); i++){
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
            if(trk.type[i]==this->GetPartType(0)) nElectron++;
            if(trk.type[i]==this->GetPartType(1)) nPip++;
            if(trk.type[i]==this->GetPartType(2)) nPim++;
            if(trk.type[i]==this->GetPartType(3)) nGamma++;
        }
    }
    ret = (nElectron>=this->GetPartQty(0) && nPip>=this->GetPartQty(1) && nPim>=this->GetPartQty(2) && nGamma>=this->GetPartQty(3));
    return ret;
}

void Filter::SetKScut(int ks)
{
    this->KScut = ks;
}

