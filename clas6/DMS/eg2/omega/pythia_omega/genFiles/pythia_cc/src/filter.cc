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
    
    unsigned int n = partType.size();
    partCtr.assign(n,0);
    
    if(this->CheckPartSize()){ // check that the particle lists have the same sizes
        cout<<"*********************************************"<<endl;
        cout<<"Initialize the event filter"<<endl;
        cout<<"Particle list:"<<endl;
        cout<<"Type\t Quantity"<<endl;
    
        for(unsigned int i=0; i<partType.size(); i++){
            cout<<partType[i]<<"\t"<<partQty[i]<<endl;
        }
        cout<<"*********************************************"<<endl;
    }else{
        cout<<"Filter::init, Mismatch in partType and partQty vectors"<<endl;
        exit(0);
    }
}

bool Filter::CheckPartSize()
{
    bool ret = (this->Get_nPartType()==this->Get_nPartQty()) ? true : false;
    return ret;
}

bool Filter::Cut()
{
    int i, j;
    bool ret = true;
    
    if(!this->CheckPartSize()){ // check that the particle lists have the same sizes
        cout<<"Filter::Cut, Mismatch in partType and partQty vectors"<<endl;
        exit(0);
    }

    this->ZeroPartCtr();  // zero each element of the vector for the particle counter
    
    //loop over tracks
    for(i=0; i<trk.Ntracks; i++){
        if(trk.ks[i]==this->GetKScut()){
            for(j=0; j<this->Get_nPartType(); j++){
                if(trk.type[i]==this->GetPartType(j)) partCtr[j]++;
            }
        }
    }
    
    cout<<"Filter::Cut, Start"<<endl;
    for(j=0; j<this->Get_nPartQty(); j++){
        ret = (this->GetPartCtr(j)>=this->GetPartQty(j));
        cout<<"Filter::Cut, "<<ret<<"\t"<<this->GetPartCtr(j)<<"\t"<<this->GetPartQty(j)<<endl;
        if(!ret) return ret;
    }
    
    return ret;
}

void Filter::SetKScut(int ks)
{
    this->KScut = ks;
}

void Filter::ZeroPartCtr()
{
    unsigned int n = this->Get_nPartCtr();
    partCtr.assign(n,0);
}
