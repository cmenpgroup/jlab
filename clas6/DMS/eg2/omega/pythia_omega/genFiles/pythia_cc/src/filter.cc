#include "filter.hh"
#include "funcs.hh"
#include "tracks.hh"
#include "read_optFile.hh"

Filter::Filter()
{
  init();
}

Filter::~Filter()
{;}

void Filter::init()
{

    this->SetKScut(1);
    
    partType = rof->fPartList; // list of particle for event selection
    partQty = rof->fPartQty; // list of particle quantities for event selection
    
    unsigned int n = partType.size();
    partCtr.assign(n,0); // zero the particle counters
    
    if(this->CheckPartSize()){ // check that the particle lists have the same sizes
        cout<<"*********************************************"<<endl;
        cout<<"Initialize the event filter"<<endl;
        cout<<"Particle list:"<<endl;
        cout<<"--------------"<<endl;
        cout<<"Type \t Quantity"<<endl;
    
        for(unsigned int i=0; i<partType.size(); i++){
            cout<<partType[i]<<" \t "<<partQty[i]<<endl;
        }
        cout<<"*********************************************"<<endl;
    }else{
        cout<<"Filter::init, Mismatch in partType and partQty vectors"<<endl;
        exit(0);
    }
}

// check that the particle typ and paricle quantity have the same size vectors
bool Filter::CheckPartSize()
{
    bool ret = (this->Get_nPartType()==this->Get_nPartQty()) ? true : false;
    return ret;
}

// check that the event has the particle final state topology
bool Filter::Cut()
{
    int i, j;
    bool ret = true;
    
    if(this->Get_nPartType()==0) return ret; // skip if no particles selected
        
    if(!this->CheckPartSize()){ // check that the particle lists have the same sizes
        cout<<"Filter::Cut, Mismatch in partType and partQty vectors"<<endl;
        exit(0);
    }

    this->ZeroPartCtr();  // zero each element of the vector for the particle counter
    
    //loop over tracks and count up the number of final state particles
    for(i=0; i<trk.Ntracks; i++){
        if(trk.ks[i]==this->GetKScut()){ // check for particle state
            for(j=0; j<this->Get_nPartType(); j++){  // loop over PID
                if(trk.type[i]==this->GetPartType(j)) partCtr[j]++;
            }
        }
    }
    
    // loop over the particle list
    for(j=0; j<this->Get_nPartQty(); j++){
        ret = (this->GetPartCtr(j)>=this->GetPartQty(j)); // check paricle tally against user-defined topology
        if(!ret) return ret;  // if match is false, exit
    }
    
    return ret;
}

// set the pythia KS parameter
void Filter::SetKScut(int ks)
{
    this->KScut = ks;
}

// set the particle counter to zero
void Filter::ZeroPartCtr()
{
    unsigned int n = this->Get_nPartCtr();
    partCtr.assign(n,0);
}
