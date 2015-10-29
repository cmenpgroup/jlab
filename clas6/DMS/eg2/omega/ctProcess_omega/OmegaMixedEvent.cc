#include <vector>
#include <string>
#include "OmegaMixedEvent.h"
#include <iostream>
#include "TLorentzVector.h"
OmegaMixedEvent::OmegaMixedEvent()
{
    int i;
    
    Label.push_back("Mixing Photon 1");
    Label.push_back("Mixing Photon 2");
    Label.push_back("Mixing pi+");
    Label.push_back("Mixing pi-");
    Label.push_back("Mixing pi0");
    
    Index.push_back(0);
    Index.push_back(1);
    Index.push_back(2);
    Index.push_back(3);
    Index.push_back(4);

    evtLabel.push_back("In Time Event");
    evtLabel.push_back("Out of Time Event");
    
    evtIndex.push_back(0);
    evtIndex.push_back(1);

    for(i=0; i<evtIndex.size(); i++){
        Photon1[i].SetPxPyPzE(0.,0.,0.,0.);
        Photon2[i].SetPxPyPzE(0.,0.,0.,0.);
        Pi0[i].SetPxPyPzE(0.,0.,0.,0.);
        PiPlus[i].SetPxPyPzE(0.,0.,0.,0.);
        PiMinus[i].SetPxPyPzE(0.,0.,0.,0.);
        Omega[i].SetPxPyPzE(0.,0.,0.,0.);
    }
}

void OmegaMixedEvent::Clear_TLorentzVectors()
{
    int i;
    
    for(i=0; i < this->Get_nEvtIndex(); i++){
        this->Photon1[i].SetPxPyPzE(0.,0.,0.,0.);
        this->Photon2[i].SetPxPyPzE(0.,0.,0.,0.);
        this->Pi0[i].SetPxPyPzE(0.,0.,0.,0.);
        this->PiPlus[i].SetPxPyPzE(0.,0.,0.,0.);
        this->PiMinus[i].SetPxPyPzE(0.,0.,0.,0.);
        this->Omega[i].SetPxPyPzE(0.,0.,0.,0.);
    }
}

void OmegaMixedEvent::Put_NumberOfEventsToMix(int i){
    this->nEvtToMix = i;
}

void OmegaMixedEvent::Put_OffsetOfEventsToMix(int i){
    this->EvtOffset = i;
}

// check that the event index is in the correct range
bool OmegaMixedEvent::Check_evtIndex(int iME)
{
    bool ret = (iME>=0 && iME<this->Get_nEvtIndex()) ? true : false;
    
    return ret;
}

// check that the method index is in the correct range
bool OmegaMixedEvent::Check_Index(int iMethod)
{
    bool ret = (iMethod>=0 && iMethod<this->Get_nIndex()) ? true : false;
    
    return ret;
}

void OmegaMixedEvent::Put_Photon1(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->Photon1[iME] = V;
}

void OmegaMixedEvent::Put_Photon2(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->Photon2[iME] = V;
}

void OmegaMixedEvent::Put_PiPlus(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->PiPlus[iME] = V;
}

void OmegaMixedEvent::Put_PiMinus(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->PiMinus[iME] = V;
}

void OmegaMixedEvent::Put_Pi0(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->Pi0[iME] = V;
}

void OmegaMixedEvent::Put_Omega(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->Omega[iME] = V;
}

TLorentzVector OmegaMixedEvent::Get_Photon1(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
  
    if(this->Check_evtIndex(iME)) ret = Photon1[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_Photon2(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = Photon2[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_PiPlus(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = PiPlus[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_PiMinus(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = PiMinus[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_Pi0(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = Pi0[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_Omega(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = Omega[iME];
    
    return ret;
}

void OmegaMixedEvent::Reconstruct_Pi0(int iME){
    if(this->Check_evtIndex(iME)){
        TLorentzVector tempPhoton1 = this->Get_Photon1(iME);
        TLorentzVector tempPhoton2 = this->Get_Photon2(iME);
        
        this->Put_Pi0(tempPhoton1 + tempPhoton2, iME);
    }
}

void OmegaMixedEvent::Reconstruct_Omega(int iME){
    if(this->Check_evtIndex(iME)){
        TLorentzVector tempPhoton1 = this->Get_Photon1(iME);
        TLorentzVector tempPhoton2 = this->Get_Photon2(iME);
        TLorentzVector tempPiPlus = this->Get_PiPlus(iME);
        TLorentzVector tempPiMinus = this->Get_PiMinus(iME);
        
        this->Put_Omega(tempPhoton1 + tempPhoton2 + tempPiPlus + tempPiMinus, iME);
    }
}

void OmegaMixedEvent::Mix_Omega(int iMethod){

    double theta, phi, tempPx, tempPy, tempPz, tempE;
    TLorentzVector tempPhoton1 = this->Get_Photon1(0);
    TLorentzVector tempPhoton2 = this->Get_Photon2(0);
    TLorentzVector tempPiPlus = this->Get_PiPlus(0);
    TLorentzVector tempPiMinus = this->Get_PiMinus(0);
    TLorentzVector tempPi0 = tempPhoton1 + tempPhoton2;
    TLorentzVector A;

    switch (iMethod){
        case 0:
            tempPhoton1 = this->Get_Photon1(1);
            tempPi0 = tempPhoton1 + tempPhoton2;
            break;
        case 1:
            tempPhoton2 = this->Get_Photon2(1);
            tempPi0 = tempPhoton1 + tempPhoton2;
            break;
        case 2:
            A = this->Get_PiPlus(1);
            A.SetTheta(tempPiPlus.Theta());
            A.SetPhi(tempPiPlus.Phi());
            tempPiPlus = A;
            break;
        case 3:
            A = this->Get_PiMinus(1);
            A.SetTheta(tempPiMinus.Theta());
            A.SetPhi(tempPiMinus.Phi());
            tempPiMinus = A;
            break;
        case 4:
            A = this->Get_Photon1(1) + this->Get_Photon2(1);
            A.SetTheta(tempPi0.Theta());
            A.SetPhi(tempPi0.Phi());
            tempPi0 = A;
            break;
        default:
            cout << "OmegaMixedEvent::Mix_Omega - incorrect iMethod " << iMethod <<endl;
            exit(0);
            break;
    }
    this->Put_Pi0(tempPi0, 1);
    this->Put_Omega((tempPi0 + tempPiPlus + tempPiMinus), 1);
}

void OmegaMixedEvent::Print_Info()
{
    int ii;
    cout<<"Mixed Event Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nLabel();ii++){
        cout << "Method " << ii+1 << "\t" << this->Get_Label(ii) << endl;
    }
    for(ii=0;ii<this->Get_nEvtLabel();ii++){
        cout << "Event type " << ii+1 << "\t" << this->Get_evtLabel(ii) << endl;
    }
    cout << endl;
}
