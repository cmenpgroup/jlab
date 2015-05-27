#ifndef OMEGAMIXEDEVENT_H
#define OMEGAMIXEDEVENT_H
#include <vector>
#include <string>
#include "TLorentzVector.h"

using namespace std;

class OmegaMixedEvent
{
    int nEvtToMix;
    int EvtOffset;
    vector<string> Label;
    vector<int> Index;
    vector<string> evtLabel;
    vector<int> evtIndex;
    TLorentzVector Photon1[2];
    TLorentzVector Photon2[2];
    TLorentzVector PiPlus[2];
    TLorentzVector PiMinus[2];
    TLorentzVector Pi0[2];
    TLorentzVector Omega[2];
public:
    OmegaMixedEvent();
    int Get_NumberOfEventsToMix(){return nEvtToMix;};
    void Put_NumberOfEventsToMix(int i);
    int Get_OffsetOfEventsToMix(){return EvtOffset;};
    void Put_OffsetOfEventsToMix(int i);
    bool Check_evtIndex(int iME);
    bool Check_Index(int iMethod);
    int Get_nLabel() {return Label.size();};
    int Get_nIndex() {return Index.size();};
    string Get_Label(int num) {return Label[num];};
    int Get_nEvtLabel() {return evtLabel.size();};
    int Get_nEvtIndex() {return evtIndex.size();};
    string Get_evtLabel(int num) {return evtLabel[num];};
    TLorentzVector Get_Photon1(int iME);
    TLorentzVector Get_Photon2(int iME);
    TLorentzVector Get_PiPlus(int iME);
    TLorentzVector Get_PiMinus(int iME);
    TLorentzVector Get_Pi0(int iME);
    TLorentzVector Get_Omega(int iME);
    void Reconstruct_Pi0(int iME);
    void Reconstruct_Omega(int iME);
    void Mix_Omega(int iMethod);
    void Put_Photon1(TLorentzVector V, int iME);
    void Put_Photon2(TLorentzVector V, int iME);
    void Put_PiPlus(TLorentzVector V, int iME);
    void Put_PiMinus(TLorentzVector V, int iME);
    void Put_Pi0(TLorentzVector V, int iME);
    void Put_Omega(TLorentzVector V, int iME);
    void Clear_TLorentzVectors();
    void Print_Info();
};
#endif
