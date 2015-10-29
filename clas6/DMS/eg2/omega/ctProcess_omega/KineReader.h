#ifndef KINEREADER_H
#define KINEREADER_H
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

using namespace std;

class KineReader
{
    TBranch *b_kine;
    TLeaf *l_evtnum;
    TLeaf *l_elecverttarg;
    TLeaf *l_q2;
    TLeaf *l_nu;
    TLeaf *l_xb;
    TLeaf *l_w;
    TLeaf *l_xcorr;
    TLeaf *l_ycorr;
    TLeaf *l_zcorr;
    TLeaf *l_nelec;
    TLeaf *l_npip;
    TLeaf *l_npim;
    TLeaf *l_ngam;
    
public:
    KineReader(TTree *tree);
    void ReadEntry(int num);
    float Get_EvtNum() {return l_evtnum->GetValue();};
    float Get_ElecVertTarg() {return l_elecverttarg->GetValue();};
    float Get_Q2() {return l_q2->GetValue();};
    float Get_Nu() {return l_nu->GetValue();};
    float Get_Xb() {return l_xb->GetValue();};
    float Get_W() {return l_w->GetValue();};
    float Get_Xcorr() {return l_xcorr->GetValue();};
    float Get_Ycorr() {return l_ycorr->GetValue();};
    float Get_Zcorr() {return l_zcorr->GetValue();};
    int Get_nElec() {return l_nelec->GetValue();};
    int Get_nPip() {return l_npip->GetValue();};
    int Get_nPim() {return l_npim->GetValue();};
    int Get_nGam() {return l_ngam->GetValue();};
};
#endif
