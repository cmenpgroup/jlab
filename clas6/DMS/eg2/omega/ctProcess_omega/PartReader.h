#ifndef PARTREADER_H
#define PARTREADER_H
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLorentzVector.h"

using namespace std;

class PartReader
{
    TBranch *b_part;
    TLeaf *l_sector;
    TLeaf *l_charge;
    TLeaf *l_pid;
    TLeaf *l_beta;
    TLeaf *l_px;
    TLeaf *l_py;
    TLeaf *l_pz;
    TLeaf *l_mom;
    TLeaf *l_mass2;
    TLeaf *l_x;
    TLeaf *l_y;
    TLeaf *l_z;
    TLeaf *l_ecx;
    TLeaf *l_ecy;
    TLeaf *l_ecz;
    TLeaf *l_ecu;
    TLeaf *l_ecv;
    TLeaf *l_ecw;
    TLeaf *l_ectot;
    TLeaf *l_ecin;
    TLeaf *l_ecout;
    TLeaf *l_ectime;
    TLeaf *l_ecpath;
    TLeaf *l_echit_m2;
    TLeaf *l_echit_m3;
    TLeaf *l_echit_m4;
    TLeaf *l_chi2ec;
    TLeaf *l_sctime;
    TLeaf *l_scpath;
    TLeaf *l_ccnphe;
    TLeaf *l_t;
    TLeaf *l_xf;
    TLeaf *l_mx2;
    TLeaf *l_pt;
    TLeaf *l_zh;
    TLeaf *l_thetapq;
    TLeaf *l_phipq;
    TLeaf *l_timecorr4;
public:
    PartReader(TTree *tree, string branchName);
    void ReadEntry(int num);
    int Get_Sector() {return l_sector->GetValue();};
    float Get_Charge() {return l_charge->GetValue();};
    float Get_Pid() {return l_pid->GetValue();};
    float Get_Beta() {return l_beta->GetValue();};
    float Get_Px() {return l_px->GetValue();};
    float Get_Py() {return l_py->GetValue();};
    float Get_Pz() {return l_pz->GetValue();};
    float Get_Mom() {return l_mom->GetValue();};
    float Get_Mass2() {return l_mass2->GetValue();};
    float Get_X() {return l_x->GetValue();};
    float Get_Y() {return l_y->GetValue();};
    float Get_Z() {return l_z->GetValue();};
    float Get_ECx() {return l_ecx->GetValue();};
    float Get_ECy() {return l_ecy->GetValue();};
    float Get_ECz() {return l_ecz->GetValue();};
    float Get_ECu() {return l_ecu->GetValue();};
    float Get_ECv() {return l_ecv->GetValue();};
    float Get_ECw() {return l_ecw->GetValue();};
    float Get_ECtot() {return l_ectot->GetValue();};
    float Get_ECin() {return l_ecin->GetValue();};
    float Get_ECout() {return l_ecout->GetValue();};
    float Get_ECtime() {return l_ectime->GetValue();};
    float Get_ECpath() {return l_ecpath->GetValue();};
    float Get_EChit_M2() {return l_echit_m2->GetValue();};
    float Get_EChit_M3() {return l_echit_m3->GetValue();};
    float Get_EChit_M4() {return l_echit_m4->GetValue();};
    float Get_Chi2EC() {return l_chi2ec->GetValue();};
    float Get_SCtime() {return l_sctime->GetValue();};
    float Get_SCpath() {return l_scpath->GetValue();};
    float Get_CCnphe() {return l_ccnphe->GetValue();};
    float Get_T() {return l_t->GetValue();};
    float Get_Xf() {return l_xf->GetValue();};
    float Get_Mx2() {return l_mx2->GetValue();};
    float Get_Pt() {return l_pt->GetValue();};
    float Get_Zh() {return l_zh->GetValue();};
    float Get_ThetPQ() {return l_thetapq->GetValue();};
    float Get_PhiPQ() {return l_phipq->GetValue();};
    float Get_TimeCorr4() {return l_timecorr4->GetValue();};

    double Get_TOF_MassSquared();
    double Get_BetaFromMass(double fMass);
    double Get_BetaDifference(double fMass);
    TLorentzVector GetLorentzVector(double mass);
    TVector3 GetVertex();
};
#endif
