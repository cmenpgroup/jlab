#ifndef HISTMANAGER_H
#define HISTMANAGER_H
#include <vector>
#include <string>
#include "DetectedParticles.h"
#include "ParticleList.h"
#include "EG2Target.h"
#include "OmegaMixedEvent.h"
#include "ElectronID.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

class HistManager
{
private:
    int MAX_SECTORS; // max. number of CLAS sectors

    int ID_ELECTRON; // PDG electron id
    int ID_PHOTON;  // PDG photon id
    int ID_PION_POS;  // PDG pi+ id
    int ID_PION_NEG;  // PDG pi- id
    int ID_PROTON; // PDG proton id

    double MASS_PHOTON; // mass of photon in GeV/c^2
    double MASS_ELECTRON; // mass of charged pion in GeV/c^2
    double MASS_PION_CHARGED; // mass of charged pion in GeV/c^2
    double MASS_PION_NEUTRAL; // mass of neutral pion in GeV/c^2
    double MASS_PROTON; // mass of proton in GeV/c^2
    double MASS_DEUTERIUM; // mass of deuterium in GeV/c^2

    double LIGHTSPEED; // speed of light in cm/ns

    double BEAM_ENERGY; // electron beam energy in GeV

    TH1D *q2;
    TH2D *NumDetPart;
    TH2D *q2_VS_theta;
    TH1D *StartTime;
    TH1D *nu_EnergyTransfer;
    TH1D *elecZVert;
    TH1D *OpAng_2Photons;
    TH1D *OpAng_elecPhoton1;
    TH1D *OpAng_elecPhoton2;
    TH2D *Xvert;
    TH2D *Yvert;
    TH2D *ZVertDiff;
    TH2D *Xvert_VS_Yvert[5];
    TH2D *Beta_VS_Momentum;
    TH2D *Beta_VS_Momentum_Recalc;
    TH2D *dBeta_VS_Momentum[5];
    TH2D *dBeta_VS_Momentum_EPC[5];
    TH2D *Beta_Recalc;
    TH2D *scMassSquared_NC;
    TH2D *scMassSquared_EC;
    TH2D *scMassSquared_PC;
    TH2D *scMassSquared_EPC;
    TH2D *Theta_VS_Phi[7];
    TH2D *TotalMomentum;
    TH2D *elecZVert_VS_Phi;
    TH2D *elecZVert_VS_Phi_Corr;

    TH2D *CCnphe;
    TH2D *ECu;
    TH2D *ECv;
    TH2D *ECw;
    TH2D *EChit_M2;
    TH2D *EChit_M3;
    TH2D *EChit_M4;
    TH2D *EChit_M2_cuts;
    TH2D *EChit_M3_cuts;
    TH2D *EChit_M4_cuts;
    TH2D *dtime_ECSC;
    TH2D *ECtot_VS_P[5];
    TH2D *ECtotP_VS_P[5];
    TH2D *ECin_VS_ECout[5];
    TH2D *EChit_M2_VS_scMsq[5];
    TH2D *EChit_M3_VS_scMsq[5];
    TH2D *EChit_M4_VS_scMsq[5];
    
    TH2D *Mom_elecID;
    TH2D *CCnphe_elecID;
    TH2D *ECu_elecID;
    TH2D *ECv_elecID;
    TH2D *ECw_elecID;
    TH2D *dtime_ECSC_elecID;
    TH2D *ECtot_VS_P_elecID[9];
    TH2D *ECtotP_VS_P_elecID[9];
    TH2D *ECin_VS_ECout_elecID[9];
    TH2D *Mom_VS_ECout_elecID[9];
    TH2D *ECu_VS_ECout_elecID[9];
    TH2D *ECv_VS_ECout_elecID[9];
    TH2D *ECw_VS_ECout_elecID[9];

    TH1D *q2_ECoutCut;
    TH1D *elecZVert_ECoutCut;
    TH1D *ECtotMinusECin_ECoutCut;
    TH2D *Beta_VS_Momentum_ECoutCut;
    TH2D *Theta_VS_Phi_ECoutCut;
    TH2D *ECtotP_VS_P_ECoutCut;
    TH2D *ECtot_VS_P_ECoutCut;

    TH1D *q2_AntiECoutCut;
    TH1D *elecZVert_AntiECoutCut;
    TH1D *ECtotMinusECin_AntiECoutCut;
    TH2D *Beta_VS_Momentum_AntiECoutCut;
    TH2D *Theta_VS_Phi_AntiECoutCut;
    TH2D *ECtotP_VS_P_AntiECoutCut;
    TH2D *ECtot_VS_P_AntiECoutCut;

    TH2D *ECinP_VS_ECoutP[6];
    TH2D *ECinP_VS_ECoutP_cut[6];
    TH2D *ECinP_VS_ECoutP_Range[5];

    TH2D *ECin_VS_ECout_ECfid;
    TH2D *ECin_VS_ECout_elecID_All;
    TH2D *ECtotP_VS_P_Sector[6];
    TH2D *ECtotP_VS_P_ECPCut[6];
    TH2D *EC_XvsY_local_Sector[6];
    TH2D *EC_XvsY_local_ECoutCut[6];
    TH2D *EC_XvsY_local_FidCut[6];
    TH2D *EC_XvsY_local_AntiFidCut[6];

    TH1D *hW[3];
    TH1D *hMx[3];
    TH1D *z_fracE[3];
    TH1D *LongMom[3];
    TH1D *TransMom[3];
    TH1D *MissMom[3];
    TH1D *MMsq[3];
    TH1D *PtSq_Omega_AllCuts[3];
    TH1D *PtSq_Omega_AllCuts_IMOmegaCut[3];
    TH1D *PtSq_Omega_AllCuts_IMOmegaSBCut[3];
    TH2D *elecZVertSector;
    TH2D *elecZVertSector_Corr;
    TH2D *OpAng_VS_IM2Photons[3];
    TH2D *OpAng_VS_E[3];
    TH2D *OpAng_VS_E_MassPi0Cut[3];
    TH2D *IM2Pions_VS_IMOmega[3];
    TH2D *IM2Pions_VS_IMOmega_AllCuts[3];
    TH2D *IM2Photons_VS_IMOmega[3];
    TH2D *W_VS_IMOmega_AllCuts[3];
    TH2D *Q2_VS_IMOmega[3];
    TH2D *Pt_VS_IMOmega[3];
    TH2D *Pl_VS_IMOmega[3];
    TH2D *OpAng_VS_IMOmega[3];
    TH2D *IMOmega[3];
    TH2D *IMOmega_woCut[3];
    TH2D *IMOmega_antiCut[3];
    TH2D *IM2Photons[3];
    TH2D *IM2Photons_woCut[3];
    TH2D *Xvert_VS_Yvert_AllCuts[3];
    TH2D *Xvert_VS_Yvert_Omega[3];

    TH2D *IM2Photons_ME[3];
    TH2D *IM2Photons_OpAng_ElecPhoton_Cut_ME[3];
    TH2D *IMOmega_ME[3];
    TH2D *IMOmega_OpAng_ElecPhoton_Cut_ME[3];
    TH2D *IMOmega_MassPi0Cut_ME[3];
    TH2D *IMOmega_ZVertCut_ME[3];
    TH2D *IMOmega_QsqCut_ME[3];
    TH2D *IMOmega_AllCuts_ME[3];
    
    TH2D *VirtualPhotonAngle_VS_IMOmega_AllCuts[3];
    TH2D *OpAngVPomega_VS_IMOmega_AllCuts[3];
    TH2D *Pt_VS_IMOmega_AllCuts[3];
    TH2D *Pl_VS_Pt_AllCuts[3];
    TH2D *OpAng_VS_IMOmega_AllCuts[3];
    TH2D *OpAngPairs_VS_IMOmega_AllCuts[3];
    
    TH2D *RelativityOpAngPhotonsA;
    TH2D *RelativityOpAngPhotonsB;
    TH1D *GammaPi0;
    TH1D *BetaPi0;

    // photon id
    TH1D *MomentumPhoton1;
    TH1D *MomentumPhoton2;
    TH1D *MomentumPhoton1_cut;
    TH1D *MomentumPhoton2_cut;
    TH1D *BetaPhoton1;
    TH1D *BetaPhoton2;
    TH1D *BetaPhoton1_cut;
    TH1D *BetaPhoton2_cut;
    TH1D *ECuPhoton1;
    TH1D *ECuPhoton2;
    TH1D *ECvPhoton1;
    TH1D *ECvPhoton2;
    TH1D *ECwPhoton1;
    TH1D *ECwPhoton2;
    TH1D *ECuPhoton1_cut;
    TH1D *ECuPhoton2_cut;
    TH1D *ECvPhoton1_cut;
    TH1D *ECvPhoton2_cut;
    TH1D *ECwPhoton1_cut;
    TH1D *ECwPhoton2_cut;
    TH1D *ECtime_ECl_Start_Photon1;
    TH1D *ECtime_ECl_Start_Photon2;
    TH1D *ECtime_ECl_Photon1;
    TH1D *ECtime_ECl_Photon2;
    TH1D *ECtime_ECl_Photon1_cut;
    TH1D *ECtime_ECl_Photon2_cut;
    TH1D *ECtimePhoton1;
    TH1D *ECtimePhoton2;
    TH1D *ECpathPhoton1;
    TH1D *ECpathPhoton2;
    TH1D *ECpathtimePhoton1;
    TH1D *ECpathtimePhoton2;
    TH2D *ECtotP_vs_P_Photon1;
    TH2D *ECtotP_vs_P_Photon2;
    TH2D *ECin_vs_ECout_Photon1;
    TH2D *ECin_vs_ECout_Photon2;
    TH2D *ECtotP_vs_P_InOutZeroCut_Photon1;
    TH2D *ECtotP_vs_P_InOutZeroCut_Photon2;
    TH2D *ECin_vs_ECout_InOutZeroCut_Photon1;
    TH2D *ECin_vs_ECout_InOutZeroCut_Photon2;
    TH2D *EC_XvsY_local_Sector_Photon1[6];
    TH2D *EC_XvsY_local_FidCut_Photon1[6];
    TH2D *EC_XvsY_local_AntiFidCut_Photon1[6];
    TH2D *EC_XvsY_local_Sector_Photon2[6];
    TH2D *EC_XvsY_local_FidCut_Photon2[6];
    TH2D *EC_XvsY_local_AntiFidCut_Photon2[6];

    //recon
    TH2D *mass2Pions_VS_massOmega_NC[3];
    TH2D *mass2Pions_VS_massOmega_EPC[3];
    TH2D *mass2Pions_VS_massOmega_EPOC[3];
    
public:
    HistManager();
    void BookHist();
    void WriteHist(string RootFile);
    TH1D* GetQ2() { return q2; };
    TH2D* GetNumDetPart() { return NumDetPart; };
    TH2D* GetQ2_VS_theta() { return q2_VS_theta; };
    TH1D* GetStartTime() { return StartTime; };
    TH1D* GetNu_EnergyTransfer() { return nu_EnergyTransfer; };
    TH1D* GetElecZVert() { return elecZVert; };
    TH1D* GetOpAng_2Photons() { return OpAng_2Photons; };
    TH1D* GetOpAng_elecPhoton1() { return OpAng_elecPhoton1; };
    TH1D* GetOpAng_elecPhoton2() { return OpAng_elecPhoton2; };
    TH2D* GetXvert() { return Xvert; };
    TH2D* GetYvert() { return Yvert; };
    TH2D* GetZVertDiff() { return ZVertDiff; };
    TH2D* GetXvert_VS_Yvert(int index) { return Xvert_VS_Yvert[index]; };
    TH2D* GetBeta_VS_Momentum() { return Beta_VS_Momentum; };
    TH2D* GetBeta_VS_Momentum_Recalc() { return Beta_VS_Momentum_Recalc; };
    TH2D* GetDBeta_VS_Momentum(int index) { return dBeta_VS_Momentum[index]; };
    TH2D* GetDBeta_VS_Momentum_EPC(int index) { return dBeta_VS_Momentum_EPC[index]; };
    TH2D* GetBeta_Recalc() { return Beta_Recalc; };
    TH2D* GetScMassSquared_NC() { return scMassSquared_NC; };
    TH2D* GetScMassSquared_EC() { return scMassSquared_EC; };
    TH2D* GetScMassSquared_PC() { return scMassSquared_PC; };
    TH2D* GetScMassSquared_EPC() { return scMassSquared_EPC; };
    TH2D* GetTheta_VS_Phi(int index) { return Theta_VS_Phi[index]; };
    TH2D* GetTotalMomentum() { return TotalMomentum; };
    TH2D* GetElecZVert_VS_Phi() { return elecZVert_VS_Phi; };
    TH2D* GetElecZVert_VS_Phi_Corr() { return elecZVert_VS_Phi_Corr; };

    TH2D* GetCCnphe() { return CCnphe; };
    TH2D* GetECu() { return ECu; };
    TH2D* GetECv() { return ECv; };
    TH2D* GetECw() { return ECw; };
    TH2D* GetEChit_M2() { return EChit_M2; };
    TH2D* GetEChit_M3() { return EChit_M3; };
    TH2D* GetEChit_M4() { return EChit_M4; };
    TH2D* GetEChit_M2_cuts() { return EChit_M2_cuts; };
    TH2D* GetEChit_M3_cuts() { return EChit_M3_cuts; };
    TH2D* GetEChit_M4_cuts() { return EChit_M4_cuts; };
    TH2D* GetDtime_ECSC() { return dtime_ECSC; };
    TH2D* GetECtot_VS_P(int index) { return ECtot_VS_P[index]; };
    TH2D* GetECtotP_VS_P(int index) { return ECtotP_VS_P[index]; };
    TH2D* GetECin_VS_ECout(int index) { return ECin_VS_ECout[index]; };
    TH2D* GetEChit_M2_VS_scMsq(int index) { return EChit_M2_VS_scMsq[index]; };
    TH2D* GetEChit_M3_VS_scMsq(int index) { return EChit_M3_VS_scMsq[index]; };
    TH2D* GetEChit_M4_VS_scMsq(int index) { return EChit_M4_VS_scMsq[index]; };
    
    TH2D* GetMom_elecID() { return Mom_elecID; };
    TH2D* GetCCnphe_elecID() { return CCnphe_elecID; };
    TH2D* GetECu_elecID() { return ECu_elecID; };
    TH2D* GetECv_elecID() { return ECv_elecID; };
    TH2D* GetECw_elecID() { return ECw_elecID; };
    TH2D* GetDtime_ECSC_elecID() { return dtime_ECSC_elecID; };
    TH2D* GetECtot_VS_P_elecID(int index) { return ECtot_VS_P_elecID[index]; };
    TH2D* GetECtotP_VS_P_elecID(int index) { return ECtotP_VS_P_elecID[index]; };
    TH2D* GetECin_VS_ECout_elecID(int index) { return ECin_VS_ECout_elecID[index]; };
    TH2D* GetMom_VS_ECout_elecID(int index) { return Mom_VS_ECout_elecID[index]; };
    TH2D* GetECu_VS_ECout_elecID(int index) { return ECu_VS_ECout_elecID[index]; };
    TH2D* GetECv_VS_ECout_elecID(int index) { return ECv_VS_ECout_elecID[index]; };
    TH2D* GetECw_VS_ECout_elecID(int index) { return ECw_VS_ECout_elecID[index]; };

    TH1D* GetQ2_ECoutCut() { return q2_ECoutCut; };
    TH1D* GetElecZVert_ECoutCut() { return elecZVert_ECoutCut; };
    TH1D* GetECtotMinusECin_ECoutCut() { return ECtotMinusECin_ECoutCut; };
    TH2D* GetBeta_VS_Momentum_ECoutCut() { return Beta_VS_Momentum_ECoutCut; };
    TH2D* GetTheta_VS_Phi_ECoutCut() { return Theta_VS_Phi_ECoutCut; };
    TH2D* GetECtotP_VS_P_ECoutCut() { return ECtotP_VS_P_ECoutCut; };
    TH2D* GetECtot_VS_P_ECoutCut() { return ECtot_VS_P_ECoutCut; };

    TH1D* GetQ2_AntiECoutCut() { return q2_AntiECoutCut; };
    TH1D* GetElecZVert_AntiECoutCut() { return elecZVert_AntiECoutCut; };
    TH1D* GetECtotMinusECin_AntiECoutCut() { return ECtotMinusECin_AntiECoutCut; };
    TH2D* GetBeta_VS_Momentum_AntiECoutCut() { return Beta_VS_Momentum_AntiECoutCut; };
    TH2D* GetTheta_VS_Phi_AntiECoutCut() { return Theta_VS_Phi_AntiECoutCut; };
    TH2D* GetECtotP_VS_P_AntiECoutCut() { return ECtotP_VS_P_AntiECoutCut; };
    TH2D* GetECtot_VS_P_AntiECoutCut() { return ECtot_VS_P_AntiECoutCut; };

    TH2D* GetECinP_VS_ECoutP(int index) { return ECinP_VS_ECoutP[index]; };
    TH2D* GetECinP_VS_ECoutP_cut(int index) { return ECinP_VS_ECoutP_cut[index]; };
    TH2D* GetECinP_VS_ECoutP_Range(int index) { return ECinP_VS_ECoutP_Range[index]; };

    TH2D* GetECin_VS_ECout_ECfid() { return ECin_VS_ECout_ECfid; };
    TH2D* GetECin_VS_ECout_elecID_All() { return ECin_VS_ECout_elecID_All; };
    TH2D* GetECtotP_VS_P_Sector(int index) { return ECtotP_VS_P_Sector[index]; };
    TH2D* GetECtotP_VS_P_ECPCut(int index) { return ECtotP_VS_P_ECPCut[index]; };
    TH2D* GetEC_XvsY_local_Sector(int index) { return EC_XvsY_local_Sector[index]; };
    TH2D* GetEC_XvsY_local_ECoutCut(int index) { return EC_XvsY_local_ECoutCut[index]; };
    TH2D* GetEC_XvsY_local_FidCut(int index) { return EC_XvsY_local_FidCut[index]; };
    TH2D* GetEC_XvsY_local_AntiFidCut(int index) { return EC_XvsY_local_AntiFidCut[index]; };

    TH1D* GetHW(int index) { return hW[index]; };
    TH1D* GetHMx(int index) { return hMx[index]; };
    TH1D* GetZ_fracE(int index) { return z_fracE[index]; };
    TH1D* GetLongMom(int index) { return LongMom[index]; };
    TH1D* GetTransMom(int index) { return TransMom[index]; };
    TH1D* GetMissMom(int index) { return MissMom[index]; };
    TH1D* GetMMsq(int index) { return MMsq[index]; };
    TH1D* GetPtSq_Omega_AllCuts(int index) { return PtSq_Omega_AllCuts[index]; };
    TH1D* GetPtSq_Omega_AllCuts_IMOmegaCut(int index) { return PtSq_Omega_AllCuts_IMOmegaCut[index]; };
    TH1D* GetPtSq_Omega_AllCuts_IMOmegaSBCut(int index) { return PtSq_Omega_AllCuts_IMOmegaSBCut[index]; };
    TH2D* GetElecZVertSector() { return elecZVertSector; };
    TH2D* GetElecZVertSector_Corr() { return elecZVertSector_Corr; };
    TH2D* GetOpAng_VS_IM2Photons(int index) { return OpAng_VS_IM2Photons[index]; };
    TH2D* GetOpAng_VS_E(int index) { return OpAng_VS_E[index]; };
    TH2D* GetOpAng_VS_E_MassPi0Cut(int index) { return OpAng_VS_E_MassPi0Cut[index]; };
    TH2D* GetIM2Pions_VS_IMOmega(int index) { return IM2Pions_VS_IMOmega[index]; };
    TH2D* GetIM2Pions_VS_IMOmega_AllCuts(int index) { return IM2Pions_VS_IMOmega_AllCuts[index]; };
    TH2D* GetIM2Photons_VS_IMOmega(int index) { return IM2Photons_VS_IMOmega[index]; };
    TH2D* GetW_VS_IMOmega_AllCuts(int index) { return W_VS_IMOmega_AllCuts[index]; };
    TH2D* GetQ2_VS_IMOmega(int index) { return Q2_VS_IMOmega[index]; };
    TH2D* GetPt_VS_IMOmega(int index) { return Pt_VS_IMOmega[index]; };
    TH2D* GetPl_VS_IMOmega(int index) { return Pl_VS_IMOmega[index]; };
    TH2D* GetOpAng_VS_IMOmega(int index) { return OpAng_VS_IMOmega[index]; };
    TH2D* GetIMOmega(int index) { return IMOmega[index]; };
    TH2D* GetIMOmega_woCut(int index) { return IMOmega_woCut[index]; };
    TH2D* GetIMOmega_antiCut(int index) { return IMOmega_antiCut[index]; };
    TH2D* GetIM2Photons(int index) { return IM2Photons[index]; };
    TH2D* GetIM2Photons_woCut(int index) { return IM2Photons_woCut[index]; };
    TH2D* GetXvert_VS_Yvert_AllCuts(int index) { return Xvert_VS_Yvert_AllCuts[index]; };
    TH2D* GetXvert_VS_Yvert_Omega(int index) { return Xvert_VS_Yvert_Omega[index]; };

    TH2D* GetIM2Photons_ME(int index) { return IM2Photons_ME[index]; };
    TH2D* GetIM2Photons_OpAng_ElecPhoton_Cut_ME(int index) { return IM2Photons_OpAng_ElecPhoton_Cut_ME[index]; };
    TH2D* GetIMOmega_ME(int index) { return IMOmega_ME[index]; };
    TH2D* GetIMOmega_OpAng_ElecPhoton_Cut_ME(int index) { return IMOmega_OpAng_ElecPhoton_Cut_ME[index]; };
    TH2D* GetIMOmega_MassPi0Cut_ME(int index) { return IMOmega_MassPi0Cut_ME[index]; };
    TH2D* GetIMOmega_ZVertCut_ME(int index) { return IMOmega_ZVertCut_ME[index]; };
    TH2D* GetIMOmega_QsqCut_ME(int index) { return IMOmega_QsqCut_ME[index]; };
    TH2D* GetIMOmega_AllCuts_ME(int index) { return IMOmega_AllCuts_ME[index]; };

    TH2D* GetVirtualPhotonAngle_VS_IMOmega_AllCuts(int index) { return VirtualPhotonAngle_VS_IMOmega_AllCuts[index]; };
    TH2D* GetOpAngVPomega_VS_IMOmega_AllCuts(int index) { return OpAngVPomega_VS_IMOmega_AllCuts[index]; };
    TH2D* GetPt_VS_IMOmega_AllCuts(int index) { return Pt_VS_IMOmega_AllCuts[index]; };
    TH2D* GetPl_VS_Pt_AllCuts(int index) { return Pl_VS_Pt_AllCuts[index]; };
    TH2D* GetOpAng_VS_IMOmega_AllCuts(int index) { return OpAng_VS_IMOmega_AllCuts[index]; };
    TH2D* GetOpAngPairs_VS_IMOmega_AllCuts(int index) { return OpAngPairs_VS_IMOmega_AllCuts[index]; };
    
    TH2D* GetRelativityOpAngPhotonsA() { return RelativityOpAngPhotonsA; };
    TH2D* GetRelativityOpAngPhotonsB() { return RelativityOpAngPhotonsB; };
    TH1D* GetGammaPi0() { return GammaPi0; };
    TH1D* GetBetaPi0() { return BetaPi0; };

    // photon id
    TH1D* GetMomentumPhoton1() { return MomentumPhoton1; };
    TH1D* GetMomentumPhoton2() { return MomentumPhoton2; };
    TH1D* GetMomentumPhoton1_cut() { return MomentumPhoton1_cut; };
    TH1D* GetMomentumPhoton2_cut() { return MomentumPhoton2_cut; };
    TH1D* GetBetaPhoton1() { return BetaPhoton1; };
    TH1D* GetBetaPhoton2() { return BetaPhoton2; };
    TH1D* GetBetaPhoton1_cut() { return BetaPhoton1_cut; };
    TH1D* GetBetaPhoton2_cut() { return BetaPhoton2_cut; };
    TH1D* GetECuPhoton1() { return ECuPhoton1; };
    TH1D* GetECuPhoton2() { return ECuPhoton2; };
    TH1D* GetECvPhoton1() { return ECvPhoton1; };
    TH1D* GetECvPhoton2() { return ECvPhoton2; };
    TH1D* GetECwPhoton1() { return ECwPhoton1; };
    TH1D* GetECwPhoton2() { return ECwPhoton2; };
    TH1D* GetECuPhoton1_cut() { return ECuPhoton1_cut; };
    TH1D* GetECuPhoton2_cut() { return ECuPhoton2_cut; };
    TH1D* GetECvPhoton1_cut() { return ECvPhoton1_cut; };
    TH1D* GetECvPhoton2_cut() { return ECvPhoton2_cut; };
    TH1D* GetECwPhoton1_cut() { return ECwPhoton1_cut; };
    TH1D* GetECwPhoton2_cut() { return ECwPhoton2_cut; };
    TH1D* GetECtime_ECl_Start_Photon1() { return ECtime_ECl_Start_Photon1; };
    TH1D* GetECtime_ECl_Start_Photon2() { return ECtime_ECl_Start_Photon2; };
    TH1D* GetECtime_ECl_Photon1() { return ECtime_ECl_Photon1; };
    TH1D* GetECtime_ECl_Photon2() { return ECtime_ECl_Photon2; };
    TH1D* GetECtime_ECl_Photon1_cut() { return ECtime_ECl_Photon1_cut; };
    TH1D* GetECtime_ECl_Photon2_cut() { return ECtime_ECl_Photon2_cut; };
    TH1D* GetECtimePhoton1() { return ECtimePhoton1; };
    TH1D* GetECtimePhoton2() { return ECtimePhoton2; };
    TH1D* GetECpathPhoton1() { return ECpathPhoton1; };
    TH1D* GetECpathPhoton2() { return ECpathPhoton2; };
    TH1D* GetECpathtimePhoton1() { return ECpathtimePhoton1; };
    TH1D* GetECpathtimePhoton2() { return ECpathtimePhoton2; };
    TH2D* GetECtotP_vs_P_Photon1() { return ECtotP_vs_P_Photon1; };
    TH2D* GetECtotP_vs_P_Photon2() { return ECtotP_vs_P_Photon2; };
    TH2D* GetECin_vs_ECout_Photon1() { return ECin_vs_ECout_Photon1; };
    TH2D* GetECin_vs_ECout_Photon2() { return ECin_vs_ECout_Photon2; };
    TH2D* GetECtotP_vs_P_InOutZeroCut_Photon1() { return ECtotP_vs_P_InOutZeroCut_Photon1; };
    TH2D* GetECtotP_vs_P_InOutZeroCut_Photon2() { return ECtotP_vs_P_InOutZeroCut_Photon2; };
    TH2D* GetECin_vs_ECout_InOutZeroCut_Photon1() { return ECin_vs_ECout_InOutZeroCut_Photon1; };
    TH2D* GetECin_vs_ECout_InOutZeroCut_Photon2() { return ECin_vs_ECout_InOutZeroCut_Photon2; };
    TH2D* GetEC_XvsY_local_Sector_Photon1(int index) { return EC_XvsY_local_Sector_Photon1[index]; };
    TH2D* GetEC_XvsY_local_FidCut_Photon1(int index) { return EC_XvsY_local_FidCut_Photon1[index]; };
    TH2D* GetEC_XvsY_local_AntiFidCut_Photon1(int index) { return EC_XvsY_local_AntiFidCut_Photon1[index]; };
    TH2D* GetEC_XvsY_local_Sector_Photon2(int index) { return EC_XvsY_local_Sector_Photon2[index]; };
    TH2D* GetEC_XvsY_local_FidCut_Photon2(int index) { return EC_XvsY_local_FidCut_Photon2[index]; };
    TH2D* GetEC_XvsY_local_AntiFidCut_Photon2(int index) { return EC_XvsY_local_AntiFidCut_Photon2[index]; };

    //recon
    TH2D* GetMass2Pions_VS_massOmega_NC(int index) { return mass2Pions_VS_massOmega_NC[index]; };
    TH2D* GetMass2Pions_VS_massOmega_EPC(int index) { return mass2Pions_VS_massOmega_EPC[index]; };
    TH2D* GetMass2Pions_VS_massOmega_EPOC(int index) { return mass2Pions_VS_massOmega_EPOC[index]; };

    int GetMAX_SECTORS() { return MAX_SECTORS; };
    int GetID_ELECTRON() { return ID_ELECTRON; };
    int GetID_PHOTON() { return ID_PHOTON; };
    int GetID_PION_POS() { return ID_PION_POS; };
    int GetID_PION_NEG() { return ID_PION_NEG; };
    int GetID_PROTON() { return ID_PROTON; };
    double GetMASS_PHOTON() { return MASS_PHOTON; };
    double GetMASS_ELECTRON() { return MASS_ELECTRON; };
    double GetMASS_PION_CHARGED() { return MASS_PION_CHARGED; };
    double GetMASS_PION_NEUTRAL() { return MASS_PION_NEUTRAL; };
    double GetMASS_PROTON() { return MASS_PROTON; };
    double GetMASS_DEUTERIUM() { return MASS_DEUTERIUM; };
    double GetLIGHTSPEED() { return LIGHTSPEED; };
    double GetBEAM_ENERGY() { return BEAM_ENERGY; };
};
#endif
