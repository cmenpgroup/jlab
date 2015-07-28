#include <iostream>
#include <string>
#include "TEventReader.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "DetectedParticles.h"
#include "IntermediateParticles.h"
#include "ReconstructedParticles.h"
#include "ParticleList.h"
#include "EG2Target.h"
#include "EG2Cuts.h"
#include "ElectronID.h"
#include "PhotonID.h"
#include "EC_geometry.h"
#include "Vertex_Corrections.h"
#include "OmegaMixedEvent.h"

#define DEBUG 0

using namespace std;

int MAX_SECTORS = 6; // max. number of CLAS sectors

int ID_ELECTRON = 11; // PDG electron id
int ID_PHOTON = 22;  // PDG photon id
int ID_PION_POS = 211;  // PDG pi+ id
int ID_PION_NEG = -211;  // PDG pi- id
int ID_PROTON = 2212; // PDG proton id

double MASS_PHOTON = 0.0; // mass of photon in GeV/c^2
double MASS_ELECTRON = 0.000511; // mass of charged pion in GeV/c^2
double MASS_PION_CHARGED = 0.138; // mass of charged pion in GeV/c^2
double MASS_PION_NEUTRAL = 0.135; // mass of neutral pion in GeV/c^2
double MASS_PROTON = 0.938; // mass of proton in GeV/c^2
double MASS_DEUTERIUM = 2*MASS_PROTON; // mass of deuterium in GeV/c^2

double LIGHTSPEED = 30.0; // speed of light in cm/ns

//double BEAM_ENERGY = 4.5; // electron beam energy in GeV
double BEAM_ENERGY = 5.01; // electron beam energy in GeV

TH1D *q2;
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
TH2D *Beta_Recalc;
TH2D *Theta_VS_Phi[7];
TH2D *TotalMomentum;
TH2D *elecZVert_VS_Phi;
TH2D *elecZVert_VS_Phi_Corr;

TH2D *CCnphe;
TH2D *ECu;
TH2D *ECv;
TH2D *ECw;
TH2D *dtime_ECSC;
TH2D *ECtot_VS_P[5];
TH2D *ECtotP_VS_P[5];
TH2D *ECin_VS_ECout[5];

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

TH2D *ECinP_VS_ECoutP[6]; //Andy
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
TH1D *Pi0Mass_PhotIDcuts;
TH1D *OmegaMass_AllCuts;

int GetSectorByPhi(double phi_rad);
int process (string inFile, int MaxEvents, int dEvents, int targMass);
void PrintUsage(char *processName);
void PrintAnalysisTime(float tStart, float tStop);
int CheckCut(double var, double LowerLimit, double UpperLimit);
void BookHist();
void WriteHist(string RootFile);
void PrintTLorentzVector(TLorentzVector TLV);

