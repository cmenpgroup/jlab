#include <iostream>
#include <string>
#include "TEventReader.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

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

//double BEAM_ENERGY = 4.5; // electron beam energy in GeV
double BEAM_ENERGY = 5.01; // electron beam energy in GeV

TH1D *q2;
TH2D *q2_VS_theta;
TH1D *nu_EnergyTransfer;
TH1D *elecZVert;
TH1D *OpAng_2Photons;
TH1D *OpAng_elecPhoton1;
TH1D *OpAng_elecPhoton2;
TH2D *ZVertDiff;
TH2D *Xvert_VS_Yvert[5];
TH2D *Beta_VS_Momentum;
TH2D *Theta_VS_Phi[7];
TH2D *TotalMomentum;

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
TH1D *elecZVertSector[6];
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

int GetSectorByPhi(double phi_rad);
int process (string inFile, int MaxEvents, int dEvents, int targMass);
void PrintUsage(char *processName);
void PrintAnalysisTime(float tStart, float tStop);
int CheckCut(double var, double LowerLimit, double UpperLimit);
void BookHist();
void WriteHist(string RootFile);
void PrintTLorentzVector(TLorentzVector TLV);

