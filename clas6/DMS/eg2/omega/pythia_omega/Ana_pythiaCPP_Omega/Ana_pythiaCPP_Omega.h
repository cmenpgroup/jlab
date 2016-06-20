#include <iostream>
#include <fstream>
#include <string>
#include "HistManager_pythiaCPP_Omega.h"
#include "PartComb_pythiaCPP_Omega.h"
#include "PDG_pythiaCPP_Omega.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2D.h"

using namespace std;

#define MAX_ANA 2
#define MAX_TRACKS 50
#define LIMIT_QSQ 1.0
#define LIMIT_W 2.0

HistManager_pythiaCPP_Omega myHistManager; // declare the histogram manager

int ID_ELECTRON = 11; // PDG electron id
int ID_PHOTON = 22;  // PDG photon id
int ID_PION_POS = 211;  // PDG pi+ id
int ID_PION_NEG = -211;  // PDG pi- id
int ID_PION_ZERO = 111;  // PDG pi0 id
int ID_ETA_MESON = 221;  // PDG eta id
int ID_OMEGA_MESON = 223;  // PDG omega id
int ID_ETA_PRIME_MESON = 331;  // PDG eta' (958) id
int ID_PHI_MESON = 333;  // PDG phi id
int ID_PROTON = 2212; // PDG proton id

double MASS_PHOTON = 0.0; // mass of photon in GeV/c^2
double MASS_ELECTRON = 0.000511; // mass of charged pion in GeV/c^2
double MASS_PION_CHARGED = 0.138; // mass of charged pion in GeV/c^2
double MASS_PION_NEUTRAL = 0.135; // mass of neutral pion in GeV/c^2
double MASS_PROTON = 0.938; // mass of proton in GeV/c^2

int process (string inFile, int iAna, int MaxEvents, int dEvents);
bool Cut_CLAS6_Theta_EC(double theta);
bool Cut_CLAS6_Theta_TOF(double theta);
bool Cut_CLAS6_Theta_Ana(TLorentzVector V4, int iAna, string detName);
void FillHists_omega(TLorentzVector pip, TLorentzVector pim, TLorentzVector pi0, double nu, int iPC);