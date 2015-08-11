#include <iostream>
#include <string>
#include "TEventReader.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "DetectedParticles.h"
#include "IntermediateParticles.h"
#include "ReconstructedParticles.h"
#include "ParticleList.h"
#include "EG2Target.h"
#include "EG2Cuts.h"
#include "ElectronID.h"
#include "PhotonID.h"
#include "ChargedPionID.h"
#include "EC_geometry.h"
#include "Vertex_Corrections.h"
#include "OmegaMixedEvent.h"
#include "HistManager.h"

#define DEBUG 0

using namespace std;

HistManager myHistManager;

int GetSectorByPhi(double phi_rad);
double Get_scMassSquared(double fMom, double fBeta);
double Get_BetaFromMass(double fMom, double fMass);
int process (string inFile, int MaxEvents, int dEvents, int targMass);
void PrintUsage(char *processName);
void PrintAnalysisTime(float tStart, float tStop);
int CheckCut(double var, double LowerLimit, double UpperLimit);
void PrintTLorentzVector(TLorentzVector TLV);

