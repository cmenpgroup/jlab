#ifndef HISTMANAGER_H
#define HISTMANAGER_H
#include <vector>
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "PartComb_pythiaCPP_Omega.h"

class HistManager_pythiaCPP_Omega
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
    
    TH1D *hIntermedPart;
    TH2D *hPartPerEvt;
    TH2D *hIMomega;
    TH2D *hMMsq_X;
    TH2D *hMMsq_Xrecoil;
    TH2D *hMMsq_Xpip;
    TH2D *hMMsq_Xpim;
    TH2D *hMMsq_Xpi0;
    TH2D *hOpAng_PipPim;
    TH2D *hOpAng_PipPimPi0;
    TH1D *hQsq_NoCuts;
    TH1D *hNu_NoCuts;
    TH1D *hW_NoCuts;
    TH1D *hQsq;
    TH1D *hNu;
    TH1D *hW;
    TH2D *hz_fracEnergy;
    TH2D *hOpAng_TwoPhoton;
    TH2D *hIMomega_VS_IMPipPim[16];  // one histogram for each particle combination
    TH2D *hDalitz_pip[16];  // one Dalitz histogram for each particle combination
    
public:
    HistManager_pythiaCPP_Omega();
    void BookHist();
    void WriteHist(std::string outFile);
    TH1D* Get_hIntermedPart() { return hIntermedPart; };
    TH1D* Get_hQsq_NoCuts() { return hQsq_NoCuts; };
    TH1D* Get_hNu_NoCuts() { return hNu_NoCuts; };
    TH1D* Get_hW_NoCuts() { return hW_NoCuts; };
    TH1D* Get_hQsq() { return hQsq; };
    TH1D* Get_hNu() { return hNu; };
    TH1D* Get_hW() { return hW; };
    
    TH2D* Get_hPartPerEvt() { return hPartPerEvt; };
    TH2D* Get_hIMomega() { return hIMomega; };
    TH2D* Get_hMMsq_X() { return hMMsq_X; };
    TH2D* Get_hMMsq_Xrecoil() { return hMMsq_Xrecoil; };
    TH2D* Get_hMMsq_Xpip() { return hMMsq_Xpip; };
    TH2D* Get_hMMsq_Xpim() { return hMMsq_Xpim; };
    TH2D* Get_hMMsq_Xpi0() { return hMMsq_Xpi0; };
    TH2D* Get_hOpAng_PipPim() { return hOpAng_PipPim; };
    TH2D* Get_hOpAng_PipPimPi0() { return hOpAng_PipPimPi0; };
    TH2D* Get_hz_fracEnergy() { return hz_fracEnergy; };
    TH2D* Get_hOpAng_TwoPhoton() { return hOpAng_TwoPhoton; };
    
    TH2D* Get_hIMomega_VS_IMPipPim(int index) { return hIMomega_VS_IMPipPim[index]; };
    TH2D* Get_hDalitz_pip(int index) { return hDalitz_pip[index]; };
    
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