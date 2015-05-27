#include <vector>
#include <string>
#include "EG2Cuts.h"
#include <iostream>
EG2Cuts::EG2Cuts()
{
    CutsLabel.push_back("NoCuts");
    CutsLabel.push_back("MassPi0");
    CutsLabel.push_back("Zdiff_Electron_PiMinus");
    CutsLabel.push_back("Zdiff_Electron_PiPlus");
    CutsLabel.push_back("QSquared");
    CutsLabel.push_back("OpAng_ElecPhoton");
    CutsLabel.push_back("BetaPhoton");
    CutsLabel.push_back("Wcut");
    CutsLabel.push_back("ElectronR");
    CutsLabel.push_back("MassOmega");
    CutsLabel.push_back("MassOmega_sideband");
    
    RangeZdiff_ElecPim.push_back(-2.0); // Lower limit on z vertex difference with electron (in cm)
    RangeZdiff_ElecPim.push_back(2.0); // Upper limit on z vertex difference with electron (in cm)
    
    RangeZdiff_ElecPip.push_back(-2.0); // Lower limit on z vertex difference with electron (in cm)
    RangeZdiff_ElecPip.push_back(2.0); // Upper limit on z vertex difference with electron (in cm)

    double pi0Centroid = 0.134;
    double pi0Width = 0.025;
    double pi0Nsigmas = 3.0;
    double pi0MLo = pi0Centroid - pi0Nsigmas*pi0Width;
    double pi0MHi = pi0Centroid + pi0Nsigmas*pi0Width;
    RangeMassPi0.push_back(pi0MLo); // Lower limit on pi0 mass (in Gev/c^2)
    RangeMassPi0.push_back(pi0MHi); // Upper limit on pi0 mass (in Gev/c^2)

    RangeQSquared.push_back(1.0); // Lower limit on Q^2 (in Gev^2)
    RangeQSquared.push_back(100000.0); // Upper limit on Q^2 (in Gev^2)

    RangeWcut.push_back(2.0); // Lower limit on W (in Gev)
    RangeWcut.push_back(100000.0); // Upper limit on W (in Gev)

    RangeElecR.push_back(0.0); // Lower limit on electron radial vertex (in cm)
    RangeElecR.push_back(0.1); // Upper limit on electron radial vertex (in cm)
    
    RangeBetaPhoton.push_back(0.95); // Lower limit on photon beta
    RangeBetaPhoton.push_back(1.05); // Upper limit photon beta
    
    RangeOpAng_ElecPhoton.push_back(12.0); // Lower limit on opening angle between e- and photon (in degrees)
    RangeOpAng_ElecPhoton.push_back(180.0); // Upper limit on opening angle between e- and photon (in degrees)

    RangeMassOmega.push_back(0.7); // Lower limit on omega mass (in Gev/c^2)
    RangeMassOmega.push_back(0.875); // Upper limit on omega mass (in Gev/c^2)

    RangeMassOmega_sb.push_back(0.610); // Lower limit on omega mass lower sideband (in Gev/c^2)
    RangeMassOmega_sb.push_back(0.965); // Upper limit on omega mass upper sideband (in Gev/c^2)
}

// check the cut on the difference in z vertex between e- and pi-
bool EG2Cuts::Check_Zdiff_ElecPim(double zdiff)
{
	bool ret = (zdiff >= this->Get_Zdiff_ElecPim_lo() && zdiff < this->Get_Zdiff_ElecPim_hi()) ? true : false;
	
	return ret;
}

// check the cut on the difference in z vertex between e- and pi+
bool EG2Cuts::Check_Zdiff_ElecPip(double zdiff)
{
	bool ret = (zdiff >= this->Get_Zdiff_ElecPip_lo() && zdiff < this->Get_Zdiff_ElecPip_hi()) ? true : false;
	
	return ret;
}

// check the cut on pi0 mass
bool EG2Cuts::Check_MassPi0(double mass)
{
	bool ret = (mass >= this->Get_MassPi0_lo() && mass < this->Get_MassPi0_hi()) ? true : false;
	
	return ret;
}

// check the cut on Q^2
bool EG2Cuts::Check_QSquared(double Qsq)
{
	bool ret = (Qsq >= this->Get_QSquared_lo() && Qsq < this->Get_QSquared_hi()) ? true : false;
	
	return ret;
}

// check the cut on opening angle between e- and photon
bool EG2Cuts::Check_OpAng_ElecPhoton(double OpAng)
{
	bool ret = (OpAng >= this->Get_OpAng_ElecPhoton_lo() && OpAng < this->Get_OpAng_ElecPhoton_hi()) ? true : false;
	
	return ret;
}

// check the cut on W
bool EG2Cuts::Check_Wcut(double W)
{
    bool ret = (W >= this->Get_Wcut_lo() && W < this->Get_Wcut_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron radius
bool EG2Cuts::Check_ElectronR(double vr)
{
    bool ret = (vr >= this->Get_ElectronR_lo() && vr < this->Get_ElectronR_hi()) ? true : false;
    
    return ret;
}

// check the cut on photon beta
bool EG2Cuts::Check_BetaPhoton(double beta)
{
    bool ret = (beta >= this->Get_BetaPhoton_lo() && beta < this->Get_BetaPhoton_hi()) ? true : false;
    
    return ret;
}

// check the cut on omega mass
bool EG2Cuts::Check_MassOmega(double mass)
{
	bool ret = (mass >= this->Get_MassOmega_lo() && mass < this->Get_MassOmega_hi()) ? true : false;
	
	return ret;
}

// check the cut on omega mass sidebands
bool EG2Cuts::Check_MassOmega_sb(double mass)
{
    int lower_sb = 0;
    int upper_sb = 0;

    lower_sb = (mass >= this->Get_MassOmega_sb_lo() && mass < this->Get_MassOmega_lo());
    upper_sb = (mass >= this->Get_MassOmega_hi() && mass < this->Get_MassOmega_sb_hi());
    
    bool ret = (lower_sb || upper_sb) ? true : false;
    
    return ret;
}

// print the cut information
void EG2Cuts::Print_Cuts()
{
	int ii;
    cout<<"EG2 Cut Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nCuts();ii++){
        cout << this->Get_CutsLabel(ii) << "\t";
        if (this->Get_CutsLabel(ii).compare("Zdiff_Electron_PiMinus")==0) {
            cout << "[" << this->Get_Zdiff_ElecPim_lo() << "," << this->Get_Zdiff_ElecPim_hi() << "] (cm)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("Zdiff_Electron_PiPlus")==0) {
            cout << "[" << this->Get_Zdiff_ElecPip_lo() << "," << this->Get_Zdiff_ElecPip_hi() << "] (cm)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("MassPi0")==0) {
            cout << "[" << this->Get_MassPi0_lo() << "," << this->Get_MassPi0_hi() << "] (GeV/c^2)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("QSquared")==0) {
            cout << "[" << this->Get_QSquared_lo() << "," << this->Get_QSquared_hi() << "] (GeV^2)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("Wcut")==0) {
            cout << "[" << this->Get_Wcut_lo() << "," << this->Get_Wcut_hi() << "] (GeV)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("BetaPhoton")==0) {
            cout << "[" << this->Get_BetaPhoton_lo() << "," << this->Get_BetaPhoton_hi() << "]" << endl;
        }else if (this->Get_CutsLabel(ii).compare("OpAng_ElecPhoton")==0) {
            cout << "[" << this->Get_OpAng_ElecPhoton_lo() << "," << this->Get_OpAng_ElecPhoton_hi() << "] (deg.)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("ElectronR")==0) {
            cout << "[" << this->Get_ElectronR_lo() << "," << this->Get_ElectronR_hi() << "] (cm)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("MassOmega")==0) {
            cout << "[" << this->Get_MassOmega_lo() << "," << this->Get_MassOmega_hi() << "] (GeV/c^2)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("MassOmega_sideband")==0) {
            cout << "[" << this->Get_MassOmega_sb_lo() << "," << this->Get_MassOmega_lo() << " or "<< this->Get_MassOmega_hi() << "," << this->Get_MassOmega_sb_hi() << "] (GeV/c^2)" << endl;
        }else{
            cout << endl;
        }
    }
    cout << endl;
}
