#include <vector>
#include <string>
#include "EG2Cuts.h"
#include <iostream>
EG2Cuts::EG2Cuts()
{
    CutsLabel.push_back("NoCuts");
    CutsLabel.push_back("MassPi0");
    CutsLabel.push_back("QSquared");
    CutsLabel.push_back("Wcut");
    CutsLabel.push_back("ZDiff_Electron_PiMinus");
    CutsLabel.push_back("ZDiff_Electron_PiPlus");
    CutsLabel.push_back("OpAng_ElecPhoton");
    CutsLabel.push_back("MassPipPim");
    CutsLabel.push_back("NumDetPart");
    CutsLabel.push_back("MassOmega");
    CutsLabel.push_back("MassOmega_sideband");
    
    topo_nelec = 1;
    topo_npim = 1;
    topo_npip = 1;
    topo_ngam = 2;
    
    RangeZDiff_ElecPim.push_back(-2.0); // Lower limit on z vertex difference with electron (in cm)
    RangeZDiff_ElecPim.push_back(2.0); // Upper limit on z vertex difference with electron (in cm)
    
    RangeZDiff_ElecPip.push_back(-2.0); // Lower limit on z vertex difference with electron (in cm)
    RangeZDiff_ElecPip.push_back(2.0); // Upper limit on z vertex difference with electron (in cm)

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

    RangeMassPipPim.push_back(0.48); // Lower limit on pi+ pi- inv. mass
    RangeMassPipPim.push_back(0.51); // Upper limit on pi+ pi- inv. mass
    
    RangeOpAng_ElecPhoton.push_back(12.0); // Lower limit on opening angle between e- and photon (in degrees)
    RangeOpAng_ElecPhoton.push_back(180.0); // Upper limit on opening angle between e- and photon (in degrees)

    RangeMassOmega.push_back(0.7); // Lower limit on omega mass (in Gev/c^2)
    RangeMassOmega.push_back(0.875); // Upper limit on omega mass (in Gev/c^2)

    RangeMassOmega_sb.push_back(0.610); // Lower limit on omega mass lower sideband (in Gev/c^2)
    RangeMassOmega_sb.push_back(0.965); // Upper limit on omega mass upper sideband (in Gev/c^2)
    
    this->InitCuts();
}

// initialize the cuts
void EG2Cuts::InitCuts()
{
    cuts_omega_MPi0 = false;
    cuts_omega_MPipPim = false;
    cuts_omega_ZDiff_ElecPim = false;
    cuts_omega_ZDiff_ElecPip = false;
    cuts_omega_ZDiff = false;
    cuts_omega_Q2 = false;
    cuts_omega_W = false;
    cuts_omega_OpAng_ElecPhot1 = false;
    cuts_omega_OpAng_ElecPhot2 = false;
    cuts_omega_OpAng_ElecPhot = false;
    cuts_omega_NumDetPart = false;
    cuts_omega_MPipPim = false;
    cuts_omega_MPipPimPi0 = false;
    cuts_omega_MPipPimPi0_sb = false;
    cuts_omega_All = false;

    cuts_omega_woMPi0 = false;
    cuts_omega_woMPipPim = false;
    cuts_omega_woZDiff = false;
    cuts_omega_woQ2 = false;
    cuts_omega_woW = false;
    cuts_omega_woMPipPim = false;
    cuts_omega_woNumDetPart = false;
    cuts_omega_woOpAng_ElecPhot = false;
}

// check the cut on the difference in z vertex between e- and pi-
bool EG2Cuts::Check_ZDiff_ElecPim(double zdiff)
{
	bool ret = (zdiff >= this->Get_ZDiff_ElecPim_lo() && zdiff < this->Get_ZDiff_ElecPim_hi()) ? true : false;
	
	return ret;
}

// check the cut on the difference in z vertex between e- and pi+
bool EG2Cuts::Check_ZDiff_ElecPip(double zdiff)
{
	bool ret = (zdiff >= this->Get_ZDiff_ElecPip_lo() && zdiff < this->Get_ZDiff_ElecPip_hi()) ? true : false;
	
	return ret;
}

// set the value of the z-vertex difference cut for individual pions
// num = 0 (pi-)
// num = 1 (pi+)
void EG2Cuts::SetCut_ZDiff_ElecPion(double zdiff, int num)
{
    switch (num) {
        case 0: cuts_omega_ZDiff_ElecPim = this->Check_ZDiff_ElecPim(zdiff); break;
        case 1: cuts_omega_ZDiff_ElecPip = this->Check_ZDiff_ElecPip(zdiff); break;
        default:
            cout<<"EG2Cuts::SetCut_ZDiff_ElecPion, Wrong pion number "<<num<<endl;
            break;
    }
}

// return the value of the z-vertex difference cut for individual pions
// num = 0 (pi-)
// num = 1 (pi+)
bool EG2Cuts::GetCut_ZDiff_ElecPion(int num)
{
    bool ret;
    switch (num) {
        case 0: ret = cuts_omega_ZDiff_ElecPim; break;
        case 1: ret = cuts_omega_ZDiff_ElecPip; break;
        default:
            cout<<"EG2Cuts::GetCut_ZDiff_ElecPion, Wrong pion number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the z-vertex difference cut
void EG2Cuts::SetCut_ZDiff_ElecPion_All()
{
    cuts_omega_ZDiff = (this->GetCut_ZDiff_ElecPion(0) && this->GetCut_ZDiff_ElecPion(1));
}

// check the cut on pi0 mass
bool EG2Cuts::Check_MassPi0(double mass)
{
	bool ret = (mass >= this->Get_MassPi0_lo() && mass < this->Get_MassPi0_hi()) ? true : false;
	
	return ret;
}

// set the value of the pi0 mass cut
void EG2Cuts::SetCut_MassPi0(double mass)
{
    cuts_omega_MPi0 = this->Check_MassPi0(mass);
}

// check the cut on pi0 mass
bool EG2Cuts::Check_MassPipPim(double mass)
{
    bool ret = (mass < this->Get_MassPipPim_lo() || mass > this->Get_MassPipPim_hi()) ? true : false;
    
    return ret;
}

// set the value of the pi+pi- mass cut
void EG2Cuts::SetCut_MassPipPim(double mass)
{
    cuts_omega_MPipPim = this->Check_MassPipPim(mass);
}

// check the cut on Q^2
bool EG2Cuts::Check_QSquared(double Qsq)
{
	bool ret = (Qsq >= this->Get_QSquared_lo() && Qsq < this->Get_QSquared_hi()) ? true : false;
	
	return ret;
}

// set the value of the Q2 cut
void EG2Cuts::SetCut_QSquared(double Qsq)
{
    cuts_omega_Q2 = this->Check_QSquared(Qsq);
}

// check the cut on opening angle between e- and photon
bool EG2Cuts::Check_OpAng_ElecPhoton(double OpAng)
{
	bool ret = (OpAng >= this->Get_OpAng_ElecPhoton_lo() && OpAng < this->Get_OpAng_ElecPhoton_hi()) ? true : false;
	
	return ret;
}

// set the value of the opening angle cut between the electron and individual photons
// num = 1 (photon 1)
// num = 2 (photon 2)
void EG2Cuts::SetCut_OpAng_ElecPhoton(double OpAng, int num)
{
    switch (num) {
        case 1: cuts_omega_OpAng_ElecPhot1 = this->Check_OpAng_ElecPhoton(OpAng); break;
        case 2: cuts_omega_OpAng_ElecPhot2 = this->Check_OpAng_ElecPhoton(OpAng); break;
        default:
            cout<<"EG2Cuts::SetCut_OpAng_ElecPhoton, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the opening angle cut between the electron and individual photons
// num = 1 (photon 1)
// num = 2 (photon 2)
bool EG2Cuts::GetCut_OpAng_ElecPhoton(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_omega_OpAng_ElecPhot1; break;
        case 2: ret = cuts_omega_OpAng_ElecPhot2; break;
        default:
            cout<<"EG2Cuts::GetCut_OpAng_ElecPhoton, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the z-vertex difference cut
void EG2Cuts::SetCut_OpAng_ElecPhoton_All()
{
    cuts_omega_OpAng_ElecPhot = (this->GetCut_OpAng_ElecPhoton(1) && this->GetCut_OpAng_ElecPhoton(2));
}

// check the cut on W
bool EG2Cuts::Check_Wcut(double W)
{
    bool ret = (W >= this->Get_Wcut_lo() && W < this->Get_Wcut_hi()) ? true : false;
    
    return ret;
}

// set the value of the pi0 mass cut
void EG2Cuts::SetCut_Wcut(double W)
{
    cuts_omega_W = this->Check_Wcut(W);
}

// check the cut on particle topology
bool EG2Cuts::Check_NumDetPart(int nElec, int nPim, int nPip, int nGam)
{
    bool check_elec = (nElec == this->Get_Topo_nElec());
    bool check_pim = (nPim == this->Get_Topo_nPim());
    bool check_pip = (nPip == this->Get_Topo_nPip());
    bool check_gam = (nGam == this->Get_Topo_nGam());
    
    bool ret = (check_elec && check_pim && check_pip && check_gam) ? true : false;
    
    return ret;
}

// set the value of the particle topology cut
void EG2Cuts::SetCut_NumDetPart(int nElec, int nPim, int nPip, int nGam)
{
    cuts_omega_NumDetPart = this->Check_NumDetPart(nElec,nPim,nPip,nGam);
}

// check the cut on omega mass
bool EG2Cuts::Check_MassOmega(double mass)
{
	bool ret = (mass >= this->Get_MassOmega_lo() && mass < this->Get_MassOmega_hi()) ? true : false;
	
	return ret;
}

// set the value of the omega mass cut
void EG2Cuts::SetCut_MassOmega(double mass)
{
    cuts_omega_MPipPimPi0 = this->Check_MassOmega(mass);
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

// set the value of the sideband omega mass cut
void EG2Cuts::SetCut_MassOmega_sb(double mass)
{
    cuts_omega_MPipPimPi0_sb = this->Check_MassOmega_sb(mass);
}

// set the value of the omega candidate cut
void EG2Cuts::SetCut_OmegaID()
{
    cuts_omega_All = (this->GetCut_MassPi0() && GetCut_ZDiff_ElecPion_All() && this->GetCut_QSquared() && this->GetCut_OpAng_ElecPhoton_All() && this->GetCut_Wcut() && this->GetCut_MassPipPim());
}

// set the value of the omega candidate cut
void EG2Cuts::SetCut_OmegaID_woMassPi0()
{
    cuts_omega_woMPi0 = (this->GetCut_ZDiff_ElecPion_All() && this->GetCut_QSquared() && this->GetCut_OpAng_ElecPhoton_All() && this->GetCut_Wcut() && this->GetCut_MassPipPim());
}

// set the value of the omega candidate cut
void EG2Cuts::SetCut_OmegaID_woWcut()
{
    cuts_omega_woW = (this->GetCut_MassPi0() && GetCut_ZDiff_ElecPion_All() && this->GetCut_QSquared() && this->GetCut_OpAng_ElecPhoton_All() && this->GetCut_MassPipPim());
}

// set the value of the omega candidate cut
void EG2Cuts::SetCut_OmegaID_woQSquared()
{
    cuts_omega_woQ2 = (this->GetCut_MassPi0() && GetCut_ZDiff_ElecPion_All() && this->GetCut_OpAng_ElecPhoton_All() && this->GetCut_Wcut() && this->GetCut_MassPipPim());
}

// set the value of the omega candidate cut
void EG2Cuts::SetCut_OmegaID_woZDiff()
{
    cuts_omega_woZDiff = (this->GetCut_MassPi0() && this->GetCut_QSquared() && this->GetCut_OpAng_ElecPhoton_All() && this->GetCut_Wcut() && this->GetCut_MassPipPim());
}

// set the value of the omega candidate cut
void EG2Cuts::SetCut_OmegaID_woOpAng_ElecPhoton()
{
    cuts_omega_woOpAng_ElecPhot = (this->GetCut_MassPi0() && GetCut_ZDiff_ElecPion_All() && this->GetCut_QSquared() && this->GetCut_Wcut() && this->GetCut_MassPipPim());
}

// set the value of the omega candidate cut
void EG2Cuts::SetCut_OmegaID_woMassPipPim()
{
    cuts_omega_woMPipPim = (this->GetCut_MassPi0() && GetCut_ZDiff_ElecPion_All() && this->GetCut_QSquared() && this->GetCut_OpAng_ElecPhoton_All() && this->GetCut_Wcut());
}

// set the value of the omega candidate cut
void EG2Cuts::SetCut_OmegaID_woNumDetPart()
{
    cuts_omega_woNumDetPart = (this->GetCut_MassPi0() && GetCut_ZDiff_ElecPion_All() && this->GetCut_QSquared() && this->GetCut_OpAng_ElecPhoton_All() && this->GetCut_Wcut() && this->GetCut_MassPipPim());
}

// print the cut information
void EG2Cuts::Print_Cuts()
{
	int ii;
    cout<<"EG2 Cut Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nCuts();ii++){
        cout << this->Get_CutsLabel(ii) << "\t";
        if (this->Get_CutsLabel(ii).compare("ZDiff_Electron_PiMinus")==0) {
            cout << "[" << this->Get_ZDiff_ElecPim_lo() << "," << this->Get_ZDiff_ElecPim_hi() << "] (cm)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("ZDiff_Electron_PiPlus")==0) {
            cout << "[" << this->Get_ZDiff_ElecPip_lo() << "," << this->Get_ZDiff_ElecPip_hi() << "] (cm)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("MassPi0")==0) {
            cout << "[" << this->Get_MassPi0_lo() << "," << this->Get_MassPi0_hi() << "] (GeV/c^2)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("MassPipPim")==0) {
            cout << "[" << this->Get_MassPipPim_lo() << "," << this->Get_MassPipPim_hi() << "] (GeV/c^2)"<<endl;
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
