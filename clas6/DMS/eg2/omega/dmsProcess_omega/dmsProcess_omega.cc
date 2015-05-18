/************************************************************************/
/*  dmsProcess_omega.cc                                                 */
/*                                                                      */
/*  Created by Angelo Licastro and Andy Beiter, Canisius College        */
/*  July 2014 - Modified by M. H. Wood, Canisius College                */
/*                                                                      */
/************************************************************************/

#include "dmsProcess_omega.h"

class DetectedParticles
{
public:
	vector<string> DetPartLabel;
	DetectedParticles();
	int Get_nDetPartLabel() {return DetPartLabel.size();};
	string Get_DetPartLabel(int num) {return DetPartLabel[num];};
	void Print_DetPartLabel();
};

DetectedParticles::DetectedParticles()
{
	DetPartLabel.push_back("Electron");
	DetPartLabel.push_back("Pi-");
	DetPartLabel.push_back("Pi+");
	DetPartLabel.push_back("Photon1");
	DetPartLabel.push_back("Photon2");
}

void DetectedParticles::Print_DetPartLabel()
{
	int ii;
    cout<<"Detected Particles"<<endl;
    cout<<"=================="<<endl;
	for(ii=0;ii<this->Get_nDetPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_DetPartLabel(ii) << endl;
	}
}

class IntermediateParticles
{
public:
	vector<string> IntPartLabel;
	IntermediateParticles();
	int Get_nIntPartLabel() {return IntPartLabel.size();};
	string Get_IntPartLabel(int num) {return IntPartLabel[num];};
	void Print_IntPartLabel();
};

IntermediateParticles::IntermediateParticles()
{
	IntPartLabel.push_back("Pi0");
}

void IntermediateParticles::Print_IntPartLabel()
{
	int ii;
    cout<<"Intermediate Particles"<<endl;
    cout<<"======================"<<endl;
	for(ii=0;ii<this->Get_nIntPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_IntPartLabel(ii) << endl;
	}
}

class ReconstructedParticles
{
public:
	vector<string> RecPartLabel;
	ReconstructedParticles();
	int Get_nRecPartLabel() {return RecPartLabel.size();};
	string Get_RecPartLabel(int num) {return RecPartLabel[num];};
	void Print_RecPartLabel();
};

ReconstructedParticles::ReconstructedParticles()
{
	RecPartLabel.push_back("Omega");
}

void ReconstructedParticles::Print_RecPartLabel()
{
	int ii;
    cout<<"Reconstructed Particles"<<endl;
    cout<<"======================="<<endl;
	for(ii=0;ii<this->Get_nRecPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_RecPartLabel(ii) << endl;
	}
}

class ParticleList
{
    vector<string> PartLabel;
public:
    ParticleList();
    int Get_nPartLabel() {return PartLabel.size();};
	string Get_PartLabel(int num) {return PartLabel[num];};
	void Print_PartLabel();
};

ParticleList::ParticleList()
{
    vector<string> temp;
    DetectedParticles DetList;
    IntermediateParticles IntList;
    ReconstructedParticles RecList;
 
    int ii;
	for(ii=0;ii<DetList.Get_nDetPartLabel();ii++){
        PartLabel.push_back(DetList.Get_DetPartLabel(ii));
    }
	for(ii=0;ii<IntList.Get_nIntPartLabel();ii++){
        PartLabel.push_back(IntList.Get_IntPartLabel(ii));
    }
	for(ii=0;ii<RecList.Get_nRecPartLabel();ii++){
        PartLabel.push_back(RecList.Get_RecPartLabel(ii));
    }
}
void ParticleList::Print_PartLabel()
{
	int ii;
    cout<<"All Particles in Analysis"<<endl;
    cout<<"========================="<<endl;
	for(ii=0;ii<this->Get_nPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_PartLabel(ii) << endl;
	}
}

class EG2Target
{
    vector<string> Label;
    vector<int> Index;
    vector<double> RangeLD2;
    vector<double> RangeNuc;
public:
    EG2Target();
    int Get_nLabel() {return Label.size();};
    int Get_nIndex() {return Index.size();};
	string Get_Label(int num) {return Label[num];};
    double Get_LD2_lo() {return RangeLD2[0];};
    double Get_LD2_hi() {return RangeLD2[1];};
    double Get_Nuc_lo() {return RangeNuc[0];};
    double Get_Nuc_hi() {return RangeNuc[1];};
    int Get_Index(double z);
	void Print_Info();
};

EG2Target::EG2Target()
{
    Label.push_back("NoTarget");
    Label.push_back("LD2");
    Label.push_back("Nuc");
    
    Index.push_back(0);
    Index.push_back(1);
    Index.push_back(2);
    
    RangeLD2.push_back(-32.0);
    RangeLD2.push_back(-28.0);

    RangeNuc.push_back(-26.0);
    RangeNuc.push_back(-23.0);
}

// Return the CLAS eg2 target index from vertex Z.
//
// Return 0 = outside target limits, 1 = liquid deuterium, 2 = nuclear
// z must be given in cm
//
int EG2Target::Get_Index(double z){
    
    Int_t ret = 0; // init the return variable
    
    if (z >= this->Get_LD2_lo() && z < this->Get_LD2_hi()) {
        ret = Index[1]; // deuterium target
    } else if (z >= this->Get_Nuc_lo() && z < this->Get_Nuc_hi()) {
        ret = Index[2]; // nuclear target
    } else {
        ret = Index[0]; // no target
    }
    
    return ret;
}

void EG2Target::Print_Info()
{
	int ii;
    cout<<"EG2 Target Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nLabel();ii++){
        cout << ii+1 << "\t" << this->Get_Label(ii) << endl;
    }
    
    cout << "LD2 target: " << this->Get_LD2_lo() << " , " << this->Get_LD2_hi() << " (cm)" << endl;
    cout << "Nuclear target: " << this->Get_Nuc_lo() << " , " << this->Get_Nuc_hi() << " (cm)" << endl;
}

class EG2Cuts
{
    vector<string> CutsLabel;
    vector<double> RangeZdiff_ElecPim;
    vector<double> RangeZdiff_ElecPip;
    vector<double> RangeMassPi0;
    vector<double> RangeQSquared;
    vector<double> RangeOpAng_ElecPhoton;
    vector<double> RangeBetaPhoton;
    vector<double> RangeWcut;
    vector<double> RangeElecR;
    vector<double> RangeMassOmega;
    vector<double> RangeMassOmega_sb;
public:
    EG2Cuts();
    int Get_nCuts() {return CutsLabel.size();};
	string Get_CutsLabel(int num) {return CutsLabel[num];};
    double Get_Zdiff_ElecPim_lo() {return RangeZdiff_ElecPim[0];};
    double Get_Zdiff_ElecPim_hi() {return RangeZdiff_ElecPim[1];};
    double Get_Zdiff_ElecPip_lo() {return RangeZdiff_ElecPip[0];};
    double Get_Zdiff_ElecPip_hi() {return RangeZdiff_ElecPip[1];};
    double Get_MassPi0_lo() {return RangeMassPi0[0];};
    double Get_MassPi0_hi() {return RangeMassPi0[1];};
    double Get_QSquared_lo() {return RangeQSquared[0];};
    double Get_QSquared_hi() {return RangeQSquared[1];};
    double Get_BetaPhoton_lo() {return RangeBetaPhoton[0];};
    double Get_BetaPhoton_hi() {return RangeBetaPhoton[1];};
    double Get_Wcut_lo() {return RangeWcut[0];};
    double Get_Wcut_hi() {return RangeWcut[1];};
    double Get_OpAng_ElecPhoton_lo() {return RangeOpAng_ElecPhoton[0];};
    double Get_OpAng_ElecPhoton_hi() {return RangeOpAng_ElecPhoton[1];};
    double Get_ElectronR_lo() {return RangeElecR[0];};
    double Get_ElectronR_hi() {return RangeElecR[1];};
    double Get_MassOmega_lo() {return RangeMassOmega[0];};
    double Get_MassOmega_hi() {return RangeMassOmega[1];};
    double Get_MassOmega_sb_lo() {return RangeMassOmega_sb[0];};
    double Get_MassOmega_sb_hi() {return RangeMassOmega_sb[1];};
    bool Check_Zdiff_ElecPim(double zdiff);
    bool Check_Zdiff_ElecPip(double zdiff);
    bool Check_MassPi0(double mass);
    bool Check_QSquared(double Qsq);
    bool Check_OpAng_ElecPhoton(double OpAng);
    bool Check_Wcut(double W);
    bool Check_ElectronR(double vr);
    bool Check_BetaPhoton(double beta);
    bool Check_MassOmega(double mass);
    bool Check_MassOmega_sb(double mass);
    void Print_Cuts();
};

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

class ElectronID
{
    vector<string> elecIDLabel;
    vector<double> RangeElecMom;
    vector<double> RangeECu;
    vector<double> RangeECv;
    vector<double> RangeECw;
    vector<double> RangeECin;
    vector<double> Range_dtECSC;
    vector<double> RangeCCnphe;
    
    // parameters to calculate the EC sampling fraction of total energy vs P
    double EC_SamplingFrac_C[6][5] = {{2.52E-1,1.22E-2,-7.94E-3,9.55E-3,3.41E-2},
        {2.78E-1,1.87E-2,-2.38E-3,1.399E-2,3.75E-2},
        {2.62E-1,2.31E-2,-3.54E-3,9.32E-3,2.90E-2},
        {2.51E-1,2.01E-2,-3.32E-3,8.21E-3,2.99E-2},
        {2.63E-1,9.55E-2,-1.02E-3,2.25E-2,3.06E-2},
        {2.55E-1,2.32E-2,-3.05E-3,1.17E-2,3.64E-2}};

    double EC_SamplingFrac_Fe[6][5] = {{2.22E-1,2.23E-2,-2.41E-3,9.23E-3,2.98E-2},
        {2.34E-1,1.95E-2,-2.08E-3,8.66E-3,3.09E-2},
        {2.52E-1,2.42E-2,-3.39E-3,1.08E-2,2.64E-2},
        {2.51E-1,2.08E-2,-3.27E-3,7.22E-3,2.98E-2},
        {2.72E-1,1.18E-2,-1.87E-3,1.84E-2,3.48E-2},
        {2.52E-1,2.28E-2,-3.11E-3,4.11E-3,3.55E-2}};
    
    double EC_SamplingFrac_Pb[6][5] = {{2.53E-1,1.38E-2,-1.40E-3,7.67E-3,3.54E-2},
        {2.49E-1,1.47E-2,-1.49E-3,7.53E-3,3.38E-2},
        {2.54E-1,2.26E-2,-3.05E-3,8.13E-3,2.77E-2},
        {2.55E-1,1.90E-2,-3.05E-3,7.20E-3,3.04E-2},
        {2.76E-1,1.11E-2,-1.76E-3,1.81E-2,3.53E-2},
        {2.62E-1,1.92E-2,-2.62E-3,1.99E-3,3.76E-2}};
    
public:
    ElectronID();
    int Get_nElecID() {return elecIDLabel.size();};
    string Get_elecIDLabel(int num) {return elecIDLabel[num];};
    double Get_ElecMom_lo() {return RangeElecMom[0];};
    double Get_ElecMom_hi() {return RangeElecMom[1];};
    double Get_ElecECu_lo() {return RangeECu[0];};
    double Get_ElecECu_hi() {return RangeECu[1];};
    double Get_ElecECv_lo() {return RangeECv[0];};
    double Get_ElecECv_hi() {return RangeECv[1];};
    double Get_ElecECw_lo() {return RangeECw[0];};
    double Get_ElecECw_hi() {return RangeECw[1];};
    double Get_ElecECin_lo() {return RangeECin[0];};
    double Get_ElecECin_hi() {return RangeECin[1];};
    double Get_ElecCCnphe_lo() {return RangeCCnphe[0];};
    double Get_ElecCCnphe_hi() {return RangeCCnphe[1];};
    double Get_Elec_dtECSC_lo() {return Range_dtECSC[0];};
    double Get_Elec_dtECSC_hi() {return Range_dtECSC[1];};
    double Get_EC_SamplingFraction(int coeff, int sector, int targMass);
    
    bool Check_ElecMom(double mom);
    bool Check_ElecECu(double ecu);
    bool Check_ElecECv(double ecv);
    bool Check_ElecECw(double ecw);
    bool Check_ElecECin(double ecin);
    bool Check_Elec_dtECSC(double dt);
    bool Check_ElecCCnphe(double nphe);
    bool Check_ElecECoverP(double mom, double ectot, int sector, int targMass);
    
    void Print_ElectronID();
};

ElectronID::ElectronID()
{
    elecIDLabel.push_back("No Cuts");
    elecIDLabel.push_back("Momentum");
    elecIDLabel.push_back("EC U-view");
    elecIDLabel.push_back("EC V-view");
    elecIDLabel.push_back("EC W-view");
    elecIDLabel.push_back("ECtot/P VS P");
    elecIDLabel.push_back("ECin");
    elecIDLabel.push_back("CC Nphe");
    elecIDLabel.push_back("dt(EC-SC)");
    
    RangeElecMom.push_back(0.64); // Lower limit on e- momentum (in GeV)
    RangeElecMom.push_back(1000.0); // Upper limit on e- momentum (in GeV)
    
    RangeECu.push_back(40); // Lower limit on EC U-view (in cm)
    RangeECu.push_back(400); // Upper limit on EC U-view (in cm)

    RangeECv.push_back(0); // Lower limit on EC V-view (in cm)
    RangeECv.push_back(360); // Upper limit on EC V-view (in cm)

    RangeECw.push_back(0); // Lower limit on EC W-view (in cm)
    RangeECw.push_back(390); // Upper limit on EC W-view (in cm)

    RangeECin.push_back(0.06); // Lower limit on EC inner energy (in GeV)
    RangeECin.push_back(10.0); // Upper limit on EC inner energy (in GeV)
    
    RangeCCnphe.push_back(0); // Lower limit on CC num. photo-electrons
    RangeCCnphe.push_back(1000); // Upper limit on CC num. photo-electrons

    double dtCentroid = 0.0;
    double dtWidth = 0.06;
    double dtNsigmas = 5.0;
    double dtLo = dtCentroid - dtNsigmas*dtWidth;
    double dtHi = dtCentroid + dtNsigmas*dtWidth;
    Range_dtECSC.push_back(dtLo); // Lower limit on time difference between EC and SC (in ns)
    Range_dtECSC.push_back(dtHi); // Upper limit on time difference between EC and SC (in ns)
    
}

// check the cut on electron momentum
bool ElectronID::Check_ElecMom(double mom)
{
    bool ret = (mom >= this->Get_ElecMom_lo() && mom < this->Get_ElecMom_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron EC U-view
bool ElectronID::Check_ElecECu(double ecu)
{
    bool ret = (ecu >= this->Get_ElecECu_lo() && ecu < this->Get_ElecECu_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron EC V-view
bool ElectronID::Check_ElecECv(double ecv)
{
    bool ret = (ecv >= this->Get_ElecECv_lo() && ecv < this->Get_ElecECv_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron EC W-view
bool ElectronID::Check_ElecECw(double ecw)
{
    bool ret = (ecw >= this->Get_ElecECw_lo() && ecw < this->Get_ElecECw_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron EC inner energy
bool ElectronID::Check_ElecECin(double ecin)
{
    bool ret = (ecin >= this->Get_ElecECin_lo() && ecin < this->Get_ElecECin_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron time difference between EC and SC
bool ElectronID::Check_Elec_dtECSC(double dt)
{
    bool ret = (dt >= this->Get_Elec_dtECSC_lo() && dt < this->Get_Elec_dtECSC_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron CC num. of photo-electrons
bool ElectronID::Check_ElecCCnphe(double nphe)
{
    bool ret = (nphe >= this->Get_ElecCCnphe_lo() && nphe < this->Get_ElecCCnphe_hi()) ? true : false;
    
    return ret;
}

double ElectronID::Get_EC_SamplingFraction(int coeff, int sector, int targMass)
{
    double ret = 0.0;
    
    if(sector>=1 && sector<=6){ //check that the sector is between 1 and 6
        if(coeff>=0 && coeff<5){
            switch (targMass){
                case 12: ret = this->EC_SamplingFrac_C[sector-1][coeff]; break;
                case 56: ret = this->EC_SamplingFrac_Fe[sector-1][coeff]; break;
                case 208: ret = this->EC_SamplingFrac_Pb[sector-1][coeff]; break;
                default:
                    cout<<"ElectronID::Get_EC_SamplingFraction: Target Mass "<< targMass <<" is unknown."<<endl;
                    ret = 0.0;
                    break;
            }
        }
        else{
            cout<<"ElectronID::Get_EC_SamplingFraction: Coefficient "<<coeff<<" is out of range."<<endl;
        }
    }
    else{
        cout<<"ElectronID::Get_EC_SamplingFraction: Sector "<<sector<<" is out of range."<<endl;
    }
    return ret;
}

// check the cut on electron EC inner energy
bool ElectronID::Check_ElecECoverP(double mom, double ectot, int sector, int targMass)
{
    bool ret = false; // initialize to false
    
    double a = this->Get_EC_SamplingFraction(0,sector,targMass);
    double b = this->Get_EC_SamplingFraction(1,sector,targMass);
    double c = this->Get_EC_SamplingFraction(2,sector,targMass);
    double d = this->Get_EC_SamplingFraction(3,sector,targMass);
    double f = this->Get_EC_SamplingFraction(4,sector,targMass);

    double centroid = a + b*mom + c*mom*mom;
    double sigma = sqrt(d*d + f*f/sqrt(mom));
    double Nsigma = 2.5;
    
    double diff = fabs(ectot/mom - centroid);
    
    ret = (diff < Nsigma*sigma) ? true : false;

    return ret;
}

// print the cut information
void ElectronID::Print_ElectronID()
{
    int ii;
    cout<<"Electron ID Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nElecID();ii++){
        cout << this->Get_elecIDLabel(ii) << "\t";
        if (this->Get_elecIDLabel(ii).compare("Momentum")==0) {
            cout << "[" << this->Get_ElecMom_lo() << "," << this->Get_ElecMom_hi() << "] (GeV)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("EC U-view")==0) {
            cout << "[" << this->Get_ElecECu_lo() << "," << this->Get_ElecECu_hi() << "] (cm)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("EC V-view")==0) {
            cout << "[" << this->Get_ElecECv_lo() << "," << this->Get_ElecECv_hi() << "] (cm)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("EC W-view")==0) {
            cout << "[" << this->Get_ElecECw_lo() << "," << this->Get_ElecECw_hi() << "] (cm)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("ECin")==0) {
            cout << "[" << this->Get_ElecECin_lo() << "," << this->Get_ElecECin_hi() << "] (GeV)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("CC Nphe")==0) {
            cout << "[" << this->Get_ElecCCnphe_lo() << "," << this->Get_ElecCCnphe_hi() << "]" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("dt(EC-SC)")==0) {
            cout << "[" << this->Get_Elec_dtECSC_lo() << "," << this->Get_Elec_dtECSC_hi() << "] (ns)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("ECtot/P VS P")==0) {
            cout << "Depends on the momentum dependent sampling fraction" << endl;
        }else{
            cout << endl;
        }
    }
    cout << endl;
}

class EC_geometry
{
    double EC_U;
    double EC_V;
    double EC_W;
    double EC_theta;
    double EC_phi;
    double ylow;
    double yhi;
    double rho;
public:
    EC_geometry();
    void Put_UVW(double u, double v, double w);
    double Get_EC_theta() {return EC_theta};
    double Get_ylow() {return ylow};
    double Get_yhi() {return yhi};
    double Get_rho() {return rho};
    double Get_U() {return EC_U};
    double Get_V() {return EC_V};
    double Get_W() {return EC_W};
    double Get_Xlocal();
    double Get_Ylocal();
    bool Check_U();
    bool Check_V();
    bool Check_W();
    void Print_EC_geometry();
};

EC_geometry::EC_geometry()
{
    EC_theta = 0.4363323; // radians
    ylow = -182.974;
    yhi = 189.956;
    rho = 1.097621; // radians
    
    EC_U = -1.0;
    EC_V = -1.0;
    EC_W = -1.0;
}

double EC_geometry::Get_Xlocal()
{
    double ret = 0.0;
    if(this->Check_U() && this->Check_V() && this->Check_W()) ret = this->Get_W()*cos(this->Get_rho()) - 0.5*this->Get_V();
    return ret;
}

double EC_geometry::Get_Ylocal()
{
    double ret = 0.0;
    if(this->Check_U() && this->Check_V() && this->Check_W()) ret = this->Get_ylow() + this->Get_U()*sin(this->Get_rho());
    return ret;
}

void EC_geometry::Put_UVW(double u, double v, double w){
    this->EC_U = u;
    this->EC_V = v;
    this->EC_W = w;
}

// check that the U-view is a positive number (ie has been filled)
bool EC_geometry::Check_U()
{
    bool ret = (this->Get_U() >= 0) ? true : false;
    return ret;
}

// check that the V-view is a positive number (ie has been filled)
bool EC_geometry::Check_V()
{
    bool ret = (this->Get_V() >= 0) ? true : false;
    return ret;
}

// check that the W-view is a positive number (ie has been filled)
bool EC_geometry::Check_W()
{
    bool ret = (this->Get_W() >= 0) ? true : false;
    return ret;
}

// print the cut information
void EC_geometry::Print_EC_geometry()
{
    cout<<"EC Geometry Constants"<<endl;
    cout<<"========================="<<endl;

    cout << "EC theta = " << this->Get_EC_theta()*180.0/3.14159 << endl;
    cout << "ylow = " << this->Get_ylow() << " (cm)" << endl;
    cout << "yhi = " << this->Get_yhi() << " (cm)" << endl;
    cout << "rho = " << this->Get_rho()*180.0/3.14159 << " (cm)" << endl;
    cout << endl;
}

class OmegaMixedEvent
{
    int nEvtToMix;
    int EvtOffset;
    vector<string> Label;
    vector<int> Index;
    vector<string> evtLabel;
    vector<int> evtIndex;
    TLorentzVector Photon1[2];
    TLorentzVector Photon2[2];
    TLorentzVector PiPlus[2];
    TLorentzVector PiMinus[2];
    TLorentzVector Pi0[2];
    TLorentzVector Omega[2];
public:
    OmegaMixedEvent();
    int Get_NumberOfEventsToMix(){return nEvtToMix;};
    void Put_NumberOfEventsToMix(int i);
    int Get_OffsetOfEventsToMix(){return EvtOffset;};
    void Put_OffsetOfEventsToMix(int i);
    bool Check_evtIndex(int iME);
    bool Check_Index(int iMethod);
    int Get_nLabel() {return Label.size();};
    int Get_nIndex() {return Index.size();};
    string Get_Label(int num) {return Label[num];};
    int Get_nEvtLabel() {return evtLabel.size();};
    int Get_nEvtIndex() {return evtIndex.size();};
    string Get_evtLabel(int num) {return evtLabel[num];};
    TLorentzVector Get_Photon1(int iME);
    TLorentzVector Get_Photon2(int iME);
    TLorentzVector Get_PiPlus(int iME);
    TLorentzVector Get_PiMinus(int iME);
    TLorentzVector Get_Pi0(int iME);
    TLorentzVector Get_Omega(int iME);
    void Reconstruct_Pi0(int iME);
    void Reconstruct_Omega(int iME);
    void Mix_Omega(int iMethod);
    void Put_Photon1(TLorentzVector V, int iME);
    void Put_Photon2(TLorentzVector V, int iME);
    void Put_PiPlus(TLorentzVector V, int iME);
    void Put_PiMinus(TLorentzVector V, int iME);
    void Put_Pi0(TLorentzVector V, int iME);
    void Put_Omega(TLorentzVector V, int iME);
    void Clear_TLorentzVectors();
    void Print_Info();
};

OmegaMixedEvent::OmegaMixedEvent()
{
    int i;
    
    Label.push_back("Mixing Photon 1");
    Label.push_back("Mixing Photon 2");
    Label.push_back("Mixing pi+");
    Label.push_back("Mixing pi-");
    Label.push_back("Mixing pi0");
    
    Index.push_back(0);
    Index.push_back(1);
    Index.push_back(2);
    Index.push_back(3);
    Index.push_back(4);

    evtLabel.push_back("In Time Event");
    evtLabel.push_back("Out of Time Event");
    
    evtIndex.push_back(0);
    evtIndex.push_back(1);

    for(i=0; i<evtIndex.size(); i++){
        Photon1[i].SetPxPyPzE(0.,0.,0.,0.);
        Photon2[i].SetPxPyPzE(0.,0.,0.,0.);
        Pi0[i].SetPxPyPzE(0.,0.,0.,0.);
        PiPlus[i].SetPxPyPzE(0.,0.,0.,0.);
        PiMinus[i].SetPxPyPzE(0.,0.,0.,0.);
        Omega[i].SetPxPyPzE(0.,0.,0.,0.);
    }
}

void OmegaMixedEvent::Clear_TLorentzVectors()
{
    int i;
    
    for(i=0; i < this->Get_nEvtIndex(); i++){
        this->Photon1[i].SetPxPyPzE(0.,0.,0.,0.);
        this->Photon2[i].SetPxPyPzE(0.,0.,0.,0.);
        this->Pi0[i].SetPxPyPzE(0.,0.,0.,0.);
        this->PiPlus[i].SetPxPyPzE(0.,0.,0.,0.);
        this->PiMinus[i].SetPxPyPzE(0.,0.,0.,0.);
        this->Omega[i].SetPxPyPzE(0.,0.,0.,0.);
    }
}

void OmegaMixedEvent::Put_NumberOfEventsToMix(int i){
    this->nEvtToMix = i;
}

void OmegaMixedEvent::Put_OffsetOfEventsToMix(int i){
    this->EvtOffset = i;
}

// check that the event index is in the correct range
bool OmegaMixedEvent::Check_evtIndex(int iME)
{
    bool ret = (iME>=0 && iME<this->Get_nEvtIndex()) ? true : false;
    
    return ret;
}

// check that the method index is in the correct range
bool OmegaMixedEvent::Check_Index(int iMethod)
{
    bool ret = (iMethod>=0 && iMethod<this->Get_nIndex()) ? true : false;
    
    return ret;
}

void OmegaMixedEvent::Put_Photon1(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->Photon1[iME] = V;
}

void OmegaMixedEvent::Put_Photon2(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->Photon2[iME] = V;
}

void OmegaMixedEvent::Put_PiPlus(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->PiPlus[iME] = V;
}

void OmegaMixedEvent::Put_PiMinus(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->PiMinus[iME] = V;
}

void OmegaMixedEvent::Put_Pi0(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->Pi0[iME] = V;
}

void OmegaMixedEvent::Put_Omega(TLorentzVector V, int iME){
    if(this->Check_evtIndex(iME)) this->Omega[iME] = V;
}

TLorentzVector OmegaMixedEvent::Get_Photon1(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
  
    if(this->Check_evtIndex(iME)) ret = Photon1[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_Photon2(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = Photon2[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_PiPlus(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = PiPlus[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_PiMinus(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = PiMinus[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_Pi0(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = Pi0[iME];
    
    return ret;
}

TLorentzVector OmegaMixedEvent::Get_Omega(int iME){
    
    TLorentzVector ret(0.,0.,0.,0.); // init the return variable
    
    if(this->Check_evtIndex(iME)) ret = Omega[iME];
    
    return ret;
}

void OmegaMixedEvent::Reconstruct_Pi0(int iME){
    if(this->Check_evtIndex(iME)){
        TLorentzVector tempPhoton1 = this->Get_Photon1(iME);
        TLorentzVector tempPhoton2 = this->Get_Photon2(iME);
        
        this->Put_Pi0(tempPhoton1 + tempPhoton2, iME);
    }
}

void OmegaMixedEvent::Reconstruct_Omega(int iME){
    if(this->Check_evtIndex(iME)){
        TLorentzVector tempPhoton1 = this->Get_Photon1(iME);
        TLorentzVector tempPhoton2 = this->Get_Photon2(iME);
        TLorentzVector tempPiPlus = this->Get_PiPlus(iME);
        TLorentzVector tempPiMinus = this->Get_PiMinus(iME);
        
        this->Put_Omega(tempPhoton1 + tempPhoton2 + tempPiPlus + tempPiMinus, iME);
    }
}

void OmegaMixedEvent::Mix_Omega(int iMethod){

    double theta, phi, tempPx, tempPy, tempPz, tempE;
    TLorentzVector tempPhoton1 = this->Get_Photon1(0);
    TLorentzVector tempPhoton2 = this->Get_Photon2(0);
    TLorentzVector tempPiPlus = this->Get_PiPlus(0);
    TLorentzVector tempPiMinus = this->Get_PiMinus(0);
    TLorentzVector tempPi0 = tempPhoton1 + tempPhoton2;
    TLorentzVector A;

/*    cout <<"Before " << iMethod <<endl;
    PrintTLorentzVector(tempPhoton1);
    PrintTLorentzVector(tempPhoton2);
    PrintTLorentzVector(tempPiPlus);
    PrintTLorentzVector(tempPiMinus);
    PrintTLorentzVector(tempPi0);
*/
    switch (iMethod){
        case 0:
            tempPhoton1 = this->Get_Photon1(1);
            tempPi0 = tempPhoton1 + tempPhoton2;
            break;
        case 1:
            tempPhoton2 = this->Get_Photon2(1);
            tempPi0 = tempPhoton1 + tempPhoton2;
            break;
        case 2:
            A = this->Get_PiPlus(1);
            A.SetTheta(tempPiPlus.Theta());
            A.SetPhi(tempPiPlus.Phi());
            tempPiPlus = A;
            break;
        case 3:
            A = this->Get_PiMinus(1);
            A.SetTheta(tempPiMinus.Theta());
            A.SetPhi(tempPiMinus.Phi());
            tempPiMinus = A;
            break;
        case 4:
            A = this->Get_Photon1(1) + this->Get_Photon2(1);
            A.SetTheta(tempPi0.Theta());
            A.SetPhi(tempPi0.Phi());
            tempPi0 = A;
            break;
        default:
            cout << "OmegaMixedEvent::Mix_Omega - incorrect iMethod " << iMethod <<endl;
            exit(0);
            break;
    }
/*    cout <<"After " << iMethod <<endl;
    PrintTLorentzVector(tempPhoton1);
    PrintTLorentzVector(tempPhoton2);
    PrintTLorentzVector(tempPiPlus);
    PrintTLorentzVector(tempPiMinus);
    PrintTLorentzVector(tempPi0);
*/
    this->Put_Pi0(tempPi0, 1);
    this->Put_Omega((tempPi0 + tempPiPlus + tempPiMinus), 1);
}

void OmegaMixedEvent::Print_Info()
{
    int ii;
    cout<<"Mixed Event Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nLabel();ii++){
        cout << "Method " << ii+1 << "\t" << this->Get_Label(ii) << endl;
    }
    for(ii=0;ii<this->Get_nEvtLabel();ii++){
        cout << "Event type " << ii+1 << "\t" << this->Get_evtLabel(ii) << endl;
    }
    cout << endl;
}

int process (string inFile, int MaxEvents, int dEvents, int targMass) {
    int i, ii, j, k, kk;
    
    int Sector_index;
    int Vz_index;
    int BankIndex_part[5];
    
    bool cutPi0Mass;
    bool cutZDiff_ElectronNPion;
    bool cutZDiff_ElectronPPion;
    bool cutZDiff;
    bool cutQSquared;
    bool cutOpAng_ElecPhoton1;
    bool cutOpAng_ElecPhoton2;
    bool cutOpAng_ElecPhoton;
    bool cutW;
    bool cutElecR;
    bool cutBetaPhoton1;
    bool cutBetaPhoton2;
    bool cutBetaPhoton;
    bool cutOmegaMass;
    bool cutOmegaMass_sb;
    
    bool cuts_woPi0Mass;
    bool cuts_woZDiff;
    bool cuts_woQsquared;
    bool cuts_woOpAng_ElecPhoton;
    bool cuts_woW;
    bool cuts_woElecR;
    bool cuts_woBetaPhoton;
    bool cutsAll;
    
    bool cuts_ElecID;
    bool ElecID_All;
    bool ElecID_Mom;
    bool ElecID_ECvsP;
    bool ElecID_ECfid;
    bool ElecID_dtECSC;
    
	double TwoPhotonAngle, elecPhoton1Angle, elecPhoton2Angle;
    double Qsq, nu, Mx, z_fracEnergy, W;
    double sinHalfTheta;
    double partMom;
    double timeEC, timeSC, pathEC, pathSC;
    double dt_ECminusSC[5];
    
    double emECu, emECv, emECw, emECin, emECout, emECtot, emCCnphe, emdt; // variables for electron id cuts
    
    EG2Target myTgt;
    EG2Cuts myCuts;
    OmegaMixedEvent myMixEvt;
    ElectronID myElecID;
    EC_geometry myECgeom;
    
    myMixEvt.Put_NumberOfEventsToMix(1); // add number of mixed event iterations
    myMixEvt.Put_OffsetOfEventsToMix(5); // add offset of the entry number for mixed events
    
    int NUM_MIXING_METHODS = myMixEvt.Get_nLabel(); // number of methods for mixing events
    int NUM_ENTRIES_OFFSET = myMixEvt.Get_NumberOfEventsToMix(); // retreive number of mixed event iterations
    int ENTRIES_OFFSET = myMixEvt.Get_OffsetOfEventsToMix(); // retrieve offset of the entry number for mixed events
    
    myCuts.Print_Cuts();
    myMixEvt.Print_Info();
    myElecID.Print_ElectronID();
    
    TLorentzVector BeamMinusElectron;
    TLorentzVector W_TLV;
    TLorentzVector Mx_TLV;
    TLorentzVector TwoPion;
    TLorentzVector TwoPhoton;
	TLorentzVector Omega;

    TLorentzVector photon1_MixedEvt;
    TLorentzVector photon2_MixedEvt;
    TLorentzVector pPion_MixedEvt;
    TLorentzVector nPion_MixedEvt;
    TLorentzVector TwoPhoton_MixedEvt;
    TLorentzVector Omega_MixedEvt;

    int iMixedEvt;
    double Mass_TwoPhoton_ME[NUM_MIXING_METHODS][NUM_ENTRIES_OFFSET];
    double Mass_Omega_ME[NUM_MIXING_METHODS][NUM_ENTRIES_OFFSET];
    
    TLorentzVector beam(0., 0., BEAM_ENERGY, sqrt(BEAM_ENERGY*BEAM_ENERGY+MASS_ELECTRON*MASS_ELECTRON));
	TLorentzVector target(0., 0., 0., MASS_PROTON);
	TLorentzVector nucleon(0., 0., 0., MASS_PROTON);
    
    TEventReader reader;
    reader.addFile(inFile.c_str());
    int entries = reader.getEntries();
    cout << "Entries: " << entries << endl;

    TEventReader readerMixedEvt;
    readerMixedEvt.addFile(inFile.c_str());
    
    int StopProcess;
    if(MaxEvents){
        StopProcess = MaxEvents;
    }else{
        StopProcess = entries;
    }

    int processed = 0;
    for (processed = 0; processed < StopProcess; processed = processed + 1) {
        
        // Initialize the cuts
        cutPi0Mass = false;
        cutZDiff_ElectronNPion = false;
        cutZDiff_ElectronPPion = false;
        cutZDiff = false;
        cutQSquared = false;
        cutW = false;
        cutBetaPhoton1 = false;
        cutBetaPhoton2 = false;
        cutBetaPhoton = false;
        cutOpAng_ElecPhoton1 = false;
        cutOpAng_ElecPhoton2 = false;
        cutOpAng_ElecPhoton = false;
        cutElecR = false;
        cutOmegaMass = false;
        cutsAll = false;
        cuts_woPi0Mass = false;
        cuts_woZDiff = false;
        cuts_woQsquared = false;
        cuts_woW = false;
        cuts_woBetaPhoton = false;
        cuts_woOpAng_ElecPhoton = false;
        cuts_woElecR = false;
        
        ElecID_All = false;
        ElecID_Mom =false;
        ElecID_ECvsP = false;
        ElecID_ECfid = false;
        ElecID_dtECSC = false;
        
        if (!(processed % dEvents)) cout << "Processed Entries: " << processed << endl;
        if (DEBUG) reader.printEvent();
        
        reader.readEntry(processed);
        
        // get the first electron lorentz vector and vertex
		TLorentzVector elec = reader.getLorentzVector(ID_ELECTRON, 0, MASS_ELECTRON);
		TVector3 elec_vert = reader.getVertex(ID_ELECTRON, 0);
        BankIndex_part[0] = reader.getIndexByPid(ID_ELECTRON, 0);
        
		//TLorentzVector prot = reader.getLorentzVector(ID_PROTON, 0, MASS_PROTON);
		//TVector3 prot_vert = reader.getVertex(ID_PROTON, 0);
        
		TLorentzVector nPion = reader.getLorentzVector(ID_PION_NEG, 0,MASS_PION_CHARGED);
		TVector3 nPion_vert = reader.getVertex(ID_PION_NEG, 0);
        BankIndex_part[1] = reader.getIndexByPid(ID_PION_NEG, 0);
        
		TLorentzVector pPion = reader.getLorentzVector(ID_PION_POS, 0, MASS_PION_CHARGED);
		TVector3 pPion_vert = reader.getVertex(ID_PION_POS, 0);
        BankIndex_part[2] = reader.getIndexByPid(ID_PION_POS, 0);

		TLorentzVector photon1 = reader.getLorentzVector(ID_PHOTON, 0, MASS_PHOTON);
		TVector3 photon1_vert = reader.getVertex(ID_PHOTON, 0);
        BankIndex_part[3] = reader.getIndexByPid(ID_PHOTON, 0);

		TLorentzVector photon2 = reader.getLorentzVector(ID_PHOTON, 1, MASS_PHOTON);
		TVector3 photon2_vert = reader.getVertex(ID_PHOTON, 1);
        BankIndex_part[4] = reader.getIndexByPid(ID_PHOTON, 1);

        myMixEvt.Clear_TLorentzVectors(); // initialize all particle TLorentzVectors to zero in myMixEvt
        
        myMixEvt.Put_Photon1(photon1,0);
        myMixEvt.Put_Photon2(photon2,0);
        myMixEvt.Put_PiPlus(pPion,0);
        myMixEvt.Put_PiMinus(nPion,0);
        myMixEvt.Reconstruct_Pi0(0);
        myMixEvt.Reconstruct_Omega(0);
        
        // determine which target the reaction occurred in
        Vz_index = myTgt.Get_Index(elec_vert.Z());
        
        // fill the target Lorentz vector
        switch(Vz_index){
            case 1: target.SetE(MASS_DEUTERIUM); break;
            case 2: target.SetE(targMass * MASS_PROTON); break;
            default: target.SetE(MASS_PROTON); break;
        }
        
        BeamMinusElectron = beam - elec; // Lorentz Vector Difference between beam and scattered electron
        TwoPion = pPion + nPion; // pion pair Lorentz vector
        TwoPhoton = photon1 + photon2; // Two photon Lorentz vector
		Omega = TwoPion + TwoPhoton; // omega Lorentz vector

/*        cout << "Event pi0 " << endl;
        PrintTLorentzVector(photon1);
        PrintTLorentzVector(photon2);
        PrintTLorentzVector(pPion);
        PrintTLorentzVector(nPion);
        PrintTLorentzVector(TwoPhoton);
        PrintTLorentzVector(myMixEvt.Get_Pi0(0));
*/
        if(NUM_ENTRIES_OFFSET*ENTRIES_OFFSET < entries){
            for(k=0; k<NUM_ENTRIES_OFFSET; k++){
                iMixedEvt = processed + (k+1)*ENTRIES_OFFSET;
                if(iMixedEvt >= entries){
                    iMixedEvt = (k+1)*ENTRIES_OFFSET;
                }
                readerMixedEvt.readEntry(iMixedEvt);
                photon1_MixedEvt = readerMixedEvt.getLorentzVector(ID_PHOTON, 0, MASS_PHOTON); // Photon 1 Lorentz vector from an out-of-time event
                photon2_MixedEvt = readerMixedEvt.getLorentzVector(ID_PHOTON, 1, MASS_PHOTON); // Photon 2 Lorentz vector from an out-of-time event
                nPion_MixedEvt = readerMixedEvt.getLorentzVector(ID_PION_NEG, 0,MASS_PION_CHARGED); // pi+ Lorentz vector from an out-of-time event
                pPion_MixedEvt = readerMixedEvt.getLorentzVector(ID_PION_POS, 0, MASS_PION_CHARGED); // pi- Lorentz vector from an out-of-time event
        
                myMixEvt.Put_Photon1(photon1_MixedEvt,1);
                myMixEvt.Put_Photon2(photon2_MixedEvt,1);
                myMixEvt.Put_PiPlus(pPion_MixedEvt,1);
                myMixEvt.Put_PiMinus(nPion_MixedEvt,1);
                
/*                cout << "Mixed test " << endl;
                PrintTLorentzVector(photon1_MixedEvt);
                PrintTLorentzVector(photon2_MixedEvt);
                PrintTLorentzVector(pPion_MixedEvt);
                PrintTLorentzVector(nPion_MixedEvt);
                PrintTLorentzVector(photon1_MixedEvt+photon2_MixedEvt);
                
                cout << "test 1" <<endl;
                PrintTLorentzVector((photon1_MixedEvt+photon2));
                
                cout << "test 2" <<endl;
                PrintTLorentzVector((photon2_MixedEvt+photon1));
*/
                for(kk=0; kk<NUM_MIXING_METHODS; kk++){
                    myMixEvt.Mix_Omega(kk); // run mixing routine for each method
                    
                    TwoPhoton_MixedEvt = myMixEvt.Get_Pi0(1); // Two photon Lorentz vector from an out-of-time event
                    Mass_TwoPhoton_ME[kk][k] = TwoPhoton_MixedEvt.M();

//                    cout << "Method " << myMixEvt.Get_Label(kk) << endl;
//                    PrintTLorentzVector(TwoPhoton_MixedEvt);

                    Omega_MixedEvt = myMixEvt.Get_Omega(1); // Omega Lorentz vector from an out-of-time event
                    Mass_Omega_ME[kk][k] = Omega_MixedEvt.M();
                }
            }
        }
        
        if(DEBUG) cout <<processed<<"\t"<<Omega.M()<<"\t"<<nPion.M()<<"\t"<<pPion.M()<<endl;

        double elecNPionZVertDiff = elec_vert.Z() - nPion_vert.Z(); // z vertex difference, e- and pi-
		double elecPPionZVertDiff = elec_vert.Z() - pPion_vert.Z(); // z vertex difference, e- and pi+
		double elecPhoton1ZVertDiff = elec_vert.Z() - photon1_vert.Z(); // z vertex difference, e- and photon 1
		double elecPhoton2ZVertDiff = elec_vert.Z() - photon2_vert.Z(); // z vertex difference, e- and photon 2
        
        Qsq = -1.0*BeamMinusElectron.M2(); // electron Q^2
        nu = BeamMinusElectron.E(); // energy transfered to target
        
        Mx_TLV = BeamMinusElectron + nucleon - Omega;
        Mx = Mx_TLV.M(); // reaction Mx
        W_TLV = BeamMinusElectron + nucleon;
        W = W_TLV.M(); // reaction W
        z_fracEnergy = Omega.E()/nu; // fractional energy taken by hadron
        
        // Find the electron sector
        Sector_index = GetSectorByPhi(elec.Phi());
        if(Sector_index){
            elecZVertSector[Sector_index-1]->Fill(elec_vert.Z());
        }else{
            cout << "Error in finding sector. Phi = " << elec.Phi() * TMath::RadToDeg() << endl;
        }
        
        //_________________________________
		// Fill histograms
		q2->Fill(Qsq);
        sinHalfTheta = sin(0.5*elec.Theta()); // sine of one-half the electron scattering angle theta
        q2_VS_theta->Fill(4.0*elec.E()*sinHalfTheta*sinHalfTheta,Qsq);
        
        nu_EnergyTransfer->Fill(nu);
		elecZVert->Fill(elec_vert.Z()); // fill electron z vertex histogram
		elecZVert_VS_Phi->Fill(elec.Phi() * TMath::RadToDeg(),elec_vert.Z()); // fill electron z vertex vs phi histogram
        
        hMx[Vz_index]->Fill(Mx); // histogram for Mx
        hW[Vz_index]->Fill(W); // histogram for W
        z_fracE[Vz_index]->Fill(z_fracEnergy); // histogram for fractional z
        
        // plots of z vertex difference between scattered electron and other decay particle
		ZVertDiff->Fill(elecNPionZVertDiff,1);
		ZVertDiff->Fill(elecPPionZVertDiff,2);
		ZVertDiff->Fill(elecPhoton1ZVertDiff,3);
		ZVertDiff->Fill(elecPhoton2ZVertDiff,4);
        
        // plots of x vs y vertices
        Xvert->Fill(elec_vert.X(), 0);
        Xvert->Fill(nPion_vert.X(), 1);
        Xvert->Fill(pPion_vert.X(), 2);
        Xvert->Fill(photon1_vert.X(), 3);
        Xvert->Fill(photon2_vert.X(), 4);

        Yvert->Fill(elec_vert.Y(), 0);
        Yvert->Fill(nPion_vert.Y(), 1);
        Yvert->Fill(pPion_vert.Y(), 2);
        Yvert->Fill(photon1_vert.Y(), 3);
        Yvert->Fill(photon2_vert.Y(), 4);
        
        Xvert_VS_Yvert[0]->Fill(elec_vert.X(), elec_vert.Y());
		Xvert_VS_Yvert[1]->Fill(nPion_vert.X(), nPion_vert.Y());
		Xvert_VS_Yvert[2]->Fill(pPion_vert.X(), pPion_vert.Y());
		Xvert_VS_Yvert[3]->Fill(photon1_vert.X(), photon1_vert.Y());
		Xvert_VS_Yvert[4]->Fill(photon2_vert.X(), photon2_vert.Y());
        
        // plots of angles theta vs phi
		Theta_VS_Phi[0]->Fill(elec.Theta() * TMath::RadToDeg(), elec.Phi() * TMath::RadToDeg());
		Theta_VS_Phi[1]->Fill(nPion.Theta() * TMath::RadToDeg(), nPion.Phi() * TMath::RadToDeg());
		Theta_VS_Phi[2]->Fill(pPion.Theta() * TMath::RadToDeg(), pPion.Phi() * TMath::RadToDeg());
		Theta_VS_Phi[3]->Fill(photon1.Theta() * TMath::RadToDeg(), photon1.Phi() * TMath::RadToDeg());
		Theta_VS_Phi[4]->Fill(photon2.Theta() * TMath::RadToDeg(), photon2.Phi() * TMath::RadToDeg());
		Theta_VS_Phi[5]->Fill(TwoPhoton.Theta() * TMath::RadToDeg(), TwoPhoton.Phi() * TMath::RadToDeg());
		Theta_VS_Phi[6]->Fill(Omega.Theta() * TMath::RadToDeg(), Omega.Phi() * TMath::RadToDeg());

        // plots of total momentum
		TotalMomentum->Fill(elec.P(),0);
		TotalMomentum->Fill(nPion.P(),1);
		TotalMomentum->Fill(pPion.P(),2);
		TotalMomentum->Fill(photon1.P(),3);
		TotalMomentum->Fill(photon2.P(),4);
		TotalMomentum->Fill(TwoPhoton.P(),5);
		TotalMomentum->Fill(Omega.P(),6);

		// plots of beta vs momentum
		Beta_VS_Momentum->Fill(elec.P(), elec.Beta());
		Beta_VS_Momentum->Fill(nPion.P(), nPion.Beta());
		Beta_VS_Momentum->Fill(pPion.P(), pPion.Beta());
		Beta_VS_Momentum->Fill(photon1.P(), photon1.Beta());
		Beta_VS_Momentum->Fill(photon2.P(), photon2.Beta());
     
        //
        // Start of  Electron ID
        //
        for(ii=0; ii<5; ii++){
            switch (ii) {
                case 0: partMom = elec.P(); break;
                case 1: partMom = nPion.P(); break;
                case 2: partMom = pPion.P(); break;
                case 3: partMom = photon1.P(); break;
                case 4: partMom = photon2.P(); break;
                default: partMom = -1.0; break;
            }
        
            CCnphe->Fill(reader.getProperty("ccnphe",BankIndex_part[ii]),ii);
            ECu->Fill(reader.getProperty("ecu",BankIndex_part[ii]),ii);
            ECv->Fill(reader.getProperty("ecv",BankIndex_part[ii]),ii);
            ECw->Fill(reader.getProperty("ecw",BankIndex_part[ii]),ii);
            ECtot_VS_P[ii]->Fill(partMom,reader.getProperty("ectot",BankIndex_part[ii]));
            ECtotP_VS_P[ii]->Fill(partMom,reader.getProperty("ectot",BankIndex_part[ii])/partMom);
            ECin_VS_ECout[ii]->Fill(reader.getProperty("ecin",BankIndex_part[ii]),reader.getProperty("ecout",BankIndex_part[ii]));
            
            timeEC = reader.getProperty("ectime",BankIndex_part[ii]);
            timeSC = reader.getProperty("sctime",BankIndex_part[ii]);
            pathEC = reader.getProperty("ecpath",BankIndex_part[ii]);
            pathSC = reader.getProperty("scpath",BankIndex_part[ii]);
            dt_ECminusSC[ii] = timeEC - timeSC - 0.7;
            dtime_ECSC->Fill(dt_ECminusSC[ii],ii);
        }

        // Electron ID cuts
        emECtot = reader.getProperty("ectot",BankIndex_part[0]);
        emECin = reader.getProperty("ecin",BankIndex_part[0]);
        emECout = reader.getProperty("ecout",BankIndex_part[0]);
        emECu = reader.getProperty("ecu",BankIndex_part[0]);
        emECv = reader.getProperty("ecv",BankIndex_part[0]);
        emECw = reader.getProperty("ecw",BankIndex_part[0]);
        emCCnphe = reader.getProperty("ccnphe",BankIndex_part[0]);
        emdt = reader.getProperty("ectime",BankIndex_part[0]) - reader.getProperty("sctime",BankIndex_part[0]) - 0.7;
      
        myECgeom.Put_UVW(emECu,emECv,emECw);
        EC_XvsY_local->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
        
        ElecID_Mom = myElecID.Check_ElecMom(elec.P());
        ElecID_ECvsP = myElecID.Check_ElecECoverP(elec.P(),emECtot,Sector_index,targMass);
        ElecID_dtECSC = myElecID.Check_Elec_dtECSC(emdt);
        ElecID_ECfid = myElecID.Check_ElecECu(emECu) && myElecID.Check_ElecECv(emECv) && myElecID.Check_ElecECw(emECw);
        ElecID_All = (ElecID_Mom && ElecID_ECvsP && ElecID_dtECSC && ElecID_ECfid);
        
        if(ElecID_ECfid){
            ECin_VS_ECout_ECfid->Fill(emECin,emECout);
        }
        
        if(ElecID_All){
            ECin_VS_ECout_elecID_All->Fill(emECin,emECout);
        }

        if (emECout < 0.01){
            Beta_VS_Momentum_ECoutCut->Fill(elec.P(), elec.Beta());
            Theta_VS_Phi_ECoutCut->Fill(elec.Theta() * TMath::RadToDeg(), elec.Phi() * TMath::RadToDeg());
            elecZVert_ECoutCut->Fill(elec_vert.Z());
            q2_ECoutCut->Fill(Qsq);
            
            ECtot_VS_P_ECoutCut->Fill(elec.P(),emECtot);
            ECtotP_VS_P_ECoutCut->Fill(elec.P(),emECtot/elec.P());
            ECtotMinusECin_ECoutCut->Fill(emECtot-emECin);
        }
        else {
            Beta_VS_Momentum_AntiECoutCut->Fill(elec.P(), elec.Beta());
            Theta_VS_Phi_AntiECoutCut->Fill(elec.Theta() * TMath::RadToDeg(), elec.Phi() * TMath::RadToDeg());
            elecZVert_AntiECoutCut->Fill(elec_vert.Z());
            q2_AntiECoutCut->Fill(Qsq);
            
            ECtot_VS_P_AntiECoutCut->Fill(elec.P(),emECtot);
            ECtotP_VS_P_AntiECoutCut->Fill(elec.P(),emECtot/elec.P());
            ECtotMinusECin_AntiECoutCut->Fill(emECtot-emECin);
        }
        
        // Testing the electron ID
        for(ii=0; ii<myElecID.Get_nElecID(); ii++){
            cuts_ElecID = false; // intialize the cuts
            
            if (myElecID.Get_elecIDLabel(ii).compare("No cuts")==0) {
                cuts_ElecID = true;
            }else if (myElecID.Get_elecIDLabel(ii).compare("Momentum")==0) {
                cuts_ElecID = myElecID.Check_ElecMom(elec.P());
            }else if (myElecID.Get_elecIDLabel(ii).compare("EC U-view")==0) {
                cuts_ElecID = myElecID.Check_ElecECu(emECu);
            }else if (myElecID.Get_elecIDLabel(ii).compare("EC V-view")==0) {
                cuts_ElecID = myElecID.Check_ElecECv(emECv);
            }else if (myElecID.Get_elecIDLabel(ii).compare("EC W-view")==0) {
                cuts_ElecID = myElecID.Check_ElecECw(emECw);
            }else if (myElecID.Get_elecIDLabel(ii).compare("ECin")==0) {
                cuts_ElecID = myElecID.Check_ElecECin(emECin);
            }else if (myElecID.Get_elecIDLabel(ii).compare("CC Nphe")==0) {
                cuts_ElecID = myElecID.Check_ElecCCnphe(emCCnphe);
            }else if (myElecID.Get_elecIDLabel(ii).compare("dt(EC-SC)")==0) {
                cuts_ElecID = myElecID.Check_Elec_dtECSC(emdt);
            }else if (myElecID.Get_elecIDLabel(ii).compare("ECtot/P VS P")==0) {
                cuts_ElecID = myElecID.Check_ElecECoverP(elec.P(),emECtot,Sector_index,targMass);
            }else{
                cuts_ElecID = false;
            }
            
            if(cuts_ElecID){
                CCnphe_elecID->Fill(emCCnphe,ii);
                Mom_elecID->Fill(elec.P(),ii);
                ECu_elecID->Fill(emECu,ii);
                ECv_elecID->Fill(emECv,ii);
                ECw_elecID->Fill(emECw,ii);
                dtime_ECSC_elecID->Fill(emdt,ii);
                ECtot_VS_P_elecID[ii]->Fill(elec.P(),emECtot);
                ECtotP_VS_P_elecID[ii]->Fill(elec.P(),emECtot/elec.P());
                ECin_VS_ECout_elecID[ii]->Fill(emECin,emECout);
                Mom_VS_ECout_elecID[ii]->Fill(elec.P(),emECout);
                ECu_VS_ECout_elecID[ii]->Fill(emECu,emECout);
                ECv_VS_ECout_elecID[ii]->Fill(emECv,emECout);
                ECw_VS_ECout_elecID[ii]->Fill(emECw,emECout);
            }
        }
        //
        // End of  Electron ID
        //
        
        // plot of two photon opening angle
		TwoPhotonAngle = TMath::RadToDeg()*photon1.Angle(photon2.Vect());
		OpAng_2Photons->Fill(TwoPhotonAngle);

		elecPhoton1Angle = TMath::RadToDeg()*elec.Angle(photon1.Vect());
        OpAng_elecPhoton1->Fill(elecPhoton1Angle);

		elecPhoton2Angle = TMath::RadToDeg()*elec.Angle(photon2.Vect());
        OpAng_elecPhoton2->Fill(elecPhoton2Angle);
        
        // plots by target (Vz_index)
		LongMom[Vz_index]->Fill(Omega.P()-Omega.Pt()); // omega long. mom.
		TransMom[Vz_index]->Fill(Omega.Pt()); // omega trans. mom.
        
		OpAng_VS_IM2Photons[Vz_index]->Fill(TwoPhoton.M(),TwoPhotonAngle); // opening angle vs 2 photon inv. mass
		OpAng_VS_E[Vz_index]->Fill(TwoPhoton.E(),TwoPhotonAngle); // opening angle vs 2 photon total energy
        
		MissMom[Vz_index]->Fill((beam + target - Omega).P());  // mising mom.
		MMsq[Vz_index]->Fill((beam + target - Omega).M2()); // missing mass^2
        
        // plots of variable vs the omega inv. mass
		IM2Pions_VS_IMOmega[Vz_index]->Fill(TwoPion.M(), Omega.M()); // variable = pion pair inv. mass
        IM2Photons_VS_IMOmega[Vz_index]->Fill(TwoPhoton.M(), Omega.M()); // variable = 2 photon inv. mass
		Q2_VS_IMOmega[Vz_index]->Fill(Qsq, Omega.M()); // variable = Q^2
		Pt_VS_IMOmega[Vz_index]->Fill(Omega.Pt(), Omega.M()); // variable = omega trans. mom.
		Pl_VS_IMOmega[Vz_index]->Fill(Omega.Pz(), Omega.M()); // variable = omega long. mom.
		OpAng_VS_IMOmega[Vz_index]->Fill(Omega.M(), TwoPhotonAngle); // variable = 2 photon opening angle

        // set the cuts
        cutPi0Mass = myCuts.Check_MassPi0(TwoPhoton.M()); // pi0 mass cut
        cutQSquared = myCuts.Check_QSquared(Qsq); // Q^2 cut
        cutW = myCuts.Check_Wcut(W); // W cut
        
        cutZDiff_ElectronNPion = myCuts.Check_Zdiff_ElecPim(elecNPionZVertDiff); // e-,pi- z-vertex matching cut
        cutZDiff_ElectronPPion = myCuts.Check_Zdiff_ElecPip(elecPPionZVertDiff); // e-,pi+ z-vertex matching cut
        cutZDiff = (cutZDiff_ElectronNPion && cutZDiff_ElectronPPion); // final e-,pion z-vertex matching cut

        cutOpAng_ElecPhoton1 = myCuts.Check_OpAng_ElecPhoton(elecPhoton1Angle); // e-,photon 1 opening angle cut
        cutOpAng_ElecPhoton2 = myCuts.Check_OpAng_ElecPhoton(elecPhoton2Angle); // e-,photon 2 opening angle cut
        cutOpAng_ElecPhoton = (cutOpAng_ElecPhoton1 && cutOpAng_ElecPhoton2); // final e-,photon opening angle cut
        
        cutBetaPhoton1 = myCuts.Check_BetaPhoton(photon1.Beta()); // photon 1 beta cut
        cutBetaPhoton2 = myCuts.Check_BetaPhoton(photon2.Beta()); // photon 2 beta cut
        cutBetaPhoton = (cutBetaPhoton1 && cutBetaPhoton2); // final photon beta cut
        
        cutElecR = myCuts.Check_ElectronR(elec_vert.Perp()); // electron vertex radius (x,y)
        
        cutOmegaMass = myCuts.Check_MassOmega(Omega.M()); // omega mass cut
        cutOmegaMass_sb = myCuts.Check_MassOmega_sb(Omega.M()); // sideband cuts on the omega mass
        
        // omega inv. mass histograms
        IM2Photons[Vz_index]->Fill(TwoPhoton.M(),0); // inv. mass of 2 photons
		IMOmega[Vz_index]->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons
		IMOmega_woCut[Vz_index]->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons
        IMOmega_antiCut[Vz_index]->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons
        
        if(cutPi0Mass) { // applying the pi0 mass cut
            OpAng_VS_E_MassPi0Cut[Vz_index]->Fill(TwoPhoton.E(),TwoPhotonAngle);
			IM2Photons[Vz_index]->Fill(TwoPhoton.M(),1);
            IMOmega[Vz_index]->Fill(Omega.M(),1);
        }else{
            IMOmega_antiCut[Vz_index]->Fill(Omega.M(),1);
        }
        cuts_woPi0Mass = (cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutBetaPhoton && cutW);
        if(cuts_woPi0Mass){ // applying all cuts except the pi0 mass
            IMOmega_woCut[Vz_index]->Fill(Omega.M(),1);
            IM2Photons[Vz_index]->Fill(TwoPhoton.M(),7);
        }
        
        if(cutQSquared) { // applying the Q^2 cut
            IM2Photons[Vz_index]->Fill(TwoPhoton.M(),2);
            IMOmega[Vz_index]->Fill(Omega.M(),2);
        }else{
            IMOmega_antiCut[Vz_index]->Fill(Omega.M(),2);
        }
        cuts_woQsquared = (cutPi0Mass && cutZDiff && cutOpAng_ElecPhoton && cutBetaPhoton && cutW);
        if(cuts_woQsquared){ // applying all cuts except Q^2
            IMOmega_woCut[Vz_index]->Fill(Omega.M(),2);
        }
        
        if(cutW){ // applying the W cut
			IM2Photons[Vz_index]->Fill(TwoPhoton.M(),3);
            IMOmega[Vz_index]->Fill(Omega.M(),3);
        }else{
            IMOmega_antiCut[Vz_index]->Fill(Omega.M(),3);
        }
        cuts_woW = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutBetaPhoton);
        if(cuts_woW){ // applying all cuts except W
            IMOmega_woCut[Vz_index]->Fill(Omega.M(),3);
        }
        
        if(cutZDiff){ // applying the e-,pion z-vertex matching cut
			IM2Photons[Vz_index]->Fill(TwoPhoton.M(),4);
            IMOmega[Vz_index]->Fill(Omega.M(),4);
        }else{
            IMOmega_antiCut[Vz_index]->Fill(Omega.M(),4);
        }
        cuts_woZDiff = (cutPi0Mass && cutQSquared && cutOpAng_ElecPhoton && cutBetaPhoton && cutW);
        if(cuts_woZDiff){ // applying all cuts except the e-,pion z-vertex matching
            IMOmega_woCut[Vz_index]->Fill(Omega.M(),4);
        }
        
        if(cutOpAng_ElecPhoton) { // applying the e-,photon opening angle cut
			IM2Photons[Vz_index]->Fill(TwoPhoton.M(),5);
            IMOmega[Vz_index]->Fill(Omega.M(),5);
        }else{
            IMOmega_antiCut[Vz_index]->Fill(Omega.M(),5);
        }
        cuts_woOpAng_ElecPhoton = (cutPi0Mass && cutZDiff && cutQSquared && cutBetaPhoton && cutW);
        if(cuts_woOpAng_ElecPhoton){ // applying all cuts except the e-,photon opening angle
            IMOmega_woCut[Vz_index]->Fill(Omega.M(),5);
        }
        
        if(cutBetaPhoton){ // applying the photon beta cut
			IM2Photons[Vz_index]->Fill(TwoPhoton.M(),6);
            IMOmega[Vz_index]->Fill(Omega.M(),6);
        }else{
            IMOmega_antiCut[Vz_index]->Fill(Omega.M(),6);
        }
        cuts_woBetaPhoton = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutW);
        if(cuts_woBetaPhoton){ // applying all cuts except the photon beta
            IMOmega_woCut[Vz_index]->Fill(Omega.M(),6);
        }

        if(cutElecR){ // applying the photon beta cut
            IM2Photons[Vz_index]->Fill(TwoPhoton.M(),7);
            IMOmega[Vz_index]->Fill(Omega.M(),7);
        }else{
            IMOmega_antiCut[Vz_index]->Fill(Omega.M(),7);
        }
        cuts_woElecR = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutW && cutBetaPhoton);
        if(cuts_woElecR){ // applying all cuts except the photon beta
            IMOmega_woCut[Vz_index]->Fill(Omega.M(),7);
        }
        
         // applying all cuts
        cutsAll = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutBetaPhoton && cutW);
        if(cutsAll){
            IMOmega[Vz_index]->Fill(Omega.M(),8);
            W_VS_IMOmega_AllCuts[Vz_index]->Fill(W, Omega.M()); // variable = W
            IM2Pions_VS_IMOmega_AllCuts[Vz_index]->Fill(TwoPion.M(), Omega.M()); // variable = pion pair inv. mass
            PtSq_Omega_AllCuts[Vz_index]->Fill(Omega.Perp2());

            Xvert_VS_Yvert_AllCuts[Vz_index]->Fill(elec_vert.X(), elec_vert.Y());

            if(cutOmegaMass){
                Xvert_VS_Yvert_Omega[Vz_index]->Fill(elec_vert.X(), elec_vert.Y());
                PtSq_Omega_AllCuts_IMOmegaCut[Vz_index]->Fill(Omega.Perp2());
            }
            if(cutOmegaMass_sb){
                PtSq_Omega_AllCuts_IMOmegaSBCut[Vz_index]->Fill(Omega.Perp2());
            }
        }else{
            IMOmega_antiCut[Vz_index]->Fill(Omega.M(),8);
        }
    
        for(k=0; k<NUM_ENTRIES_OFFSET; k++){ // loop over number of mixed event iterations
            for(j=0; j<NUM_MIXING_METHODS; j++){ // loop over number of mixed event methods
                IM2Photons_ME[Vz_index]->Fill(Mass_TwoPhoton_ME[j][k],j); // inv. mass of 2 photons
                IMOmega_ME[Vz_index]->Fill(Mass_Omega_ME[j][k],j); // inv. mass of pi+ pi- 2 photons
                if(cutOpAng_ElecPhoton1 && cutOpAng_ElecPhoton2) {
                    IM2Photons_OpAng_ElecPhoton_Cut_ME[Vz_index]->Fill(Mass_TwoPhoton_ME[j][k],j);
                    IMOmega_OpAng_ElecPhoton_Cut_ME[Vz_index]->Fill(Mass_Omega_ME[j][k],j);
                }
                if(cutPi0Mass) {
                    IMOmega_MassPi0Cut_ME[Vz_index]->Fill(Mass_Omega_ME[j][k],j);
                }
                if(cutZDiff_ElectronNPion && cutZDiff_ElectronPPion){
                    IMOmega_ZVertCut_ME[Vz_index]->Fill(Mass_Omega_ME[j][k],j);
                }
                if(cutQSquared) {
                    IMOmega_QsqCut_ME[Vz_index]->Fill(Mass_Omega_ME[j][k],j);
                }
                if(cutsAll){
                    IMOmega_AllCuts_ME[Vz_index]->Fill(Mass_Omega_ME[j][k],j);
                }
            }
        }
        
        //-----------------------------------------------------
        // plots to check special relativistic kinematics
		TLorentzVector pi0(0, 0, TwoPhoton.Pz(), TMath::Sqrt(MASS_PION_NEUTRAL*MASS_PION_NEUTRAL + TwoPhoton.Pz() * TwoPhoton.Pz()));
		BetaPi0->Fill(pi0.Beta());
		GammaPi0->Fill(pi0.Gamma());
        
        TLorentzRotation lbr(pi0.BoostVector());
        
		TLorentzVector photon1_pi0Rest_caseA(0, 0, 0.5 * MASS_PION_NEUTRAL, 0.5 * MASS_PION_NEUTRAL);
        TLorentzVector photon2_pi0Rest_caseA(0, 0, -0.5 * MASS_PION_NEUTRAL, 0.5 * MASS_PION_NEUTRAL);
        photon1_pi0Rest_caseA.Transform(lbr);
        photon2_pi0Rest_caseA.Transform(lbr);
        
        TLorentzVector photon1_pi0Rest_caseB(0, 0.5 * MASS_PION_NEUTRAL, 0, 0.5 * MASS_PION_NEUTRAL);
        TLorentzVector photon2_pi0Rest_caseB(0, -0.5 * MASS_PION_NEUTRAL, 0, 0.5 * MASS_PION_NEUTRAL);
        photon1_pi0Rest_caseB.Transform(lbr);
        photon2_pi0Rest_caseB.Transform(lbr);
        
        RelativityOpAngPhotonsA->Fill(pi0.Pz(),photon1_pi0Rest_caseA.Angle(photon2_pi0Rest_caseA.Vect())*TMath::RadToDeg());
        RelativityOpAngPhotonsB->Fill(pi0.Pz(),photon1_pi0Rest_caseB.Angle(photon2_pi0Rest_caseB.Vect())*TMath::RadToDeg());
        
        //-----------------------------------------------------
    }
    
    return processed;
}

// Return the CLAS sector number from the azimuthal angle phi
//
// Angle phi must be given in radians
//
int GetSectorByPhi(Double_t phi_rad){
    
    Int_t ret = 0; // init the return variable
    Double_t phi_deg = phi_rad * TMath::RadToDeg(); // convert to degrees
    
    if(CheckCut(phi_deg,-30,30)) {
        ret = 1;
    } else if(CheckCut(phi_deg,-90,-30)) {
        ret = 2;
    } else if(CheckCut(phi_deg,-150,-90)) {
        ret = 3;
    } else if(phi_deg > 150 || phi_deg < -150) {
        ret = 4;
    } else if(CheckCut(phi_deg,90,150)) {
        ret = 5;
    } else if(CheckCut(phi_deg,30,90)) {
        ret = 6;
    }
    
    return ret;
}

//
// CheckCut - return 1 for true and 0 for false cut
//
//				x = variable
//				LowerLimit = lower limit on the range
//				UpperLimit = upper limit on the range
//
int CheckCut(double var, double LowerLimit, double UpperLimit)
{
	int ret = (var >= LowerLimit && var < UpperLimit) ? 1 : 0;
	
	return ret;
}

void PrintUsage(char *processName)
{
    cerr << processName << " <options> <filename>\n";
    cerr << "\toptions are:\n";
    cerr << "\t-o<filename>\tROOT output file (def. = K0from2Pions.root).\n";
    cerr << "\t-M#\t\tprocess maximum # of events.\n";
    cerr << "\t-D#\t\tinform user when # of events have been processed (def. = 1000).\n";
    cerr << "\t-T#\t\tTarget mass number\n";
    cerr << "\t-i\t\tquiet mode (no counter).\n";
    cerr << "\t-h\t\tprint the above" << endl;
}


void PrintAnalysisTime(float tStart, float tStop){
    //time to complete function
    float minutes = 0;
    float seconds = 0;
    minutes = (tStop - tStart)/1000000;
    minutes = (minutes)/60;
    seconds = fmod(minutes,1);
    minutes = minutes-seconds;
    seconds = seconds*60;
    
    if (minutes==0){
        cout<<endl<<"Completed in "<<seconds<<" seconds."<<endl<<endl;
    }
    else{
        cout<<endl<<"Completed in "<<minutes<<" minutes and "<<seconds<<" seconds."<<endl<<endl;
    }
}

//
// BookHist - routine to set up histograms
//
void BookHist(){

    int i, j;
    
    char hname[100];
	char htitle[100];
    
    int nIMomega = 125;
    double IMomegaLo = 0.0;
    double IMomegaHi = 2.5;

    int nPtSq_omega = 100;
    double PtSq_omegaLo = 0.0;
    double PtSq_omegaHi = 2.0;
    
    DetectedParticles myDetPart;
    ParticleList myPartList;
    EG2Target myTgt;
    OmegaMixedEvent myMixEvt;
    ElectronID myElecID;
    
    int nME_Methods = myMixEvt.Get_nLabel();
    int nElecID = myElecID.Get_nElecID();
    
    sprintf(hname,"q2");
    sprintf(htitle,"Q^{2}");
    q2 = new TH1D(hname,htitle, 100, 0., 4.);

    sprintf(hname,"q2_VS_theta");
    sprintf(htitle,"Q^{2} vs. 4E_{e'}sin^{2}(0.5*#theta_{e'})");
    q2_VS_theta = new TH2D(hname,htitle, 200, 0., 1.0, 200, 0., 4.);
    
    sprintf(hname,"nu_EnergyTransfer");
    sprintf(htitle,"\nu");
    nu_EnergyTransfer = new TH1D(hname,htitle, 100, 0., 5.);

    sprintf(hname,"elecZVert");
    sprintf(htitle,"Z Vertex of Electron");
	elecZVert = new TH1D(hname,htitle, 300, -35, -20);

    sprintf(hname,"elecZVert_VS_Phi");
    sprintf(htitle,"Z Vertex  vs. #phi, Electrons");
    elecZVert_VS_Phi = new TH2D(hname,htitle, 360, -180., 180., 300, -35., -20.);

    sprintf(hname,"Xvert");
    sprintf(htitle,"X vertex");
    Xvert = new TH2D(hname,htitle, 300, -10, 10,5,-0.5,4.5);

    sprintf(hname,"Yvert");
    sprintf(htitle,"Y vertex");
    Yvert = new TH2D(hname,htitle, 300, -10, 10,5,-0.5,4.5);
    
    sprintf(hname,"ZVertDiff");
    sprintf(htitle,"Difference Between Z Vertices of electron and other particle");
    ZVertDiff = new TH2D(hname,htitle, 300, -10, 10,4,0.5,4.5);
    
    sprintf(hname,"Beta_VS_Momentum");
    sprintf(htitle,"Beta vs Momentum");
	Beta_VS_Momentum = new TH2D(hname,htitle, 500, 0, 5, 115, 0, 1.15);

    sprintf(hname,"TotalMomentum");
    sprintf(htitle,"Total Momentum");
	TotalMomentum = new TH2D(hname,htitle, 600, 0, 6, 7, -0.5, 6.5);
    
    sprintf(hname,"OpAng_2Photons");
    sprintf(htitle,"Opening Angle Between Photons");
	OpAng_2Photons = new TH1D(hname,htitle, 180, 0, 180);

    sprintf(hname,"OpAng_elecPhoton1");
    sprintf(htitle,"Opening Angle Between e^{-} and #gamma_{1}");
	OpAng_elecPhoton1 = new TH1D(hname,htitle, 180, 0, 180);

    sprintf(hname,"OpAng_elecPhoton2");
    sprintf(htitle,"Opening Angle Between e^{-} and #gamma_{2}");
	OpAng_elecPhoton2 = new TH1D(hname,htitle, 180, 0, 180);
    
    // particle ID histogram
    sprintf(hname,"CCnphe");
    sprintf(htitle,"CC Number of Photo-electrons");
    CCnphe = new TH2D(hname,htitle, 100, 0, 100, 5, -0.5, 4.5);

    sprintf(hname,"ECu");
    sprintf(htitle,"EC U-view");
    ECu = new TH2D(hname,htitle, 450, 0, 450, 5, -0.5, 4.5);
    
    sprintf(hname,"ECv");
    sprintf(htitle,"EC V-view");
    ECv = new TH2D(hname,htitle, 450, 0, 450, 5, -0.5, 4.5);
    
    sprintf(hname,"ECw");
    sprintf(htitle,"EC W-view");
    ECw = new TH2D(hname,htitle, 450, 0, 450, 5, -0.5, 4.5);
    
    sprintf(hname,"dtime_ECSC");
    sprintf(htitle,"#Delta t(EC-SC)");
    dtime_ECSC = new TH2D(hname,htitle, 100, -5.0, 5.0, 5, -0.5, 4.5);
    
    for(i=0; i<myPartList.Get_nPartLabel(); i++){
        sprintf(hname,"Theta_VS_Phi_%s",myPartList.Get_PartLabel(i).c_str());
        sprintf(htitle,"Theta vs Phi for %s",myPartList.Get_PartLabel(i).c_str());
        Theta_VS_Phi[i] = new TH2D(hname,htitle, 180, 0, 180, 360, -180, 180);
    }
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        sprintf(hname,"Xvert_VS_Yvert_%s",myDetPart.Get_DetPartLabel(i).c_str());
        sprintf(htitle,"X Vertex vs Y Vertex, %s",myDetPart.Get_DetPartLabel(i).c_str());
        Xvert_VS_Yvert[i] = new TH2D(hname,htitle, 100, -0.05, 0.05, 100, -0.05, 0.05);

    	sprintf(hname,"ECtot_VS_P_%s",myDetPart.Get_DetPartLabel(i).c_str());
    	sprintf(htitle,"ECtot vs P, %s",myDetPart.Get_DetPartLabel(i).c_str());
    	ECtot_VS_P[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 1.0);

    	sprintf(hname,"ECin_VS_ECout_%s",myDetPart.Get_DetPartLabel(i).c_str());
    	sprintf(htitle,"ECin vs ECout, %s",myDetPart.Get_DetPartLabel(i).c_str());
    	ECin_VS_ECout[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.25);

        sprintf(hname,"ECtotP_VS_P_%s",myDetPart.Get_DetPartLabel(i).c_str());
        sprintf(htitle,"ECtot/P vs P, %s",myDetPart.Get_DetPartLabel(i).c_str());
        ECtotP_VS_P[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);
    }

    // electron ID histogram
    sprintf(hname,"Mom_elecID");
    sprintf(htitle,"Momentum");
    Mom_elecID = new TH2D(hname,htitle, 500, 0, 5.0, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"CCnphe_elecID");
    sprintf(htitle,"CC Number of Photo-electrons");
    CCnphe_elecID = new TH2D(hname,htitle, 100, 0, 100, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"ECu_elecID");
    sprintf(htitle,"EC U-view");
    ECu_elecID = new TH2D(hname,htitle, 450, 0, 450, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"ECv_elecID");
    sprintf(htitle,"EC V-view");
    ECv_elecID = new TH2D(hname,htitle, 450, 0, 450, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"ECw_elecID");
    sprintf(htitle,"EC W-view");
    ECw_elecID = new TH2D(hname,htitle, 450, 0, 450, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"dtime_ECSC_elecID");
    sprintf(htitle,"#Delta t(EC-SC)");
    dtime_ECSC_elecID = new TH2D(hname,htitle, 100, -5.0, 5.0, nElecID, -0.5, nElecID - 0.5);
    
    for(i=0; i<myElecID.Get_nElecID(); i++){
        sprintf(hname,"ECtot_VS_P_elecID_0%i",i);
        sprintf(htitle,"ECtot vs P, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECtot_VS_P_elecID[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 1.0);
        
        sprintf(hname,"ECin_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"ECin vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECin_VS_ECout_elecID[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.25);
        
        sprintf(hname,"ECtotP_VS_P_elecID_0%i",i);
        sprintf(htitle,"ECtot/P vs P, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECtotP_VS_P_elecID[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);

        sprintf(hname,"Mom_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"Mom. vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        Mom_VS_ECout_elecID[i] = new TH2D(hname,htitle, 500, 0, 5.0, 100, 0, 0.25);

        sprintf(hname,"ECu_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"EC U-view vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECu_VS_ECout_elecID[i] = new TH2D(hname,htitle, 450, 0, 450, 100, 0, 0.25);

        sprintf(hname,"ECv_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"EC V-view vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECv_VS_ECout_elecID[i] = new TH2D(hname,htitle, 450, 0, 450, 100, 0, 0.25);
        
        sprintf(hname,"ECw_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"EC W-view vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECw_VS_ECout_elecID[i] = new TH2D(hname,htitle, 450, 0, 450, 100, 0, 0.25);
    }
    
    sprintf(hname,"ECin_VS_ECout_elecID_All");
    sprintf(htitle,"ECin vs ECout, all e- cuts");
    ECin_VS_ECout_elecID_All = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.25);
    
    sprintf(hname,"Theta_VS_Phi_ECoutCut");
    sprintf(htitle,"Theta vs Phi for e-, EC_{out} < 0.01 GeV");
    Theta_VS_Phi_ECoutCut = new TH2D(hname,htitle, 80, 0, 80, 360, -180, 180);
    
    sprintf(hname,"q2_ECoutCut");
    sprintf(htitle,"Q^{2}, EC_{out} < 0.01 GeV");
    q2_ECoutCut = new TH1D(hname,htitle, 100, 0., 4.);
    
    sprintf(hname,"elecZVert_ECoutCut");
    sprintf(htitle,"Z Vertex of Electron, EC_{out} < 0.01 GeV");
    elecZVert_ECoutCut = new TH1D(hname,htitle, 300, -35, -20);
    
    sprintf(hname,"Beta_VS_Momentum_ECoutCut");
    sprintf(htitle,"Beta vs Momentum, EC_{out} < 0.01 GeV");
    Beta_VS_Momentum_ECoutCut = new TH2D(hname,htitle, 500, 0, 5, 115, 0, 1.15);
    
    sprintf(hname,"ECtot_VS_P_ECoutCut");
    sprintf(htitle,"ECtot vs P, EC_{out} < 0.01 GeV");
    ECtot_VS_P_ECoutCut = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 1.0);
    
    sprintf(hname,"ECtotP_VS_P_ECoutCut");
    sprintf(htitle,"ECtot/P vs P, EC_{out} < 0.01 GeV");
    ECtotP_VS_P_ECoutCut = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);

    sprintf(hname,"ECtotMinusECin_ECoutCut");
    sprintf(htitle,"ECtot - ECin, EC_{out} < 0.01 GeV");
    ECtotMinusECin_ECoutCut = new TH1D(hname,htitle, 100,-1.0,1.0);
    
    sprintf(hname,"Theta_VS_Phi_AntiECoutCut");
    sprintf(htitle,"Theta vs Phi for e-, EC_{out} >= 0.01 GeV");
    Theta_VS_Phi_AntiECoutCut = new TH2D(hname,htitle, 80, 0, 80, 360, -180, 180);
    
    sprintf(hname,"q2_AntiECoutCut");
    sprintf(htitle,"Q^{2}, EC_{out} >= 0.01 GeV");
    q2_AntiECoutCut = new TH1D(hname,htitle, 100, 0., 4.);
    
    sprintf(hname,"elecZVert_AntiECoutCut");
    sprintf(htitle,"Z Vertex of Electron, EC_{out} >= 0.01 GeV");
    elecZVert_AntiECoutCut = new TH1D(hname,htitle, 300, -35, -20);
    
    sprintf(hname,"Beta_VS_Momentum_AntiECoutCut");
    sprintf(htitle,"Beta vs Momentum, EC_{out} >= 0.01 GeV");
    Beta_VS_Momentum_AntiECoutCut = new TH2D(hname,htitle, 500, 0, 5, 115, 0, 1.15);
    
    sprintf(hname,"ECtot_VS_P_AntiECoutCut");
    sprintf(htitle,"ECtot vs P, EC_{out} >= 0.01 GeV");
    ECtot_VS_P_AntiECoutCut = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 1.0);
    
    sprintf(hname,"ECtotP_VS_P_AntiECoutCut");
    sprintf(htitle,"ECtot/P vs P, EC_{out} >= 0.01 GeV");
    ECtotP_VS_P_AntiECoutCut = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);
    
    sprintf(hname,"ECtotMinusECin_AntiECoutCut");
    sprintf(htitle,"ECtot - ECin, EC_{out} >= 0.01 GeV");
    ECtotMinusECin_AntiECoutCut = new TH1D(hname,htitle, 100,-1.0,1.0);
    
    sprintf(hname,"ECin_VS_ECout_ECfid");
    sprintf(htitle,"ECin vs ECout, EC fid. cut");
    ECin_VS_ECout_ECfid = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.25);
    
    sprintf(hname,"EC_XvsY_local");
    sprintf(htitle,"EC local X vs local Y");
    EC_XvsY_local = new TH2D(hname,htitle, 100, -100, 100, 100, 0, -100,100);
    
    for(i=0; i<myTgt.Get_nIndex(); i++){
        sprintf(hname,"Xvert_VS_Yvert_AllCuts_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"X Vertex vs Y Vertex, All Cuts, %s",myTgt.Get_Label(i).c_str());
        Xvert_VS_Yvert_AllCuts[i] = new TH2D(hname,htitle, 100, -0.05, 0.05, 100, -0.05, 0.05);

        sprintf(hname,"Xvert_VS_Yvert_Omega_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"X Vertex vs Y Vertex, #omega, %s",myTgt.Get_Label(i).c_str());
        Xvert_VS_Yvert_Omega[i] = new TH2D(hname,htitle, 100, -0.05, 0.05, 100, -0.05, 0.05);

		sprintf(hname,"hW_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"W of Reaction, %s",myTgt.Get_Label(i).c_str());
		hW[i] = new TH1D(hname, htitle, 250, 0, 5);
        
        sprintf(hname,"hMx_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"M_{x} of Reaction, %s",myTgt.Get_Label(i).c_str());
        hMx[i] = new TH1D(hname, htitle, 250, 0, 5);
        
		sprintf(hname,"z_fracE_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Fractional Energy, %s",myTgt.Get_Label(i).c_str());
		z_fracE[i] = new TH1D(hname, htitle, 150, 0, 1.5);
        
		sprintf(hname,"LongMom_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Longitudinal Momentum of Reconstructed Particle, %s",myTgt.Get_Label(i).c_str());
		LongMom[i] = new TH1D(hname, htitle, 500, 0, 5);
        
		sprintf(hname,"TransMom_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Transverse Momentum of Reconstructed Particle, %s",myTgt.Get_Label(i).c_str());
		TransMom[i] = new TH1D(hname, htitle, 500, 0, 5);
        
		sprintf(hname,"IM2Photons_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of Pi0 - Single Cut, %s",myTgt.Get_Label(i).c_str());
		IM2Photons[i] = new TH2D(hname, htitle, 100, 0., 1., 9, -0.5, 8.5);
        
        sprintf(hname,"OpAng_VS_IM2Photons_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. Reconstructed Mass of Pi0, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_IM2Photons[i] = new TH2D(hname, htitle, 100, 0., 1., 100, 0, 100.);
        
		sprintf(hname,"OpAng_VS_E_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. #pi^{0} Total Energy, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_E[i] = new TH2D(hname, htitle, 350, 0., 3.5, 100, 0, 100.);

        sprintf(hname,"OpAng_VS_E_MassPi0Cut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. #pi^{0} Total Energy with IM(#pi^{0}) Cut, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_E_MassPi0Cut[i] = new TH2D(hname, htitle, 350, 0., 3.5, 100, 0, 100.);

        sprintf(hname,"IM2Pions_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #pi^{+}#pi^{-} vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
        IM2Pions_VS_IMOmega[i] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

        sprintf(hname,"IM2Pions_VS_IMOmega_AllCuts_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #pi^{+}#pi^{-} vs Reconstructed Mass of #omega, all cuts, %s",myTgt.Get_Label(i).c_str());
        IM2Pions_VS_IMOmega_AllCuts[i] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"IM2Photons_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #pi^{0} vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		IM2Photons_VS_IMOmega[i] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);
        
        sprintf(hname,"W_VS_IMOmega_AllCuts_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"W vs Reconstructed Mass of #omega, All Cuts except W, %s",myTgt.Get_Label(i).c_str());
        W_VS_IMOmega_AllCuts[i] = new TH2D(hname, htitle, 250, 0., 5.0, nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"Q2_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Q2 vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		Q2_VS_IMOmega[i] = new TH2D(hname, htitle, 100, 0., 4.0, nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"Pt_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Omega Trans. Mom. vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		Pt_VS_IMOmega[i] = new TH2D(hname, htitle, 500, 0., 5., nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"Pl_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Omega Long. Mom. vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		Pl_VS_IMOmega[i] = new TH2D(hname, htitle, 500, 0., 5., nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"OpAng_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_IMOmega[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, 100, 0., 100.);
        
		sprintf(hname,"MissMom_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Missing Momentum, %s",myTgt.Get_Label(i).c_str());
		MissMom[i] = new TH1D(hname, htitle, 600, 0, 6);
        
		sprintf(hname,"MMsq_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Missing Mass Squared, %s",myTgt.Get_Label(i).c_str());
		MMsq[i] = new TH1D(hname, htitle, 700, 0, 7);

        sprintf(hname,"IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - Single Cut, %s",myTgt.Get_Label(i).c_str());
		IMOmega[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, 9, -0.5, 8.5);
        
        sprintf(hname,"IMOmega_woCut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - All Cuts, %s",myTgt.Get_Label(i).c_str());
		IMOmega_woCut[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, 9, -0.5, 8.5);

        sprintf(hname,"IMOmega_antiCut_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - Anti-Cuts, %s",myTgt.Get_Label(i).c_str());
        IMOmega_antiCut[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, 9, -0.5, 8.5);
        
        sprintf(hname,"PtSq_Omega_AllCuts_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"#omega Meson - All ID Cuts, %s",myTgt.Get_Label(i).c_str());
		PtSq_Omega_AllCuts[i] = new TH1D(hname, htitle, nPtSq_omega, PtSq_omegaLo, PtSq_omegaHi);

        sprintf(hname,"PtSq_Omega_AllCuts_IMOmegaCut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"#omega Meson - All ID Cuts & IM(#omega) Cut, %s",myTgt.Get_Label(i).c_str());
		PtSq_Omega_AllCuts_IMOmegaCut[i] = new TH1D(hname, htitle, nPtSq_omega, PtSq_omegaLo, PtSq_omegaHi);

        sprintf(hname,"PtSq_Omega_AllCuts_IMOmegaSBCut_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"#omega Meson - All ID Cuts & IM(#omega) Cut, sideband, %s",myTgt.Get_Label(i).c_str());
        PtSq_Omega_AllCuts_IMOmegaSBCut[i] = new TH1D(hname, htitle, nPtSq_omega, PtSq_omegaLo, PtSq_omegaHi);
            
        sprintf(hname,"IM2Photons_ME%i_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of Pi0, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IM2Photons_ME[i] = new TH2D(hname, htitle, 100, 0., 1., nME_Methods, -0.5, nME_Methods-0.5);
            
        sprintf(hname,"IM2Photons_OpAng_ElecPhoton_Cut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of Pi0 with e^{-} #gamma Opening Angle Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IM2Photons_OpAng_ElecPhoton_Cut_ME[i] = new TH2D(hname, htitle, 100, 0., 1., nME_Methods, -0.5, nME_Methods-0.5);
            
        sprintf(hname,"IMOmega_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_OpAng_ElecPhoton_Cut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - e^{-} #gamma Opening Angle Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_OpAng_ElecPhoton_Cut_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_MassPi0Cut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - Pi0 Mass Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_MassPi0Cut_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_ZVertCut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - Z Vertex Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_ZVertCut_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_QsqCut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - Q^{2} Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_QsqCut_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_AllCuts_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - All Cuts, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_AllCuts_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
    }
    
    for(i=0; i<MAX_SECTORS; i++){
		sprintf(hname,"elecZVertSector%i",i+1);
		sprintf(htitle,"Z Vertex of Electron in Sector %i",i+1);
        elecZVertSector[i] = new TH1D(hname, htitle, 300, -40, -10);
    }
    
    sprintf(hname,"RelativityOpAngPhotonsA");
    sprintf(htitle,"0/180 Degree Decay Photons");
	RelativityOpAngPhotonsA = new TH2D(hname,htitle, 500, 0, 5, 200, 0, 200);

    sprintf(hname,"RelativityOpAngPhotonsB");
    sprintf(htitle,"90/-90 Degree Decay Photons");
	RelativityOpAngPhotonsB = new TH2D(hname,htitle, 500, 0, 5, 180, 0, 180);
    
    sprintf(hname,"GammaPi0");
    sprintf(htitle,"Gamma of Pi0");
	GammaPi0 = new TH1D(hname,htitle, 240, 1, 25);

    sprintf(hname,"BetaPi0");
    sprintf(htitle,"Beta of Pi0");
	BetaPi0 = new TH1D(hname,htitle, 100, 0, 1);
}

//
// WriteHist - routine to write histograms to the output file
//
void WriteHist(string RootFile){
    
    int i, j;
    
    DetectedParticles myDetPart;
    ParticleList myPartList;
    EG2Target myTgt;
    ElectronID myElecID;
    
	TFile *out = new TFile(RootFile.c_str(), "recreate");
	out->cd();
  
    // create a directory for check on kinematics
    TDirectory *cdKine = out->mkdir("Kinematics");
    cdKine->cd();
    
    q2->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    q2->GetYaxis()->SetTitle("Counts");
	q2->Write();

    q2_VS_theta->GetYaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    q2_VS_theta->GetXaxis()->SetTitle("4E_{e'}sin^{2}(0.5*#theta_{e'})");
	q2_VS_theta->Write();
    
    nu_EnergyTransfer->GetXaxis()->SetTitle("\nu (GeV)");
    nu_EnergyTransfer->GetYaxis()->SetTitle("Counts");
	nu_EnergyTransfer->Write();
    
    elecZVert->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert->GetYaxis()->SetTitle("Counts");
	elecZVert->Write();

    elecZVert_VS_Phi->GetXaxis()->SetTitle("#phi (deg.)");
    elecZVert_VS_Phi->GetYaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_VS_Phi->Write();
    
    Xvert->GetXaxis()->SetTitle("X vertex (cm)");
    Xvert->GetYaxis()->SetTitle("Particle");
    Xvert->Write();

    Yvert->GetXaxis()->SetTitle("Y vertex (cm)");
    Yvert->GetYaxis()->SetTitle("Particle");
    Yvert->Write();
    
    ZVertDiff->GetXaxis()->SetTitle(" #Delta Z (cm)");
    ZVertDiff->GetYaxis()->SetTitle("Particle");
	ZVertDiff->Write();
    
    Beta_VS_Momentum->GetXaxis()->SetTitle("Momentum (GeV/c)");
    Beta_VS_Momentum->GetYaxis()->SetTitle("#beta");
	Beta_VS_Momentum->Write();

    TotalMomentum->GetXaxis()->SetTitle("Momentum (GeV/c)");
    TotalMomentum->GetYaxis()->SetTitle("Particle");
	TotalMomentum->Write();
    
    OpAng_2Photons->GetXaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
    OpAng_2Photons->GetYaxis()->SetTitle("Counts");
	OpAng_2Photons->Write();

    OpAng_elecPhoton1->GetXaxis()->SetTitle("Opening Angle between e^{-} and #gamma_{1} (deg.)");
    OpAng_elecPhoton1->GetYaxis()->SetTitle("Counts");
	OpAng_elecPhoton1->Write();

    OpAng_elecPhoton2->GetXaxis()->SetTitle("Opening Angle between e^{-} and #gamma_{2} (deg.)");
    OpAng_elecPhoton2->GetYaxis()->SetTitle("Counts");
    OpAng_elecPhoton2->Write();
 
    for(i=0; i<myPartList.Get_nPartLabel(); i++){
        Theta_VS_Phi[i]->GetXaxis()->SetTitle("#theta (deg.)");
        Theta_VS_Phi[i]->GetYaxis()->SetTitle("#phi (deg.)");
        Theta_VS_Phi[i]->Write();
    }
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        Xvert_VS_Yvert[i]->GetXaxis()->SetTitle("X vertex (cm)");
        Xvert_VS_Yvert[i]->GetYaxis()->SetTitle("Y vertex (cm)");
        Xvert_VS_Yvert[i]->Write();
    }
    
    // create a directory for check on detector info
    TDirectory *cdDetectors = out->mkdir("Detectors");
    cdDetectors->cd();
    
    CCnphe->GetXaxis()->SetTitle("Number of Photo-electrons");
    CCnphe->GetYaxis()->SetTitle("Particle");
    CCnphe->Write();

    ECu->GetXaxis()->SetTitle("EC U (cm)");
    ECu->GetYaxis()->SetTitle("Particle");
    ECu->Write();
    
    ECv->GetXaxis()->SetTitle("EC V (cm)");
    ECv->GetYaxis()->SetTitle("Particle");
    ECv->Write();
    
    ECw->GetXaxis()->SetTitle("EC W (cm)");
    ECw->GetYaxis()->SetTitle("Particle");
    ECw->Write();

    dtime_ECSC->GetXaxis()->SetTitle("#Delta t(EC-SC) (ns)");
    dtime_ECSC->GetYaxis()->SetTitle("Particle");
    dtime_ECSC->Write();
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        ECtot_VS_P[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
    	ECtot_VS_P[i]->GetYaxis()->SetTitle("EC total energy");
        ECtot_VS_P[i]->Write();

        ECtotP_VS_P[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P[i]->Write();
        
        ECin_VS_ECout[i]->GetXaxis()->SetTitle("EC inner energy");
    	ECin_VS_ECout[i]->GetYaxis()->SetTitle("EC outer energy");
        ECin_VS_ECout[i]->Write();
    }

    // create a directory for electron ID
    TDirectory *cdTgt[myTgt.Get_nIndex()];
    
    for(i=0; i<myTgt.Get_nIndex(); i++){
        cdTgt[i] = out->mkdir(myTgt.Get_Label(i).c_str());
        cdTgt[i]->cd();

        Xvert_VS_Yvert_AllCuts[i]->GetXaxis()->SetTitle("X vertex (cm)");
        Xvert_VS_Yvert_AllCuts[i]->GetYaxis()->SetTitle("Y vertex (cm)");
        Xvert_VS_Yvert_AllCuts[i]->Write();

        Xvert_VS_Yvert_Omega[i]->GetXaxis()->SetTitle("X vertex (cm)");
        Xvert_VS_Yvert_Omega[i]->GetYaxis()->SetTitle("Y vertex (cm)");
        Xvert_VS_Yvert_Omega[i]->Write();
        
        hMx[i]->GetXaxis()->SetTitle("M_{x} (GeV)");
        hMx[i]->GetYaxis()->SetTitle("Counts");
        hMx[i]->Write();
        
        hW[i]->GetXaxis()->SetTitle("W (GeV)");
        hW[i]->GetYaxis()->SetTitle("Counts");
        hW[i]->Write();

        z_fracE[i]->GetXaxis()->SetTitle("z");
        z_fracE[i]->GetYaxis()->SetTitle("Counts");
        z_fracE[i]->Write();
        
        LongMom[i]->GetXaxis()->SetTitle("#omega Longitudinal Momentum (GeV/c)");
        LongMom[i]->GetYaxis()->SetTitle("Counts");
		LongMom[i]->Write();

        TransMom[i]->GetXaxis()->SetTitle("#omega Transverse Momentum (GeV/c)");
        TransMom[i]->GetYaxis()->SetTitle("Counts");
		TransMom[i]->Write();

        IM2Photons[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons[i]->GetYaxis()->SetTitle("Cut Index");
        IM2Photons[i]->Write();

        OpAng_VS_IM2Photons[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        OpAng_VS_IM2Photons[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
		OpAng_VS_IM2Photons[i]->Write();
        
        OpAng_VS_E[i]->GetXaxis()->SetTitle("#pi^{0} Total Energy (GeV)");
        OpAng_VS_E[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
        OpAng_VS_E[i]->Write();
        
        OpAng_VS_E_MassPi0Cut[i]->GetXaxis()->SetTitle("#pi^{0} Total Energy (GeV)");
        OpAng_VS_E_MassPi0Cut[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
        OpAng_VS_E_MassPi0Cut[i]->Write();

        IM2Pions_VS_IMOmega[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} Inv. Mass (GeV/c^{2})");
        IM2Pions_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        IM2Pions_VS_IMOmega[i]->Write();

        IM2Pions_VS_IMOmega_AllCuts[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} Inv. Mass (GeV/c^{2})");
        IM2Pions_VS_IMOmega_AllCuts[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        IM2Pions_VS_IMOmega_AllCuts[i]->Write();
        
        IM2Photons_VS_IMOmega[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        IM2Photons_VS_IMOmega[i]->Write();

        W_VS_IMOmega_AllCuts[i]->GetXaxis()->SetTitle("W (GeV)");
        W_VS_IMOmega_AllCuts[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        W_VS_IMOmega_AllCuts[i]->Write();
        
        Q2_VS_IMOmega[i]->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
        Q2_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        Q2_VS_IMOmega[i]->Write();
        
        Pt_VS_IMOmega[i]->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");
        Pt_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
		Pt_VS_IMOmega[i]->Write();

        Pl_VS_IMOmega[i]->GetXaxis()->SetTitle("Longitudinal Momentum (GeV/c)");
        Pl_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
		Pl_VS_IMOmega[i]->Write();
		
        OpAng_VS_IMOmega[i]->GetXaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        OpAng_VS_IMOmega[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
        OpAng_VS_IMOmega[i]->Write();
        
        MissMom[i]->GetXaxis()->SetTitle("Missing Momentum (GeV/c)");
        MissMom[i]->GetYaxis()->SetTitle("Counts");
		MissMom[i]->Write();
        
        MMsq[i]->GetXaxis()->SetTitle("Missing Mass Squared (GeV/c)^{2}");
        MMsq[i]->GetYaxis()->SetTitle("Counts");
		MMsq[i]->Write();

        IMOmega[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega[i]->GetYaxis()->SetTitle("Cut index");
		IMOmega[i]->Write();
        
        IMOmega_woCut[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_woCut[i]->GetYaxis()->SetTitle("Cut Index");
        IMOmega_woCut[i]->Write();

        IMOmega_antiCut[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_antiCut[i]->GetYaxis()->SetTitle("Cut Index");
        IMOmega_antiCut[i]->Write();
        
        PtSq_Omega_AllCuts[i]->GetXaxis()->SetTitle("P^{2}_{T}  (GeV/c)^{2}");
        PtSq_Omega_AllCuts[i]->GetYaxis()->SetTitle("Counts");
        PtSq_Omega_AllCuts[i]->Write();

        PtSq_Omega_AllCuts_IMOmegaCut[i]->GetXaxis()->SetTitle("P^{2}_{T}  (GeV/c)^{2}");
        PtSq_Omega_AllCuts_IMOmegaCut[i]->GetYaxis()->SetTitle("Counts");
        PtSq_Omega_AllCuts_IMOmegaCut[i]->Write();

        PtSq_Omega_AllCuts_IMOmegaSBCut[i]->GetXaxis()->SetTitle("P^{2}_{T}  (GeV/c)^{2}");
        PtSq_Omega_AllCuts_IMOmegaSBCut[i]->GetYaxis()->SetTitle("Counts");
        PtSq_Omega_AllCuts_IMOmegaSBCut[i]->Write();
        
        IM2Photons_ME[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IM2Photons_ME[i]->Write();
            
        IM2Photons_OpAng_ElecPhoton_Cut_ME[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons_OpAng_ElecPhoton_Cut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IM2Photons_OpAng_ElecPhoton_Cut_ME[i]->Write();
            
        IMOmega_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_ME[i]->Write();
            
        IMOmega_OpAng_ElecPhoton_Cut_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_OpAng_ElecPhoton_Cut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_OpAng_ElecPhoton_Cut_ME[i]->Write();
            
        IMOmega_MassPi0Cut_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_MassPi0Cut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_MassPi0Cut_ME[i]->Write();
            
        IMOmega_ZVertCut_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_ZVertCut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_ZVertCut_ME[i]->Write();
            
        IMOmega_QsqCut_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_QsqCut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_QsqCut_ME[i]->Write();
            
        IMOmega_AllCuts_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_AllCuts_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_AllCuts_ME[i]->Write();
    }
    
    out->cd();
    for(i=0; i<MAX_SECTORS; i++){
        elecZVertSector[i]->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
        elecZVertSector[i]->GetYaxis()->SetTitle("Counts");
        elecZVertSector[i]->Write();
    }

    // create a directory for check on relativisitic kinematics
    TDirectory *cdRel = out->mkdir("Relativity");
    cdRel->cd();
    
    RelativityOpAngPhotonsA->Write();
	RelativityOpAngPhotonsB->Write();
	BetaPi0->Write();
	GammaPi0->Write();

    // create a directory for electron ID
    TDirectory *cdElecID = out->mkdir("ElectronID");
    cdElecID->cd();
    
    Mom_elecID->GetXaxis()->SetTitle("Momentum (GeV)");
    Mom_elecID->GetYaxis()->SetTitle("e- ID Cut");
    Mom_elecID->Write();
    
    CCnphe_elecID->GetXaxis()->SetTitle("Number of Photo-electrons");
    CCnphe_elecID->GetYaxis()->SetTitle("e- ID Cut");
    CCnphe_elecID->Write();
    
    ECu_elecID->GetXaxis()->SetTitle("EC U (cm)");
    ECu_elecID->GetYaxis()->SetTitle("e- ID Cut");
    ECu_elecID->Write();
    
    ECv_elecID->GetXaxis()->SetTitle("EC V (cm)");
    ECv_elecID->GetYaxis()->SetTitle("e- ID Cut");
    ECv_elecID->Write();
    
    ECw_elecID->GetXaxis()->SetTitle("EC W (cm)");
    ECw_elecID->GetYaxis()->SetTitle("e- ID Cut");
    ECw_elecID->Write();
    
    dtime_ECSC_elecID->GetXaxis()->SetTitle("#Delta t(EC-SC) (ns)");
    dtime_ECSC_elecID->GetYaxis()->SetTitle("e- ID Cut");
    dtime_ECSC_elecID->Write();
    
    for(i=0; i<myElecID.Get_nElecID(); i++){
        ECtot_VS_P_elecID[i]->GetXaxis()->SetTitle("Momentum (GeV)");
        ECtot_VS_P_elecID[i]->GetYaxis()->SetTitle("EC total energy");
        ECtot_VS_P_elecID[i]->Write();
        
        ECtotP_VS_P_elecID[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P_elecID[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P_elecID[i]->Write();
        
        ECin_VS_ECout_elecID[i]->GetXaxis()->SetTitle("EC inner energy");
        ECin_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        ECin_VS_ECout_elecID[i]->Write();
        
        Mom_VS_ECout_elecID[i]->GetXaxis()->SetTitle("Momentum (GeV)");
        Mom_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        Mom_VS_ECout_elecID[i]->Write();

        ECu_VS_ECout_elecID[i]->GetXaxis()->SetTitle("EC U (cm)");
        ECu_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        ECu_VS_ECout_elecID[i]->Write();

        ECv_VS_ECout_elecID[i]->GetXaxis()->SetTitle("EC V (cm)");
        ECv_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        ECv_VS_ECout_elecID[i]->Write();
        
        ECw_VS_ECout_elecID[i]->GetXaxis()->SetTitle("EC W (cm)");
        ECw_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        ECw_VS_ECout_elecID[i]->Write();
    }
    
    q2_ECoutCut->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    q2_ECoutCut->GetYaxis()->SetTitle("Counts");
    q2_ECoutCut->Write();
    
    elecZVert_ECoutCut->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_ECoutCut->GetYaxis()->SetTitle("Counts");
    elecZVert_ECoutCut->Write();
    
    Beta_VS_Momentum_ECoutCut->GetXaxis()->SetTitle("Momentum (GeV/c)");
    Beta_VS_Momentum_ECoutCut->GetYaxis()->SetTitle("#beta");
    Beta_VS_Momentum_ECoutCut->Write();
    
    Theta_VS_Phi_ECoutCut->GetXaxis()->SetTitle("#theta (deg.)");
    Theta_VS_Phi_ECoutCut->GetYaxis()->SetTitle("#phi (deg.)");
    Theta_VS_Phi_ECoutCut->Write();

    ECtot_VS_P_ECoutCut->GetXaxis()->SetTitle("Momentum (GeV)");
    ECtot_VS_P_ECoutCut->GetYaxis()->SetTitle("EC total energy");
    ECtot_VS_P_ECoutCut->Write();
    
    ECtotP_VS_P_ECoutCut->GetXaxis()->SetTitle("Momentum (GeV/c)");
    ECtotP_VS_P_ECoutCut->GetYaxis()->SetTitle("EC_{total}/Mom.");
    ECtotP_VS_P_ECoutCut->Write();

    ECtotMinusECin_ECoutCut->GetXaxis()->SetTitle("EC_{total} - EC_{in} (GeV)");
    ECtotMinusECin_ECoutCut->GetYaxis()->SetTitle("Counts");
    ECtotMinusECin_ECoutCut->Write();
    
    q2_AntiECoutCut->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    q2_AntiECoutCut->GetYaxis()->SetTitle("Counts");
    q2_AntiECoutCut->Write();
    
    elecZVert_AntiECoutCut->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_AntiECoutCut->GetYaxis()->SetTitle("Counts");
    elecZVert_AntiECoutCut->Write();
    
    Beta_VS_Momentum_AntiECoutCut->GetXaxis()->SetTitle("Momentum (GeV/c)");
    Beta_VS_Momentum_AntiECoutCut->GetYaxis()->SetTitle("#beta");
    Beta_VS_Momentum_AntiECoutCut->Write();
    
    Theta_VS_Phi_AntiECoutCut->GetXaxis()->SetTitle("#theta (deg.)");
    Theta_VS_Phi_AntiECoutCut->GetYaxis()->SetTitle("#phi (deg.)");
    Theta_VS_Phi_AntiECoutCut->Write();
    
    ECtot_VS_P_AntiECoutCut->GetXaxis()->SetTitle("Momentum (GeV)");
    ECtot_VS_P_AntiECoutCut->GetYaxis()->SetTitle("EC total energy");
    ECtot_VS_P_AntiECoutCut->Write();
    
    ECtotP_VS_P_AntiECoutCut->GetXaxis()->SetTitle("Momentum (GeV/c)");
    ECtotP_VS_P_AntiECoutCut->GetYaxis()->SetTitle("EC_{total}/Mom.");
    ECtotP_VS_P_AntiECoutCut->Write();
    
    ECtotMinusECin_AntiECoutCut->GetXaxis()->SetTitle("EC_{total} - EC_{in} (GeV)");
    ECtotMinusECin_AntiECoutCut->GetYaxis()->SetTitle("Counts");
    ECtotMinusECin_AntiECoutCut->Write();
    
    ECin_VS_ECout_ECfid->GetXaxis()->SetTitle("EC inner energy");
    ECin_VS_ECout_ECfid->GetYaxis()->SetTitle("EC outer energy");
    ECin_VS_ECout_ECfid->Write();

    ECin_VS_ECout_elecID_All->GetXaxis()->SetTitle("EC inner energy");
    ECin_VS_ECout_elecID_All->GetYaxis()->SetTitle("EC outer energy");
    ECin_VS_ECout_elecID_All->Write();

    EC_XvsY_local->GetXaxis()->SetTitle("EC X_{local} (cm)");
    EC_XvsY_local->GetYaxis()->SetTitle("EC Y_{local} (cm)");
    EC_XvsY_local->Write();
    
    out->Close();
}

void PrintTLorentzVector(TLorentzVector TLV){
    cout <<"Px "<<TLV.Px()<<"\t";
    cout <<"Py "<<TLV.Py()<<"\t";
    cout <<"Pz "<<TLV.Pz()<<"\t";
    cout <<"E " <<TLV.E() <<"\t";
    cout <<"M " <<TLV.M() <<"\t";
    cout <<"M^2 " <<TLV.M2() <<endl;
}

#ifndef __CINT__
int main (int argc, char **argv) {

    extern char *optarg;
    int c;
    extern int optind;
    
    int i;
    int dEvents = 1000; // increment of events for processing print statement
    int MaxEvents = 0; // max. number of events to process
    int TotEvents;
    
    int targMass = 2; // mass number for target
    
    bool bBatchMode = false;    // events counter is on by default
    
    string inFile;
    string outFile = "dmsProcess_omega.root";

    float timeStart = clock(); // start time
    
    for (i = 0; i < argc; ++i) cerr << argv[i] << " "; cerr << endl;
    while ((c = getopt(argc,argv, "o:M:D:T:ih")) != -1 ) {
        switch (c) {
            case 'o': outFile = optarg; break;
            case 'M': MaxEvents = atoi(optarg); break;
            case 'D': dEvents = atoi(optarg); break;
            case 'T': targMass = atoi(optarg); break;
            case 'i': bBatchMode = true; break;
            case 'h':
                PrintUsage(argv[0]);
                exit(0);
                break;
                
            default:
                cerr << "Unrecognized argument: " << optarg << endl;
                PrintUsage(argv[0]);
                exit(0);
                break;
        }
    }
    
    BookHist(); // declare histograms
    
    for (i = optind; i < argc; ++i) {
        inFile = argv[i]; // process all arguments on command line.
        if (inFile != '-') { // we have a file to process
            
            cout << "Analyzing file " << inFile << endl; // let user know which file is being processed
            
            // process the root file and return number of processed events
            TotEvents = process(inFile,MaxEvents,dEvents,targMass);
        }
    }
    cout<<TotEvents<<" events processed"<<endl; // print out stats
    
    WriteHist(outFile); // write histograms to a file
    
    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);
}
#endif
