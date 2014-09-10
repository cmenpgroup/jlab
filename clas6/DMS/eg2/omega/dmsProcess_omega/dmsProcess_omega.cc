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
    vector<double> RangeMassOmega;
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
    double Get_OpAng_ElecPhoton_lo() {return RangeOpAng_ElecPhoton[0];};
    double Get_OpAng_ElecPhoton_hi() {return RangeOpAng_ElecPhoton[1];};
    double Get_MassOmega_lo() {return RangeMassOmega[0];};
    double Get_MassOmega_hi() {return RangeMassOmega[1];};
    bool Check_Zdiff_ElecPim(double zdiff);
    bool Check_Zdiff_ElecPip(double zdiff);
    bool Check_MassPi0(double mass);
    bool Check_QSquared(double Qsq);
    bool Check_OpAng_ElecPhoton(double OpAng);
    bool Check_MassOmega(double mass);
	void Print_Cuts();
};

EG2Cuts::EG2Cuts()
{
    CutsLabel.push_back("NoCuts");
    CutsLabel.push_back("Zdiff_Electron_PiMinus");
    CutsLabel.push_back("Zdiff_Electron_PiPlus");
    CutsLabel.push_back("MassPi0");
    CutsLabel.push_back("QSquared");
    CutsLabel.push_back("OpAng_ElecPhoton");
    CutsLabel.push_back("MassOmega");
    
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

    RangeQSquared.push_back(-100000.0); // Lower limit on Q^2 (in Gev^2)
    RangeQSquared.push_back(-1.0); // Upper limit on Q^2 (in Gev^2)
    
    RangeOpAng_ElecPhoton.push_back(12.0); // Lower limit on opening angle between e- and photon (in degrees)
    RangeOpAng_ElecPhoton.push_back(180.0); // Upper limit on opening angle between e- and photon (in degrees)

    RangeMassOmega.push_back(0.7); // Lower limit on omega mass (in Gev/c^2)
    RangeMassOmega.push_back(0.875); // Upper limit on omega mass (in Gev/c^2)
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

// check the cut on omega mass
bool EG2Cuts::Check_MassOmega(double mass)
{
	bool ret = (mass >= this->Get_MassOmega_lo() && mass < this->Get_MassOmega_hi()) ? true : false;
	
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
        }else if (this->Get_CutsLabel(ii).compare("OpAng_ElecPhoton")==0) {
            cout << "[" << this->Get_OpAng_ElecPhoton_lo() << "," << this->Get_OpAng_ElecPhoton_hi() << "] (deg.)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("MassOmega")==0) {
            cout << "[" << this->Get_MassOmega_lo() << "," << this->Get_MassOmega_hi() << "] (GeV/c^2)" << endl;
        }else{
            cout << endl;
        }
    }
}

int process (string inFile, int MaxEvents, int dEvents, int targMass) {
    int i, j, k;
    
    int Sector_index;
    int Vz_index;
    bool cutPi0Mass;
    bool cutZDiff_ElectronNPion;
    bool cutZDiff_ElectronPPion;
    bool cutQSquared;
    bool cutOpAng_ElecPhoton1;
    bool cutOpAng_ElecPhoton2;
    bool cutsAll;
    
	double TwoPhotonAngle, elecPhoton1Angle, elecPhoton2Angle;
    double Qsq, nu, dubU, z_fracEnergy;
    double sinHalfTheta;
    
    EG2Target myTgt;
    EG2Cuts myCuts;
    
    myCuts.Print_Cuts();
    
    TLorentzVector BeamMinusElectron;
    TLorentzVector TwoPhoton;
	TLorentzVector Omega;

    TLorentzVector photon1_MixedEvt;
    TLorentzVector photon2_MixedEvt;
    TLorentzVector TwoPhoton_MixedEvt;
	TLorentzVector Omega_MixedEvt;

    int iMixedEvt;
    double Mass_TwoPhoton_ME[2][NUM_ENTRIES_OFFSET];
    double Mass_Omega_ME[2][NUM_ENTRIES_OFFSET];
    
    TLorentzVector beam(0., 0., BEAM_ENERGY, sqrt(BEAM_ENERGY*BEAM_ENERGY+MASS_ELECTRON*MASS_ELECTRON));
	TLorentzVector target(0., 0., 0., MASS_PROTON);

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
        cutQSquared = false;
        cutOpAng_ElecPhoton1 = false;
        cutOpAng_ElecPhoton2 = false;
        cutsAll = false;
        
        if (!(processed % dEvents)) cout << "Processed Entries: " << processed << endl;
        if (DEBUG) reader.printEvent();
        
        reader.readEntry(processed);
        
        // get the first electron lorentz vector and vertex
		TLorentzVector elec = reader.getLorentzVector(ID_ELECTRON, 0, MASS_ELECTRON);
		TVector3 elec_vert = reader.getVertex(ID_ELECTRON, 0);
        
		//TLorentzVector prot = reader.getLorentzVector(ID_PROTON, 0, MASS_PROTON);
		//TVector3 prot_vert = reader.getVertex(ID_PROTON, 0);
        
		TLorentzVector nPion = reader.getLorentzVector(ID_PION_NEG, 0,MASS_PION_CHARGED);
		TVector3 nPion_vert = reader.getVertex(ID_PION_NEG, 0);
        
		TLorentzVector pPion = reader.getLorentzVector(ID_PION_POS, 0, MASS_PION_CHARGED);
		TVector3 pPion_vert = reader.getVertex(ID_PION_POS, 0);
        
		TLorentzVector photon1 = reader.getLorentzVector(ID_PHOTON, 0, MASS_PHOTON);
		TVector3 photon1_vert = reader.getVertex(ID_PHOTON, 0);
        
		TLorentzVector photon2 = reader.getLorentzVector(ID_PHOTON, 1, MASS_PHOTON);
		TVector3 photon2_vert = reader.getVertex(ID_PHOTON, 1);
        
        // determine which target the reaction occurred in
        Vz_index = myTgt.Get_Index(elec_vert.Z());
        
        // fill the target Lorentz vector
        switch(Vz_index){
            case 1: target.SetE(MASS_DEUTERIUM); break;
            case 2: target.SetE(targMass * MASS_PROTON); break;
            default: target.SetE(MASS_PROTON); break;
        }
        
        BeamMinusElectron = beam - elec; // Lorentz Vector Difference between beam and scattered electron
		TwoPhoton = photon1 + photon2; // Two photon Lorentz vector
		Omega = pPion + nPion + TwoPhoton; // omega Lorentz vector
        
        if(NUM_ENTRIES_OFFSET*ENTRIES_OFFSET < entries){
            for(k=0; k<NUM_ENTRIES_OFFSET; k++){
                iMixedEvt = processed + (k+1)*ENTRIES_OFFSET;
                if(iMixedEvt >= entries){
                    iMixedEvt = (k+1)*ENTRIES_OFFSET;
                }
                readerMixedEvt.readEntry(iMixedEvt);
                photon1_MixedEvt = readerMixedEvt.getLorentzVector(ID_PHOTON, 0, MASS_PHOTON);
                photon2_MixedEvt = readerMixedEvt.getLorentzVector(ID_PHOTON, 1, MASS_PHOTON);
        
                Mass_TwoPhoton_ME[0][k] = (photon1_MixedEvt + photon2).M();
                Mass_Omega_ME[0][k] = (pPion + nPion + photon1_MixedEvt + photon2).M();

                Mass_TwoPhoton_ME[1][k] = (photon1_MixedEvt + photon2).M();
                Mass_Omega_ME[1][k] = (pPion + nPion + photon1 + photon2_MixedEvt).M();
            }
        }
        
        if(DEBUG) cout <<processed<<"\t"<<Omega.M()<<"\t"<<nPion.M()<<"\t"<<pPion.M()<<endl;

        double elecNPionZVertDiff = elec_vert.Z() - nPion_vert.Z(); // z vertex difference, e- and pi-
		double elecPPionZVertDiff = elec_vert.Z() - pPion_vert.Z(); // z vertex difference, e- and pi+
		double elecPhoton1ZVertDiff = elec_vert.Z() - photon1_vert.Z(); // z vertex difference, e- and photon 1
		double elecPhoton2ZVertDiff = elec_vert.Z() - photon2_vert.Z(); // z vertex difference, e- and photon 2
        
        Qsq = BeamMinusElectron.M2(); // electron Q^2
        nu = BeamMinusElectron.E(); // energy transfered to target
        dubU = sqrt(MASS_PROTON*MASS_PROTON + Qsq + 2*MASS_PROTON*nu); // reaction W
        z_fracEnergy = Omega.E()/nu; // fractional energy taken by hadron
        
        //_________________________________
		// Fill histograms
		q2->Fill(Qsq);
        sinHalfTheta = sin(0.5*elec.Theta()); // sine of one-half the electron scattering angle theta
        q2_VS_theta->Fill(4.0*elec.E()*sinHalfTheta*sinHalfTheta,Qsq);
        
        nu_EnergyTransfer->Fill(nu);
		elecZVert->Fill(elec_vert.Z());
        
        W[Vz_index]->Fill(dubU);
        z_fracE[Vz_index]->Fill(z_fracEnergy);
        
        // plots of z vertex difference between scattered electron and other decay particle
		ZVertDiff->Fill(elecNPionZVertDiff,1);
		ZVertDiff->Fill(elecPPionZVertDiff,2);
		ZVertDiff->Fill(elecPhoton1ZVertDiff,3);
		ZVertDiff->Fill(elecPhoton2ZVertDiff,4);
        
		// plots of x vs y vertices
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

        Sector_index = GetSectorByPhi(elec.Phi());
        if(Sector_index){
            elecZVertSector[Sector_index-1]->Fill(elec_vert.Z());
        }else{
            cout << "Error in finding sector. Phi = " << elec.Phi() * TMath::RadToDeg() << endl;
        }
        
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
		IM2Photons_VS_IMOmega[Vz_index]->Fill(TwoPhoton.M(), Omega.M()); // variable = 2 photon inv.
		Q2_VS_IMOmega[Vz_index]->Fill(Qsq, Omega.M()); // variable = Q^2
		Pt_VS_IMOmega[Vz_index]->Fill(Omega.Pt(), Omega.M()); // variable = omega trans. mom.
		Pl_VS_IMOmega[Vz_index]->Fill(Omega.Pz(), Omega.M()); // variable = omega long. mom.
		OpAng_VS_IMOmega[Vz_index]->Fill(Omega.M(), TwoPhotonAngle); // variable = 2 photon opening angle

        // plots of omega inv. mass after cuts
        cutPi0Mass = myCuts.Check_MassPi0(TwoPhoton.M());
        cutZDiff_ElectronNPion = myCuts.Check_Zdiff_ElecPim(elecNPionZVertDiff);
        cutZDiff_ElectronPPion = myCuts.Check_Zdiff_ElecPip(elecPPionZVertDiff);
        cutQSquared = myCuts.Check_QSquared(Qsq);
        cutOpAng_ElecPhoton1 = myCuts.Check_OpAng_ElecPhoton(elecPhoton1Angle);
        cutOpAng_ElecPhoton2 = myCuts.Check_OpAng_ElecPhoton(elecPhoton2Angle);
        cutOmegaMass = myCuts.Check_MassOmega(Omega.M());
        
        IM2Photons[Vz_index]->Fill(TwoPhoton.M()); // inv. mass of 2 photons
		IMOmega[Vz_index]->Fill(Omega.M()); // inv. mass of pi+ pi- 2 photons
		if(cutOpAng_ElecPhoton1 && cutOpAng_ElecPhoton2) {
			IM2Photons_OpAng_ElecPhoton_Cut[Vz_index]->Fill(TwoPhoton.M());
            IMOmega_OpAng_ElecPhoton_Cut[Vz_index]->Fill(Omega.M());
		}
		if(cutPi0Mass) {
            OpAng_VS_E_MassPi0Cut[Vz_index]->Fill(TwoPhoton.E(),TwoPhotonAngle);
			IMOmega_MassPi0Cut[Vz_index]->Fill(Omega.M());
		}
		if(cutZDiff_ElectronNPion && cutZDiff_ElectronPPion){
			IMOmega_ZVertCut[Vz_index]->Fill(Omega.M());
		}
		if(cutQSquared) {
			IMOmega_QsqCut[Vz_index]->Fill(Omega.M());
		}
        cutsAll = (cutZDiff_ElectronNPion && cutZDiff_ElectronPPion && cutPi0Mass && cutQSquared && cutOpAng_ElecPhoton1 && cutOpAng_ElecPhoton2);
		if(cutsAll){
			IMOmega_AllCuts[Vz_index]->Fill(Omega.M());
			PtSq_Omega_AllCuts[Vz_index]->Fill(Omega.Perp2());
            if(){
                PtSq_Omega_AllCuts_IMOmegaCut[Vz_index]->Fill(Omega.Perp2());
            }
		}
        
        for(k=0; k<NUM_ENTRIES_OFFSET; k++){
            for(j=0; j<2; j++){
                IM2Photons_ME[Vz_index][j]->Fill(Mass_TwoPhoton_ME[j][k]); // inv. mass of 2 photons
                IMOmega_ME[Vz_index][j]->Fill(Mass_Omega_ME[j][k]); // inv. mass of pi+ pi- 2 photons
                if(cutOpAng_ElecPhoton1 && cutOpAng_ElecPhoton2) {
                    IM2Photons_OpAng_ElecPhoton_Cut_ME[Vz_index][j]->Fill(Mass_TwoPhoton_ME[j][k]);
                    IMOmega_OpAng_ElecPhoton_Cut_ME[Vz_index][j]->Fill(Mass_Omega_ME[j][k]);
                }
                if(cutPi0Mass) {
                    IMOmega_MassPi0Cut_ME[Vz_index][j]->Fill(Mass_Omega_ME[j][k]);
                }
                if(cutZDiff_ElectronNPion && cutZDiff_ElectronPPion){
                    IMOmega_ZVertCut_ME[Vz_index][j]->Fill(Mass_Omega_ME[j][k]);
                }
                if(cutQSquared) {
                    IMOmega_QsqCut_ME[Vz_index][j]->Fill(Mass_Omega_ME[j][k]);
                }
                if(cutsAll){
                    IMOmega_AllCuts_ME[Vz_index][j]->Fill(Mass_Omega_ME[j][k]);
                }
            }
        }
        
        //-----------------------------------------------------
        // plots to check special relativistic kinematics
		TLorentzVector pi0(0, 0, (photon1+photon2).Pz(), TMath::Sqrt(MASS_PION_NEUTRAL*MASS_PION_NEUTRAL + (photon1+photon2).Pz() * (photon1+photon2).Pz()));
		double gamma = pi0.Gamma();
		double beta = pi0.Beta();
        
		BetaPi0->Fill(beta);
		GammaPi0->Fill(gamma);
        
		TLorentzVector photon1_pi0Rest_caseA(0, 0, 0.5 * MASS_PION_NEUTRAL, 0.5 * MASS_PION_NEUTRAL);
		TLorentzVector photon2_pi0Rest_caseA(0, 0, -0.5 * MASS_PION_NEUTRAL, 0.5 * MASS_PION_NEUTRAL);
		TLorentzVector photon1_pi0Rest_caseB(0, 0.5 * MASS_PION_NEUTRAL, 0, 0.5 * MASS_PION_NEUTRAL);
		TLorentzVector photon2_pi0Rest_caseB(0, -0.5 * MASS_PION_NEUTRAL, 0, 0.5 * MASS_PION_NEUTRAL);
        
		TLorentzVector photon1_pi0Moving_caseA(photon1_pi0Rest_caseA.Px(), photon1_pi0Rest_caseA.Py(), gamma * (photon1_pi0Rest_caseA.Pz() - beta * photon1_pi0Rest_caseA.E()), gamma * (photon1_pi0Rest_caseA.E() - beta * photon1_pi0Rest_caseA.Pz()));
        
		TLorentzVector photon2_pi0Moving_caseA(photon2_pi0Rest_caseA.Px(), photon2_pi0Rest_caseA.Py(), gamma * (photon2_pi0Rest_caseA.Pz() - beta * photon2_pi0Rest_caseA.E()), gamma * (photon2_pi0Rest_caseA.E() - beta * photon2_pi0Rest_caseA.Pz()));
        
		TLorentzVector photon1_pi0Moving_caseB(photon1_pi0Rest_caseB.Px(), photon1_pi0Rest_caseB.Py(), gamma * (photon1_pi0Rest_caseB.Pz() - beta * photon1_pi0Rest_caseB.E()), gamma * (photon1_pi0Rest_caseB.E() - beta * photon1_pi0Rest_caseB.Pz()));
        
		TLorentzVector photon2_pi0Moving_caseB(photon2_pi0Rest_caseB.Px(), photon2_pi0Rest_caseB.Py(), gamma * (photon2_pi0Rest_caseB.Pz() - beta * photon2_pi0Rest_caseB.E()), gamma * (photon2_pi0Rest_caseB.E() - beta * photon2_pi0Rest_caseB.Pz()));
        
		RelativityOpAngPhotonsA->Fill(pi0.Pz(), TMath::ACos((photon1_pi0Moving_caseA.Px() * photon2_pi0Moving_caseA.Px() + photon1_pi0Moving_caseA.Py() * photon2_pi0Moving_caseA.Py() + photon1_pi0Moving_caseA.Pz() * photon2_pi0Moving_caseA.Pz())/(photon1_pi0Moving_caseA.P() * photon2_pi0Moving_caseA.P())) * TMath::RadToDeg());
        
		RelativityOpAngPhotonsB->Fill(pi0.Pz(), TMath::ACos((photon1_pi0Moving_caseB.Px() * photon2_pi0Moving_caseB.Px() + photon1_pi0Moving_caseB.Py() * photon2_pi0Moving_caseB.Py() + photon1_pi0Moving_caseB.Pz() * photon2_pi0Moving_caseB.Pz())/(photon1_pi0Moving_caseB.P() * photon2_pi0Moving_caseB.P())) * TMath::RadToDeg());

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
    double PtSq_omegaHi = 3.0;
    
    DetectedParticles myDetPart;
    ParticleList myPartList;
    EG2Target myTgt;
    
    sprintf(hname,"q2");
    sprintf(htitle,"Q^{2}");
    q2 = new TH1D(hname,htitle, 100, -4., 0.);

    sprintf(hname,"q2_VS_theta");
    sprintf(htitle,"Q^{2} vs. 4E_{e'}sin^{2}(0.5*#theta_{e'})");
    q2_VS_theta = new TH2D(hname,htitle, 200, 0., 1.0, 200, -4., 0.);
    
    sprintf(hname,"nu_EnergyTransfer");
    sprintf(htitle,"\nu");
    nu_EnergyTransfer = new TH1D(hname,htitle, 100, 0., 5.);
    
    sprintf(hname,"elecZVert");
    sprintf(htitle,"Z Vertex of Electron");
	elecZVert = new TH1D(hname,htitle, 300, -40, -10);

    sprintf(hname,"ZVertDiff");
    sprintf(htitle,"Difference Between Z Vertices of electron and other particle");
    ZVertDiff = new TH2D(hname,htitle, 300, -10, 10,4,0.5,4.5);

    sprintf(hname,"Beta_VS_Momentum");
    sprintf(htitle,"Beta vs Momentum");
	Beta_VS_Momentum = new TH2D(hname,htitle, 500, 0, 5, 105, 0, 1.05);

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
    
    for(i=0; i<myPartList.Get_nPartLabel(); i++){
        sprintf(hname,"Theta_VS_Phi_%s",myPartList.Get_PartLabel(i).c_str());
        sprintf(htitle,"Theta vs Phi for %s",myPartList.Get_PartLabel(i).c_str());
        Theta_VS_Phi[i] = new TH2D(hname,htitle, 180, 0, 180, 360, -180, 180);
    }
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        sprintf(hname,"Xvert_VS_Yvert_%s",myDetPart.Get_DetPartLabel(i).c_str());
        sprintf(htitle,"X Vertex vs Y Vertex, %s",myDetPart.Get_DetPartLabel(i).c_str());
        Xvert_VS_Yvert[i] = new TH2D(hname,htitle, 100, -5, 5, 100, -5, 5);
    }
    
	for(i=0; i<myTgt.Get_nIndex(); i++){
		sprintf(hname,"W_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"W of Reaction, %s",myTgt.Get_Label(i).c_str());
		W[i] = new TH1D(hname, htitle, 250, 0, 5);
        
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
		sprintf(htitle,"Reconstructed Mass of Pi0, %s",myTgt.Get_Label(i).c_str());
		IM2Photons[i] = new TH1D(hname, htitle, 100, 0., 1.);

        sprintf(hname,"IM2Photons_OpAng_ElecPhoton_Cut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of Pi0 with e^{-} #gamma Opening Angle Cut, %s",myTgt.Get_Label(i).c_str());
		IM2Photons_OpAng_ElecPhoton_Cut[i] = new TH1D(hname, htitle, 100, 0., 1.);
        
		sprintf(hname,"OpAng_VS_IM2Photons_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. Reconstructed Mass of Pi0, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_IM2Photons[i] = new TH2D(hname, htitle, 100, 0., 1., 100, 0, 100.);
        
		sprintf(hname,"OpAng_VS_E_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. #pi^{0} Total Energy, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_E[i] = new TH2D(hname, htitle, 350, 0., 3.5, 100, 0, 100.);

        sprintf(hname,"OpAng_VS_E_MassPi0Cut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. #pi^{0} Total Energy with IM(#pi^{0}) Cut, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_E_MassPi0Cut[i] = new TH2D(hname, htitle, 350, 0., 3.5, 100, 0, 100.);
        
		sprintf(hname,"IM2Photons_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #pi^{0} vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		IM2Photons_VS_IMOmega[i] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"Q2_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Q2 vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		Q2_VS_IMOmega[i] = new TH2D(hname, htitle, 100, -4., 0., nIMomega, IMomegaLo, IMomegaHi);
        
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
		sprintf(htitle,"Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		IMOmega[i] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);

        sprintf(hname,"IMOmega_OpAng_ElecPhoton_Cut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - e^{-} #gamma Opening Angle Cut, %s",myTgt.Get_Label(i).c_str());
		IMOmega_OpAng_ElecPhoton_Cut[i] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);

		sprintf(hname,"IMOmega_MassPi0Cut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - Pi0 Mass Cut, %s",myTgt.Get_Label(i).c_str());
		IMOmega_MassPi0Cut[i] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"IMOmega_ZVertCut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - Z Vertex Cut, %s",myTgt.Get_Label(i).c_str());
		IMOmega_ZVertCut[i] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);

        sprintf(hname,"IMOmega_QsqCut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - Q^{2} Cut, %s",myTgt.Get_Label(i).c_str());
		IMOmega_QsqCut[i] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);
        
        sprintf(hname,"IMOmega_AllCuts_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - All Cuts, %s",myTgt.Get_Label(i).c_str());
		IMOmega_AllCuts[i] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);

        sprintf(hname,"PtSq_Omega_AllCuts_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"#omega P^{2}_{T} - All ID Cuts, %s",myTgt.Get_Label(i).c_str());
		PtSq_Omega_AllCuts[i] = new TH1D(hname, htitle, nPtSq_omega, PtSq_omegaLo, PtSq_omegaHi);

        sprintf(hname,"PtSq_Omega_AllCuts_IMOmegaCut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"#omega P^{2}_{T} - All ID Cuts & IM(#omega) Cut, %s",myTgt.Get_Label(i).c_str());
		PtSq_Omega_AllCuts_IMOmegaCut[i] = new TH1D(hname, htitle, nPtSq_omega, PtSq_omegaLo, PtSq_omegaHi);
        
        for(j=0; j<2; j++){
            
            sprintf(hname,"IM2Photons_ME%i_%s",j+1,myTgt.Get_Label(i).c_str());
            sprintf(htitle,"Reconstructed Mass of Pi0, Mixed Evt %i, %s",j+1,myTgt.Get_Label(i).c_str());
            IM2Photons_ME[i][j] = new TH1D(hname, htitle, 100, 0., 1.);
            
            sprintf(hname,"IM2Photons_OpAng_ElecPhoton_Cut_ME%i_%s",j+1,myTgt.Get_Label(i).c_str());
            sprintf(htitle,"Reconstructed Mass of Pi0 with e^{-} #gamma Opening Angle Cut, Mixed Evt %i, %s",j+1,myTgt.Get_Label(i).c_str());
            IM2Photons_OpAng_ElecPhoton_Cut_ME[i][j] = new TH1D(hname, htitle, 100, 0., 1.);
            
            sprintf(hname,"IMOmega_ME%i_%s",j+1,myTgt.Get_Label(i).c_str());
            sprintf(htitle,"Reconstructed Mass of #omega, Mixed Evt %i, %s",j+1,myTgt.Get_Label(i).c_str());
            IMOmega_ME[i][j] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);
        
            sprintf(hname,"IMOmega_OpAng_ElecPhoton_Cut_ME%i_%s",j+1,myTgt.Get_Label(i).c_str());
            sprintf(htitle,"Reconstructed Mass of #omega - e^{-} #gamma Opening Angle Cut, Mixed Evt %i, %s",j+1,myTgt.Get_Label(i).c_str());
            IMOmega_OpAng_ElecPhoton_Cut_ME[i][j] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);
        
            sprintf(hname,"IMOmega_MassPi0Cut_ME%i_%s",j+1,myTgt.Get_Label(i).c_str());
            sprintf(htitle,"Reconstructed Mass of #omega - Pi0 Mass Cut, Mixed Evt %i, %s",j+1,myTgt.Get_Label(i).c_str());
            IMOmega_MassPi0Cut_ME[i][j] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);
        
            sprintf(hname,"IMOmega_ZVertCut_ME%i_%s",j+1,myTgt.Get_Label(i).c_str());
            sprintf(htitle,"Reconstructed Mass of #omega - Z Vertex Cut, Mixed Evt %i, %s",j+1,myTgt.Get_Label(i).c_str());
            IMOmega_ZVertCut_ME[i][j] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);
        
            sprintf(hname,"IMOmega_QsqCut_ME%i_%s",j+1,myTgt.Get_Label(i).c_str());
            sprintf(htitle,"Reconstructed Mass of #omega - Q^{2} Cut, Mixed Evt %i, %s",j+1,myTgt.Get_Label(i).c_str());
            IMOmega_QsqCut_ME[i][j] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);
        
            sprintf(hname,"IMOmega_AllCuts_ME%i_%s",j+1,myTgt.Get_Label(i).c_str());
            sprintf(htitle,"Reconstructed Mass of #omega - All Cuts, Mixed Evt %i, %s",j+1,myTgt.Get_Label(i).c_str());
            IMOmega_AllCuts_ME[i][j] = new TH1D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi);
        }
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
    
	TFile *out = new TFile(RootFile.c_str(), "recreate");
	out->cd();
    
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

    for(i=0; i<myTgt.Get_nIndex(); i++){
        W[i]->GetXaxis()->SetTitle("W (GeV)");
        W[i]->GetYaxis()->SetTitle("Counts");
        W[i]->Write();

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
        IM2Photons[i]->GetYaxis()->SetTitle("Counts");
		IM2Photons[i]->Write();

        IM2Photons_OpAng_ElecPhoton_Cut[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons_OpAng_ElecPhoton_Cut[i]->GetYaxis()->SetTitle("Counts");
        IM2Photons_OpAng_ElecPhoton_Cut[i]->Write();

        OpAng_VS_IM2Photons[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        OpAng_VS_IM2Photons[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
		OpAng_VS_IM2Photons[i]->Write();
        
        OpAng_VS_E[i]->GetXaxis()->SetTitle("#pi^{0} Total Energy (GeV)");
        OpAng_VS_E[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
        OpAng_VS_E[i]->Write();
        
        OpAng_VS_E_MassPi0Cut[i]->GetXaxis()->SetTitle("#pi^{0} Total Energy (GeV)");
        OpAng_VS_E_MassPi0Cut[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
        OpAng_VS_E_MassPi0Cut[i]->Write();
		
        IM2Photons_VS_IMOmega[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        IM2Photons_VS_IMOmega[i]->Write();
		
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
        IMOmega[i]->GetYaxis()->SetTitle("Counts");
		IMOmega[i]->Write();
        
        IMOmega_OpAng_ElecPhoton_Cut[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_OpAng_ElecPhoton_Cut[i]->GetYaxis()->SetTitle("Counts");
		IMOmega_OpAng_ElecPhoton_Cut[i]->Write();
        
        IMOmega_MassPi0Cut[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_MassPi0Cut[i]->GetYaxis()->SetTitle("Counts");
		IMOmega_MassPi0Cut[i]->Write();

        IMOmega_ZVertCut[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_ZVertCut[i]->GetYaxis()->SetTitle("Counts");
		IMOmega_ZVertCut[i]->Write();

        IMOmega_QsqCut[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_QsqCut[i]->GetYaxis()->SetTitle("Counts");
		IMOmega_QsqCut[i]->Write();
        
        IMOmega_AllCuts[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_AllCuts[i]->GetYaxis()->SetTitle("Counts");
        IMOmega_AllCuts[i]->Write();

        PtSq_Omega_AllCuts[i]->GetXaxis()->SetTitle("#omega P^{2}_{T}  (GeV/c)^{2}");
        PtSq_Omega_AllCuts[i]->GetYaxis()->SetTitle("Counts");
        PtSq_Omega_AllCuts[i]->Write();

        PtSq_Omega_AllCuts_IMOmegaCut[i]->GetXaxis()->SetTitle("#omega P^{2}_{T}  (GeV/c)^{2}");
        PtSq_Omega_AllCuts_IMOmegaCut[i]->GetYaxis()->SetTitle("Counts");
        PtSq_Omega_AllCuts_IMOmegaCut[i]->Write();
        
        for(j=0; j<2; j++){
            IM2Photons_ME[i][j]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
            IM2Photons_ME[i][j]->GetYaxis()->SetTitle("Counts");
            IM2Photons_ME[i][j]->Write();
            
            IM2Photons_OpAng_ElecPhoton_Cut_ME[i][j]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
            IM2Photons_OpAng_ElecPhoton_Cut_ME[i][j]->GetYaxis()->SetTitle("Counts");
            IM2Photons_OpAng_ElecPhoton_Cut_ME[i][j]->Write();
            
            IMOmega_ME[i][j]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
            IMOmega_ME[i][j]->GetYaxis()->SetTitle("Counts");
            IMOmega_ME[i][j]->Write();
            
            IMOmega_OpAng_ElecPhoton_Cut_ME[i][j]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
            IMOmega_OpAng_ElecPhoton_Cut_ME[i][j]->GetYaxis()->SetTitle("Counts");
            IMOmega_OpAng_ElecPhoton_Cut_ME[i][j]->Write();
            
            IMOmega_MassPi0Cut_ME[i][j]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
            IMOmega_MassPi0Cut_ME[i][j]->GetYaxis()->SetTitle("Counts");
            IMOmega_MassPi0Cut_ME[i][j]->Write();
            
            IMOmega_ZVertCut_ME[i][j]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
            IMOmega_ZVertCut_ME[i][j]->GetYaxis()->SetTitle("Counts");
            IMOmega_ZVertCut_ME[i][j]->Write();
            
            IMOmega_QsqCut_ME[i][j]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
            IMOmega_QsqCut_ME[i][j]->GetYaxis()->SetTitle("Counts");
            IMOmega_QsqCut_ME[i][j]->Write();
            
            IMOmega_AllCuts_ME[i][j]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
            IMOmega_AllCuts_ME[i][j]->GetYaxis()->SetTitle("Counts");
            IMOmega_AllCuts_ME[i][j]->Write();
        }
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        elecZVertSector[i]->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
        elecZVertSector[i]->GetYaxis()->SetTitle("Counts");
        elecZVertSector[i]->Write();
    }
    
	RelativityOpAngPhotonsA->Write();
	RelativityOpAngPhotonsB->Write();
	BetaPi0->Write();
	GammaPi0->Write();

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
