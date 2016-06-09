/************************************************************************/
/*   ctProcess_omega.cc                                                 */
/*                                                                      */
/*  Created by Angelo Licastro and Andy Beiter, Canisius College        */
/*  July 2014 - Modified by M. H. Wood, Canisius College                */
/*                                                                      */
/************************************************************************/

#include "ctProcess_omega.h"

int process (string inFile, int MaxEvents, int dEvents, int targMass, bool printCuts) {
    int i, ii, j, k, kk;
    
    int Sector_index, Vz_index;
    
    bool cuts_Electron;
    bool cuts_PipPim;
    bool cuts_Photons;
    
	double TwoPhotonAngle, elecPhoton1Angle, elecPhoton2Angle, PairsAngle;
    double Qsq, nu, Mx, z_fracEnergy, W;
    double sinHalfTheta;
    
    double eventStartTime; // event start time from HEAD bank
    
    PhotonID tempPhotID;
    DetectedParticles myDetPart;
    EG2Target myTgt;
    EG2Cuts myCuts;
    OmegaMixedEvent myMixEvt;
    
    if(printCuts){ // print out the cut parameters
        myCuts.Print_Cuts();
        myMixEvt.Print_Info();
    }
    
    myMixEvt.Put_NumberOfEventsToMix(1); // add number of mixed event iterations
    myMixEvt.Put_OffsetOfEventsToMix(5); // add offset of the entry number for mixed events
    
    int NUM_MIXING_METHODS = myMixEvt.Get_nLabel(); // number of methods for mixing events
    int NUM_ENTRIES_OFFSET = myMixEvt.Get_NumberOfEventsToMix(); // retreive number of mixed event iterations
    int ENTRIES_OFFSET = myMixEvt.Get_OffsetOfEventsToMix(); // retrieve offset of the entry number for mixed events
    
    TVector3 TargetV3(0.043,-0.33,0);
    
    Vertex_Corrections myVertCorr; // create the vertex correction object
    myVertCorr.Put_Target_Vertex(TargetV3); // set the target vertex positions
    
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
    
    TVector3 elec_vert_corr;

    int iMixedEvt;
    double Mass_TwoPhoton_ME[NUM_MIXING_METHODS][NUM_ENTRIES_OFFSET];
    double Mass_Omega_ME[NUM_MIXING_METHODS][NUM_ENTRIES_OFFSET];
    int ID_ELECTRON = myHistManager.GetID_ELECTRON();
    int ID_PION_NEG = myHistManager.GetID_PION_NEG();
    int ID_PION_POS = myHistManager.GetID_PION_POS();
    int ID_PHOTON = myHistManager.GetID_PHOTON();
    double MASS_PHOTON = myHistManager.GetMASS_PHOTON();
    double MASS_DEUTERIUM = myHistManager.GetMASS_DEUTERIUM();
    double BEAM_ENERGY = myHistManager.GetBEAM_ENERGY();
    double MASS_ELECTRON = myHistManager.GetMASS_ELECTRON();
    double MASS_PROTON = myHistManager.GetMASS_PROTON();
    double MASS_PION_CHARGED = myHistManager.GetMASS_PION_CHARGED();
    double MASS_PION_NEUTRAL = myHistManager.GetMASS_PION_NEUTRAL();
    double LIGHTSPEED = myHistManager.GetLIGHTSPEED();

    TLorentzVector beam(0., 0., BEAM_ENERGY, sqrt(BEAM_ENERGY*BEAM_ENERGY+MASS_ELECTRON*MASS_ELECTRON));
	TLorentzVector target(0., 0., 0., MASS_PROTON);
	TLorentzVector nucleon(0., 0., 0., MASS_PROTON);
    
    TFile *rootfile = new TFile(inFile.c_str(),"READ");
    TTree *myTree = (TTree*)rootfile->Get("Data");
    int entries = myTree->GetEntries();
    
    KineReader myKineReader(myTree);
    PartReader elecReader(myTree,"Electron");
    PartReader pipReader(myTree,"PiPlus");
    PartReader pimReader(myTree,"PiMinus");
    PartReader photon1Reader(myTree,"Photon1");
    PartReader photon2Reader(myTree,"Photon2");
    
    cout << "Entries: " << entries << endl;

//    TEventReader readerMixedEvt;
//    readerMixedEvt.addFile(inFile.c_str());
    
    int StopProcess;
    if(MaxEvents){
        StopProcess = MaxEvents;
    }else{
        StopProcess = entries;
    }

    int processed;
    for (processed = 0; processed < StopProcess; processed = processed + 1) {
        
        myCounter.Increment("Total Events"); // total events counter
        
        myCuts.InitCuts(); // initialize omega cuts
        tempPhotID.InitCuts();
        
        if (!(processed % dEvents)) cout << "Processed Entries: " << processed << endl;

        myKineReader.ReadEntry(processed);
        elecReader.ReadEntry(processed);
        pipReader.ReadEntry(processed);
        pimReader.ReadEntry(processed);
        photon1Reader.ReadEntry(processed);
        photon2Reader.ReadEntry(processed);

//        eventStartTime = reader.getStartTime(); // evetn start time
//        myHistManager.GetStartTime()->Fill(eventStartTime);
        
        // get the first electron lorentz vector and vertex
		TLorentzVector elec = elecReader.GetLorentzVector(MASS_ELECTRON);
		TVector3 elec_vert = elecReader.GetVertex();
        
        // get the pi- lorentz vector and vertex
        TLorentzVector nPion = pimReader.GetLorentzVector(MASS_PION_CHARGED);
		TVector3 nPion_vert = pimReader.GetVertex();
        
        // get the pi+ lorentz vector and vertex
		TLorentzVector pPion = pipReader.GetLorentzVector(MASS_PION_CHARGED);
		TVector3 pPion_vert = pipReader.GetVertex();

        // get the first photon lorentz vector and vertex
		TLorentzVector photon1 = photon1Reader.GetLorentzVector(MASS_PHOTON);
		TVector3 photon1_vert = photon1Reader.GetVertex();

        // get the second photon lorentz vector and vertex
		TLorentzVector photon2 = photon2Reader.GetLorentzVector(MASS_PHOTON);
		TVector3 photon2_vert = photon2Reader.GetVertex();

/*        myMixEvt.Clear_TLorentzVectors(); // initialize all particle TLorentzVectors to zero in myMixEvt
        myMixEvt.Put_Photon1(photon1,0);
        myMixEvt.Put_Photon2(photon2,0);
        myMixEvt.Put_PiPlus(pPion,0);
        myMixEvt.Put_PiMinus(nPion,0);
        myMixEvt.Reconstruct_Pi0(0);
        myMixEvt.Reconstruct_Omega(0);
*/
        Vz_index = myTgt.Get_Index(elec_vert.Z()); // determine which target the reaction occurred in
        
        switch(Vz_index){ // fill the target Lorentz vector
            case 1: target.SetE(MASS_DEUTERIUM); break;
            case 2: target.SetE(targMass * MASS_PROTON); break;
            default: target.SetE(MASS_PROTON); break;
        }
        
        BeamMinusElectron = beam - elec; // Lorentz Vector Difference between beam and scattered electron
        TwoPion = pPion + nPion; // pion pair Lorentz vector
        TwoPhoton = photon1 + photon2; // Two photon Lorentz vector
		Omega = TwoPion + TwoPhoton; // omega Lorentz vector
/*
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
                
                for(kk=0; kk<NUM_MIXING_METHODS; kk++){
                    myMixEvt.Mix_Omega(kk); // run mixing routine for each method
                    
                    TwoPhoton_MixedEvt = myMixEvt.Get_Pi0(1); // Two photon Lorentz vector from an out-of-time event
                    Mass_TwoPhoton_ME[kk][k] = TwoPhoton_MixedEvt.M();

                    Omega_MixedEvt = myMixEvt.Get_Omega(1); // Omega Lorentz vector from an out-of-time event
                    Mass_Omega_ME[kk][k] = Omega_MixedEvt.M();
                }
            }
        }
*/
        if(DEBUG) cout <<processed<<"\t"<<Omega.M()<<"\t"<<nPion.M()<<"\t"<<pPion.M()<<endl;
        
        double elecNPionZVertDiff = elec_vert.Z() - nPion_vert.Z(); // z vertex difference, e- and pi-
		double elecPPionZVertDiff = elec_vert.Z() - pPion_vert.Z(); // z vertex difference, e- and pi+
		double elecPhoton1ZVertDiff = elec_vert.Z() - photon1_vert.Z(); // z vertex difference, e- and photon 1
		double elecPhoton2ZVertDiff = elec_vert.Z() - photon2_vert.Z(); // z vertex difference, e- and photon 2
        
        
        Qsq = myKineReader.Get_Q2(); // electron Q^2
        nu = myKineReader.Get_Nu(); // energy transfered to target
        W = myKineReader.Get_W(); // reaction W
        
        Mx_TLV = BeamMinusElectron + nucleon - Omega;
        Mx = Mx_TLV.M(); // reaction Mx
        z_fracEnergy = Omega.E()/nu; // fractional energy taken by hadron

        // correct the electron vertex
        myVertCorr.Put_Particle_Vertex(elec_vert);
        myVertCorr.Put_Particle_Dir(elec.Vect());
        myVertCorr.Put_Particle_Phi(elec.Phi());
        myVertCorr.Correct_Vertex();
        elec_vert_corr = myVertCorr.Get_Particle_Vertex_Corrected();

        // Find the electron sector
        Sector_index = GetSectorByPhi(elec.Phi());
        if(Sector_index){
            myHistManager.GetElecZVertSector()->Fill(elec_vert.Z(),Sector_index);
            myHistManager.GetElecZVertSector_Corr()->Fill(elec_vert_corr.Z(),Sector_index);
        }else{
            cout << "Error in finding sector. Phi = " << elec.Phi() * TMath::RadToDeg() << endl;
        }
        
        //_________________________________
		// Fill histograms
		myHistManager.GetQ2()->Fill(Qsq);
        sinHalfTheta = sin(0.5*elec.Theta()); // sine of one-half the electron scattering angle theta
        myHistManager.GetQ2_VS_theta()->Fill(4.0*elec.E()*sinHalfTheta*sinHalfTheta,Qsq);
        
        myHistManager.GetNu_EnergyTransfer()->Fill(nu);
		myHistManager.GetElecZVert()->Fill(elec_vert.Z()); // fill electron z vertex histogram
		myHistManager.GetElecZVert_VS_Phi()->Fill(elec.Phi() * TMath::RadToDeg(),elec_vert.Z()); // fill electron z vertex vs phi histogram
        myHistManager.GetElecZVert_VS_Phi_Corr()->Fill(elec.Phi() * TMath::RadToDeg(),elec_vert_corr.Z()); // fill electron z vertex vs phi histogram
        
        myHistManager.GetHMx(Vz_index)->Fill(Mx); // histogram for Mx
        myHistManager.GetHW(Vz_index)->Fill(W); // histogram for W
        myHistManager.GetZ_fracE(Vz_index)->Fill(z_fracEnergy); // histogram for fractional z

        myHistManager.GetNumDetPart()->Fill(myKineReader.Get_nElec(), myDetPart.Get_PartIndex("Electron"));
        myHistManager.GetNumDetPart()->Fill(myKineReader.Get_nPim(), myDetPart.Get_PartIndex("Pi-"));
        myHistManager.GetNumDetPart()->Fill(myKineReader.Get_nPip(), myDetPart.Get_PartIndex("Pi+"));
        myHistManager.GetNumDetPart()->Fill(myKineReader.Get_nGam(), myDetPart.Get_PartIndex("Photon1"));
        
        // plots of z vertex difference between scattered electron and other decay particle
		myHistManager.GetZVertDiff()->Fill(elecNPionZVertDiff,myDetPart.Get_PartIndex("Pi-"));
		myHistManager.GetZVertDiff()->Fill(elecPPionZVertDiff,myDetPart.Get_PartIndex("Pi+"));
		myHistManager.GetZVertDiff()->Fill(elecPhoton1ZVertDiff,myDetPart.Get_PartIndex("Photon1"));
		myHistManager.GetZVertDiff()->Fill(elecPhoton2ZVertDiff,myDetPart.Get_PartIndex("Photon2"));
        
        // plots of x vs y vertices
        myHistManager.GetXvert()->Fill(elec_vert.X(), myDetPart.Get_PartIndex("Electron"));
        myHistManager.GetXvert()->Fill(nPion_vert.X(), myDetPart.Get_PartIndex("Pi-"));
        myHistManager.GetXvert()->Fill(pPion_vert.X(), myDetPart.Get_PartIndex("Pi+"));
        myHistManager.GetXvert()->Fill(photon1_vert.X(), myDetPart.Get_PartIndex("Photon1"));
        myHistManager.GetXvert()->Fill(photon2_vert.X(), myDetPart.Get_PartIndex("Photon2"));

        myHistManager.GetYvert()->Fill(elec_vert.Y(), myDetPart.Get_PartIndex("Electron"));
        myHistManager.GetYvert()->Fill(nPion_vert.Y(), myDetPart.Get_PartIndex("Pi-"));
        myHistManager.GetYvert()->Fill(pPion_vert.Y(), myDetPart.Get_PartIndex("Pi+"));
        myHistManager.GetYvert()->Fill(photon1_vert.Y(), myDetPart.Get_PartIndex("Photon1"));
        myHistManager.GetYvert()->Fill(photon2_vert.Y(), myDetPart.Get_PartIndex("Photon2"));
        
        myHistManager.GetXvert_VS_Yvert(myDetPart.Get_PartIndex("Electron"))->Fill(elec_vert.X(), elec_vert.Y());
		myHistManager.GetXvert_VS_Yvert(myDetPart.Get_PartIndex("Pi-"))->Fill(nPion_vert.X(), nPion_vert.Y());
		myHistManager.GetXvert_VS_Yvert(myDetPart.Get_PartIndex("Pi+"))->Fill(pPion_vert.X(), pPion_vert.Y());
		myHistManager.GetXvert_VS_Yvert(myDetPart.Get_PartIndex("Photon1"))->Fill(photon1_vert.X(), photon1_vert.Y());
		myHistManager.GetXvert_VS_Yvert(myDetPart.Get_PartIndex("Photon2"))->Fill(photon2_vert.X(), photon2_vert.Y());
        
        // plots of angles theta vs phi
		myHistManager.GetTheta_VS_Phi(0)->Fill(elec.Theta() * TMath::RadToDeg(), elec.Phi() * TMath::RadToDeg());
		myHistManager.GetTheta_VS_Phi(1)->Fill(nPion.Theta() * TMath::RadToDeg(), nPion.Phi() * TMath::RadToDeg());
		myHistManager.GetTheta_VS_Phi(2)->Fill(pPion.Theta() * TMath::RadToDeg(), pPion.Phi() * TMath::RadToDeg());
		myHistManager.GetTheta_VS_Phi(3)->Fill(photon1.Theta() * TMath::RadToDeg(), photon1.Phi() * TMath::RadToDeg());
		myHistManager.GetTheta_VS_Phi(4)->Fill(photon2.Theta() * TMath::RadToDeg(), photon2.Phi() * TMath::RadToDeg());
		myHistManager.GetTheta_VS_Phi(5)->Fill(TwoPhoton.Theta() * TMath::RadToDeg(), TwoPhoton.Phi() * TMath::RadToDeg());
		myHistManager.GetTheta_VS_Phi(6)->Fill(Omega.Theta() * TMath::RadToDeg(), Omega.Phi() * TMath::RadToDeg());

        // plots of total momentum
		myHistManager.GetTotalMomentum()->Fill(elec.P(),0);
		myHistManager.GetTotalMomentum()->Fill(nPion.P(),1);
		myHistManager.GetTotalMomentum()->Fill(pPion.P(),2);
		myHistManager.GetTotalMomentum()->Fill(photon1.P(),3);
		myHistManager.GetTotalMomentum()->Fill(photon2.P(),4);
		myHistManager.GetTotalMomentum()->Fill(TwoPhoton.P(),5);
		myHistManager.GetTotalMomentum()->Fill(Omega.P(),6);

		// plots of beta vs momentum
		myHistManager.GetBeta_VS_Momentum()->Fill(elec.P(), elec.Beta());
		myHistManager.GetBeta_VS_Momentum()->Fill(nPion.P(), nPion.Beta());
		myHistManager.GetBeta_VS_Momentum()->Fill(pPion.P(), pPion.Beta());
		myHistManager.GetBeta_VS_Momentum()->Fill(photon1.P(), photon1.Beta());
		myHistManager.GetBeta_VS_Momentum()->Fill(photon2.P(), photon2.Beta());
        
        Fill_EC_Histograms(elecReader,myDetPart.Get_PartIndex("Electron"));
        Fill_EC_Histograms(pimReader,myDetPart.Get_PartIndex("Pi-"));
        Fill_EC_Histograms(pipReader,myDetPart.Get_PartIndex("Pi+"));
        Fill_EC_Histograms(photon1Reader,myDetPart.Get_PartIndex("Photon1"));
        Fill_EC_Histograms(photon2Reader,myDetPart.Get_PartIndex("Photon2"));

        cuts_PipPim = false; // initialize the charged pion cut
        cuts_PipPim = Analyze_PipPim(pimReader, pipReader); // charged pion pair analysis
        
        cuts_Electron = false; // initialize the electron cut
        cuts_Electron = Analyze_Electron(elecReader, Qsq, targMass); // electron analysis

        cuts_Photons = false; // initialize the photon pair cut
        cuts_Photons = Analyze_Photons(photon1Reader, photon2Reader); // photon pair analysis
        
        double emSCMassSq = elecReader.Get_TOF_MassSquared(); //calculate the TOF mass-squared for e-
        double pimSCMassSq = pimReader.Get_TOF_MassSquared(); // calculate the TOF mass-squared for pi-
        double pipSCMassSq = pipReader.Get_TOF_MassSquared(); // calculate the TOF mass-squared for pi+
        double scMassSq_phot1 = photon1Reader.Get_TOF_MassSquared(); // calculate the TOF mass-squared for photon 1
        double scMassSq_phot2 = photon2Reader.Get_TOF_MassSquared(); // calculate the TOF mass-squared for photon 2
        
        // no cuts, TOF mass
        myHistManager.GetScMassSquared_NC()->Fill(emSCMassSq,myDetPart.Get_PartIndex("Electron"));
        myHistManager.GetScMassSquared_NC()->Fill(pimSCMassSq,myDetPart.Get_PartIndex("Pi-"));
        myHistManager.GetScMassSquared_NC()->Fill(pipSCMassSq,myDetPart.Get_PartIndex("Pi+"));
        myHistManager.GetScMassSquared_NC()->Fill(scMassSq_phot1,myDetPart.Get_PartIndex("Photon1"));
        myHistManager.GetScMassSquared_NC()->Fill(scMassSq_phot2,myDetPart.Get_PartIndex("Photon2"));
        
        if(cuts_Electron){ // electron ID
            myHistManager.GetEChit_M2_cuts()->Fill(elecReader.Get_EChit_M2(),myDetPart.Get_PartIndex("Electron"));
            myHistManager.GetEChit_M3_cuts()->Fill(elecReader.Get_EChit_M3(),myDetPart.Get_PartIndex("Electron"));
            myHistManager.GetEChit_M4_cuts()->Fill(elecReader.Get_EChit_M4(),myDetPart.Get_PartIndex("Electron"));
            
            myHistManager.GetScMassSquared_EC()->Fill(emSCMassSq,myDetPart.Get_PartIndex("Electron"));
            myHistManager.GetScMassSquared_EC()->Fill(pimSCMassSq,myDetPart.Get_PartIndex("Pi-"));
            myHistManager.GetScMassSquared_EC()->Fill(pipSCMassSq,myDetPart.Get_PartIndex("Pi+"));
            myHistManager.GetScMassSquared_EC()->Fill(scMassSq_phot1,myDetPart.Get_PartIndex("Photon1"));
            myHistManager.GetScMassSquared_EC()->Fill(scMassSq_phot2,myDetPart.Get_PartIndex("Photon2"));
        }

        if(cuts_Photons){ // photon ID
            myHistManager.GetEChit_M2_cuts()->Fill(photon1Reader.Get_EChit_M2(),myDetPart.Get_PartIndex("Photon1"));
            myHistManager.GetEChit_M3_cuts()->Fill(photon1Reader.Get_EChit_M3(),myDetPart.Get_PartIndex("Photon1"));
            myHistManager.GetEChit_M4_cuts()->Fill(photon1Reader.Get_EChit_M4(),myDetPart.Get_PartIndex("Photon1"));

            myHistManager.GetEChit_M2_cuts()->Fill(photon2Reader.Get_EChit_M2(),myDetPart.Get_PartIndex("Photon2"));
            myHistManager.GetEChit_M3_cuts()->Fill(photon2Reader.Get_EChit_M3(),myDetPart.Get_PartIndex("Photon2"));
            myHistManager.GetEChit_M4_cuts()->Fill(photon2Reader.Get_EChit_M4(),myDetPart.Get_PartIndex("Photon2"));
            
            myHistManager.GetScMassSquared_PC()->Fill(emSCMassSq,myDetPart.Get_PartIndex("Electron"));
            myHistManager.GetScMassSquared_PC()->Fill(pimSCMassSq,myDetPart.Get_PartIndex("Pi-"));
            myHistManager.GetScMassSquared_PC()->Fill(pipSCMassSq,myDetPart.Get_PartIndex("Pi+"));
            myHistManager.GetScMassSquared_PC()->Fill(scMassSq_phot1,myDetPart.Get_PartIndex("Photon1"));
            myHistManager.GetScMassSquared_PC()->Fill(scMassSq_phot2,myDetPart.Get_PartIndex("Photon2"));
        }
        
        if(cuts_PipPim){ // charged pion pair ID
            myHistManager.GetEChit_M2_cuts()->Fill(pimReader.Get_EChit_M2(),myDetPart.Get_PartIndex("Pi-"));
            myHistManager.GetEChit_M3_cuts()->Fill(pimReader.Get_EChit_M3(),myDetPart.Get_PartIndex("Pi-"));
            myHistManager.GetEChit_M4_cuts()->Fill(pimReader.Get_EChit_M4(),myDetPart.Get_PartIndex("Pi-"));

            myHistManager.GetEChit_M2_cuts()->Fill(pipReader.Get_EChit_M2(),myDetPart.Get_PartIndex("Pi+"));
            myHistManager.GetEChit_M3_cuts()->Fill(pipReader.Get_EChit_M3(),myDetPart.Get_PartIndex("Pi+"));
            myHistManager.GetEChit_M4_cuts()->Fill(pipReader.Get_EChit_M4(),myDetPart.Get_PartIndex("Pi+"));
        }
        
        // 1 - pi+ pi-
        TLorentzVector pipPi0 = pPion + TwoPhoton; // 2 - pi+ pi0
        TLorentzVector pimPi0 = nPion + TwoPhoton; // 3 - pi- pi0

        // NC - no cuts
        myHistManager.GetMass2Pions_VS_massOmega_NC(0)->Fill(TwoPion.M(), Omega.M()); // 1st pion pair
        myHistManager.GetMass2Pions_VS_massOmega_NC(1)->Fill(pipPi0.M(), Omega.M()); // 2nd pion pair
        myHistManager.GetMass2Pions_VS_massOmega_NC(2)->Fill(pimPi0.M(), Omega.M()); // 3rd pion pair

        //
        // Start omega ID
        //
        if(cuts_Electron && cuts_Photons && cuts_PipPim){
            // EPC - electron, photon, pion cuts
            myHistManager.GetDBeta_VS_Momentum_EPC(0)->Fill(elec.P(), elecReader.Get_BetaDifference(MASS_ELECTRON));
            myHistManager.GetDBeta_VS_Momentum_EPC(1)->Fill(nPion.P(), pimReader.Get_BetaDifference(MASS_PION_CHARGED));
            myHistManager.GetDBeta_VS_Momentum_EPC(2)->Fill(pPion.P(), pipReader.Get_BetaDifference(MASS_PION_CHARGED));
            myHistManager.GetDBeta_VS_Momentum_EPC(3)->Fill(photon1.P(), photon1Reader.Get_BetaDifference(MASS_PHOTON));
            myHistManager.GetDBeta_VS_Momentum_EPC(4)->Fill(photon2.P(), photon2Reader.Get_BetaDifference(MASS_PHOTON));

            // electron, photon, pion cuts, TOF mass
            myHistManager.GetScMassSquared_EPC()->Fill(emSCMassSq,myDetPart.Get_PartIndex("Electron"));
            myHistManager.GetScMassSquared_EPC()->Fill(pimSCMassSq,myDetPart.Get_PartIndex("Pi-"));
            myHistManager.GetScMassSquared_EPC()->Fill(pipSCMassSq,myDetPart.Get_PartIndex("Pi+"));
            myHistManager.GetScMassSquared_EPC()->Fill(scMassSq_phot1,myDetPart.Get_PartIndex("Photon1"));
            myHistManager.GetScMassSquared_EPC()->Fill(scMassSq_phot2,myDetPart.Get_PartIndex("Photon2"));
            
            // plot of two photon opening angle
            TwoPhotonAngle = TMath::RadToDeg()*photon1.Angle(photon2.Vect());
            myHistManager.GetOpAng_2Photons()->Fill(TwoPhotonAngle);

            elecPhoton1Angle = TMath::RadToDeg()*elec.Angle(photon1.Vect());
            myHistManager.GetOpAng_elecPhoton1()->Fill(elecPhoton1Angle);

            elecPhoton2Angle = TMath::RadToDeg()*elec.Angle(photon2.Vect());
            myHistManager.GetOpAng_elecPhoton2()->Fill(elecPhoton2Angle);
            
            // plots by target (Vz_index)
            myHistManager.GetLongMom(Vz_index)->Fill(Omega.P()-Omega.Pt()); // omega long. mom.
            myHistManager.GetTransMom(Vz_index)->Fill(Omega.Pt()); // omega trans. mom.
            
            myHistManager.GetOpAng_VS_IM2Photons(Vz_index)->Fill(TwoPhoton.M(),TwoPhotonAngle); // opening angle vs 2 photon inv. mass
            myHistManager.GetOpAng_VS_E(Vz_index)->Fill(TwoPhoton.E(),TwoPhotonAngle); // opening angle vs 2 photon total energy
            
            myHistManager.GetMissMom(Vz_index)->Fill((beam + target - Omega).P());  // mising mom.
            myHistManager.GetMMsq(Vz_index)->Fill((beam + target - Omega).M2()); // missing mass^2
            
            // plots of variable vs the omega inv. mass
            myHistManager.GetIM2Pions_VS_IMOmega(Vz_index)->Fill(TwoPion.M(), Omega.M()); // variable = pion pair inv. mass
            myHistManager.GetIM2Photons_VS_IMOmega(Vz_index)->Fill(TwoPhoton.M(), Omega.M()); // variable = 2 photon inv. mass
            myHistManager.GetQ2_VS_IMOmega(Vz_index)->Fill(Qsq, Omega.M()); // variable = Q^2
            myHistManager.GetPt_VS_IMOmega(Vz_index)->Fill(Omega.Pt(), Omega.M()); // variable = omega trans. mom.
            myHistManager.GetPl_VS_IMOmega(Vz_index)->Fill(Omega.Pz(), Omega.M()); // variable = omega long. mom.
            myHistManager.GetOpAng_VS_IMOmega(Vz_index)->Fill(TwoPhotonAngle, Omega.M()); // variable = 2 photon opening angle

            // set the cuts
            myCuts.SetCut_MassPi0(TwoPhoton.M()); // pi0 mass cut
            myCuts.SetCut_QSquared(Qsq); // Q^2 cut
            myCuts.SetCut_Wcut(W); // W cut
            myCuts.SetCut_MassPipPim(TwoPion.M()); // pi+ pi- inv. mass cut
            myCuts.SetCut_NumDetPart(myKineReader.Get_nElec(),myKineReader.Get_nPim(),myKineReader.Get_nPip(),myKineReader.Get_nGam()); // topology cut
            
            myCuts.SetCut_ZDiff_ElecPion(elecNPionZVertDiff,0); // e-,pi- z-vertex matching cut
            myCuts.SetCut_ZDiff_ElecPion(elecPPionZVertDiff,1); // e-,pi+ z-vertex matching cut
            myCuts.SetCut_ZDiff_ElecPion_All(); // final e-,pion z-vertex matching cut

            myCuts.SetCut_OpAng_ElecPhoton(elecPhoton1Angle,1); // e-,photon 1 opening angle cut
            myCuts.SetCut_OpAng_ElecPhoton(elecPhoton2Angle,2); // e-,photon 2 opening angle cut
            myCuts.SetCut_OpAng_ElecPhoton_All(); // final e-,photon opening angle cut
            
            myCuts.SetCut_MassOmega(Omega.M()); // omega mass cut
            myCuts.SetCut_MassOmega_sb(Omega.M()); // sideband cuts on the omega mass
            
            // omega inv. mass histograms
            myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),0); // inv. mass of 2 photons
            myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons
            myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons
            myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),0); // inv. mass of 2 photons
            myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons

            if(myCuts.GetCut_MassPi0()) { // applying the pi0 mass cut
                myCounter.Increment("Omega ID (Mpi0)");
                myHistManager.GetOpAng_VS_E_MassPi0Cut(Vz_index)->Fill(TwoPhoton.E(),TwoPhotonAngle);
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),2);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),2);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),2);
            }
            myCuts.SetCut_OmegaID_woMassPi0();
            if(myCuts.GetCut_OmegaID_woMassPi0()){ // applying all cuts except the pi0 mass
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),2);
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),2);
            }
            
            if(myCuts.GetCut_QSquared()) { // applying the Q^2 cut
                myCounter.Increment("Omega ID (Q2)");
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),3);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),3);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),3);
            }
            myCuts.SetCut_OmegaID_woQSquared();
            if(myCuts.GetCut_OmegaID_woQSquared()){ // applying all cuts except Q^2
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),3);
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),3);
            }
            
            if(myCuts.GetCut_Wcut()){ // applying the W cut
                myCounter.Increment("Omega ID (W)");
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),4);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),4);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),4);
            }
            myCuts.SetCut_OmegaID_woWcut();
            if(myCuts.GetCut_OmegaID_woWcut()){ // applying all cuts except W
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),4);
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),4);
            }
            
            if(myCuts.GetCut_ZDiff_ElecPion_All()){ // applying the e-,pion z-vertex matching cut
                myCounter.Increment("Omega ID (Zmatch)");
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),5);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),5);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),5);
            }
            myCuts.SetCut_OmegaID_woZDiff();
            if(myCuts.GetCut_OmegaID_woZDiff()){ // applying all cuts except the e-,pion z-vertex matching
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),5);
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),5);
            }
            
            if(myCuts.GetCut_OpAng_ElecPhoton_All()) { // applying the e-,photon opening angle cut
                myCounter.Increment("Omega ID (OpAngElectron)");
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),6);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),6);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),6);
            }
            myCuts.SetCut_OmegaID_woOpAng_ElecPhoton();
            if(myCuts.GetCut_OmegaID_woOpAng_ElecPhoton()){ // applying all cuts except the e-,photon opening angle
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),6);
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),6);
            }
            
            if(myCuts.GetCut_MassPipPim()) { // applying the pi+pi- mass cut
                myCounter.Increment("Omega ID (MPipPim)");
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),7);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),7);
                myHistManager.GetMass2Pions_VS_massOmega_EPC(0)->Fill(TwoPion.M(), Omega.M());
                myHistManager.GetMass2Pions_VS_massOmega_EPC(1)->Fill(pipPi0.M(), Omega.M());
                myHistManager.GetMass2Pions_VS_massOmega_EPC(2)->Fill(pimPi0.M(), Omega.M());
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),7);
            }
            myCuts.SetCut_OmegaID_woMassPipPim();
            if(myCuts.GetCut_OmegaID_woMassPipPim()){ // applying all cuts except the e-,photon opening angle
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),7);
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),7);
            }

            if(myCuts.GetCut_NumDetPart()) { // applying the particle topology cut
                myCounter.Increment("Omega ID (NumDetPart)");
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),8);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),8);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),8);
            }
            myCuts.SetCut_OmegaID_woNumDetPart();
            if(myCuts.GetCut_OmegaID_woNumDetPart()){ // applying all cuts except the e-,photon opening angle
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),8);
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),8);
            }
            
            if((photon1Reader.Get_EChit_M2() < 20) || (photon2Reader.Get_EChit_M2() < 20)){
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),9);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),9);
            }else{
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),9);
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),9);
            }

            if((photon1Reader.Get_EChit_M2()>=20 && photon1Reader.Get_EChit_M2()<70) && (photon2Reader.Get_EChit_M2()>=20 && photon2Reader.Get_EChit_M2()<70)){
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),10);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),10);
            }else{
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),10);
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),10);
            }
            
            if((photon1Reader.Get_EChit_M2() >= 70) && (photon2Reader.Get_EChit_M2() >= 70)){
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),11);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),11);
            }else{
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),11);
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),11);
            }

            if((photon1Reader.Get_EChit_M3() < 20) || (photon2Reader.Get_EChit_M3() < 20)){
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),12);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),12);
            }else{
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),12);
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),12);
            }
            
            if((photon1Reader.Get_EChit_M3()>=20 && photon1Reader.Get_EChit_M3()<200) && (photon2Reader.Get_EChit_M3()>=20 && photon2Reader.Get_EChit_M3()<200)){
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),13);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),13);
            }else{
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),13);
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),13);
            }
            
            if((photon1Reader.Get_EChit_M3() >= 200) && (photon2Reader.Get_EChit_M3() >= 200)){
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),14);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),14);
            }else{
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),14);
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),14);
            }
            
            tempPhotID.SetCut_PhotonSCMsq(scMassSq_phot1,1);
            tempPhotID.SetCut_PhotonSCMsq(scMassSq_phot2,2);
            tempPhotID.SetCut_PhotonSCMsq_All();
            if(tempPhotID.GetCut_PhotonSCMsq_All()){
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),15);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),15);
            }else{
                myHistManager.GetIM2Photons_woCut(Vz_index)->Fill(TwoPhoton.M(),15);
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),15);
            }
            
             // applying all cuts
            myCuts.SetCut_OmegaID();
            if(myCuts.GetCut_OmegaID()){
                myCounter.Increment("Omega ID (All)");
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),1);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),1);
                myHistManager.GetW_VS_IMOmega_AllCuts(Vz_index)->Fill(W, Omega.M()); // variable = W
                myHistManager.GetIM2Pions_VS_IMOmega_AllCuts(Vz_index)->Fill(TwoPion.M(), Omega.M()); // variable = pion pair inv. mass
                myHistManager.GetPtSq_Omega_AllCuts(Vz_index)->Fill(Omega.Perp2());

                myHistManager.GetXvert_VS_Yvert_AllCuts(Vz_index)->Fill(elec_vert.X(), elec_vert.Y());

                myHistManager.GetVirtualPhotonAngle_VS_IMOmega_AllCuts(Vz_index)->Fill(BeamMinusElectron.Theta()*TMath::RadToDeg(),Omega.M());
                myHistManager.GetOpAngVPomega_VS_IMOmega_AllCuts(Vz_index)->Fill(Omega.Angle(BeamMinusElectron.Vect())*TMath::RadToDeg(),Omega.M());
                myHistManager.GetPt_VS_IMOmega_AllCuts(Vz_index)->Fill(Omega.Pt(),Omega.M());
                myHistManager.GetPl_VS_Pt_AllCuts(Vz_index)->Fill(Omega.Pz(),Omega.Pt());
                
                PairsAngle = TMath::RadToDeg()*TwoPhoton.Angle(TwoPion.Vect());
                
                myHistManager.GetOpAng_VS_IMOmega_AllCuts(Vz_index)->Fill(TwoPhotonAngle, Omega.M()); // variable = 2 photon opening angle
                myHistManager.GetOpAngPairs_VS_IMOmega_AllCuts(Vz_index)->Fill(PairsAngle, Omega.M()); // variable = angular difference between 2 photons and 2 pions
                
                if(myCuts.GetCut_MassPi0()) {
                    // EPC - electron, photon, pion, omega cuts
                    myHistManager.GetMass2Pions_VS_massOmega_EPOC(0)->Fill(TwoPion.M(), Omega.M()); //
                    myHistManager.GetMass2Pions_VS_massOmega_EPOC(1)->Fill(pipPi0.M(), Omega.M()); //
                    myHistManager.GetMass2Pions_VS_massOmega_EPOC(2)->Fill(pimPi0.M(), Omega.M()); //
                }
                if(myCuts.GetCut_MassOmega()){
                    myHistManager.GetXvert_VS_Yvert_Omega(Vz_index)->Fill(elec_vert.X(), elec_vert.Y());
                    myHistManager.GetPtSq_Omega_AllCuts_IMOmegaCut(Vz_index)->Fill(Omega.Perp2());
                }
                if(myCuts.GetCut_MassOmega_sb()){
                    myHistManager.GetPtSq_Omega_AllCuts_IMOmegaSBCut(Vz_index)->Fill(Omega.Perp2());
                }
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),1);
            }
/*
            for(k=0; k<NUM_ENTRIES_OFFSET; k++){ // loop over number of mixed event iterations
                for(j=0; j<NUM_MIXING_METHODS; j++){ // loop over number of mixed event methods
                    myHistManager.GetIM2Photons_ME(Vz_index)->Fill(Mass_TwoPhoton_ME[j][k],j); // inv. mass of 2 photons
                    myHistManager.GetIMOmega_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j); // inv. mass of pi+ pi- 2 photons
                    if(myCuts.GetCut_OpAng_ElecPhoton_All()) {
                        myHistManager.GetIM2Photons_OpAng_ElecPhoton_Cut_ME(Vz_index)->Fill(Mass_TwoPhoton_ME[j][k],j);
                        myHistManager.GetIMOmega_OpAng_ElecPhoton_Cut_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                    if(myCuts.GetCut_MassPi0()) {
                        myHistManager.GetIMOmega_MassPi0Cut_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                    if(myCuts.GetCut_ZDiff_ElecPion_All()){
                        myHistManager.GetIMOmega_ZVertCut_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                    if(myCuts.GetCut_QSquared()) {
                        myHistManager.GetIMOmega_QsqCut_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                    if(cutsAll){
                        myHistManager.GetIMOmega_AllCuts_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                }
            }
 */
        }
        //
        // End omega ID
        //
        
        //-----------------------------------------------------
        // plots to check special relativistic kinematics
		TLorentzVector pi0(0, 0, TwoPhoton.Pz(), TMath::Sqrt(MASS_PION_NEUTRAL*MASS_PION_NEUTRAL + TwoPhoton.Pz() * TwoPhoton.Pz()));
		myHistManager.GetBetaPi0()->Fill(pi0.Beta());
		myHistManager.GetGammaPi0()->Fill(pi0.Gamma());
        
        TLorentzRotation lbr(pi0.BoostVector());
        
		TLorentzVector photon1_pi0Rest_caseA(0, 0, 0.5 * MASS_PION_NEUTRAL, 0.5 * MASS_PION_NEUTRAL);
        TLorentzVector photon2_pi0Rest_caseA(0, 0, -0.5 * MASS_PION_NEUTRAL, 0.5 * MASS_PION_NEUTRAL);
        photon1_pi0Rest_caseA.Transform(lbr);
        photon2_pi0Rest_caseA.Transform(lbr);
        
        TLorentzVector photon1_pi0Rest_caseB(0, 0.5 * MASS_PION_NEUTRAL, 0, 0.5 * MASS_PION_NEUTRAL);
        TLorentzVector photon2_pi0Rest_caseB(0, -0.5 * MASS_PION_NEUTRAL, 0, 0.5 * MASS_PION_NEUTRAL);
        photon1_pi0Rest_caseB.Transform(lbr);
        photon2_pi0Rest_caseB.Transform(lbr);
        
        myHistManager.GetRelativityOpAngPhotonsA()->Fill(pi0.Pz(),photon1_pi0Rest_caseA.Angle(photon2_pi0Rest_caseA.Vect())*TMath::RadToDeg());
        myHistManager.GetRelativityOpAngPhotonsB()->Fill(pi0.Pz(),photon1_pi0Rest_caseB.Angle(photon2_pi0Rest_caseB.Vect())*TMath::RadToDeg());
        //-----------------------------------------------------
    }
    
    rootfile->Close(); // close the root TFile
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

void Fill_EC_Histograms(PartReader Rdr, int ip)
{
    myHistManager.GetCCnphe()->Fill(Rdr.Get_CCnphe(),ip);
    myHistManager.GetECu()->Fill(Rdr.Get_ECu(),ip);
    myHistManager.GetECv()->Fill(Rdr.Get_ECv(),ip);
    myHistManager.GetECw()->Fill(Rdr.Get_ECw(),ip);
    myHistManager.GetEChit_M2()->Fill(Rdr.Get_EChit_M2(),ip);
    myHistManager.GetEChit_M3()->Fill(Rdr.Get_EChit_M3(),ip);
    myHistManager.GetEChit_M4()->Fill(Rdr.Get_EChit_M4(),ip);
    myHistManager.GetECtot_VS_P(ip)->Fill(Rdr.Get_Mom(),Rdr.Get_ECtot());
    myHistManager.GetECin_VS_ECout(ip)->Fill(Rdr.Get_ECin(),Rdr.Get_ECout());
    
    if(Rdr.Get_Mom()!=0.0){
        myHistManager.GetECtotP_VS_P(ip)->Fill(Rdr.Get_Mom(),Rdr.Get_ECtot()/Rdr.Get_Mom());
    }
    
    myHistManager.GetDtime_ECSC()->Fill(Rdr.Get_ECtime() - Rdr.Get_SCtime() - ECSC_TimeOffset,ip);
}

// Charged pion ID
bool Analyze_PipPim(PartReader pimReader, PartReader pipReader)
{
    ChargedPionID myChPionID; // create the charged pion pair object, which will initialize the cuts to false
    
    double MASS_PION_CHARGED = myHistManager.GetMASS_PION_CHARGED();
    
    double pimBeta = pimReader.Get_Beta(); // get pi- beta
    double pimMom = pimReader.Get_Mom(); // get pi- total momentum
    double pimdBeta = pimReader.Get_BetaDifference(MASS_PION_CHARGED); // calculate the beta difference
    double pimSCMassSq = pimReader.Get_TOF_MassSquared(); // calculate the TOF mass-squared for pi-
    double pimEChit_M2 = pimReader.Get_EChit_M2();
    double pimEChit_M3 = pimReader.Get_EChit_M3();
    double pimEChit_M4 = pimReader.Get_EChit_M4();
    
    myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(pimMom, pimBeta);
    myHistManager.GetBeta_Recalc()->Fill(pimBeta,1);
    myHistManager.GetDBeta_VS_Momentum(1)->Fill(pimMom, pimdBeta);
    
    myChPionID.SetCut_NegPionDiffBeta(pimdBeta); // cut on pi- beta difference
    if(myChPionID.GetCut_NegPionDiffBeta()) myCounter.Increment("Neg. Pion ID (dBeta)");
    
    double pipBeta = pipReader.Get_Beta(); // get pi+ beta
    double pipMom = pipReader.Get_Mom(); // get pi+ total momentum
    double pipdBeta = pipReader.Get_BetaDifference(MASS_PION_CHARGED); // calculate the beta difference
    double pipSCMassSq = pipReader.Get_TOF_MassSquared(); // calculate the TOF mass-squared for pi+
    double pipEChit_M2 = pipReader.Get_EChit_M2();
    double pipEChit_M3 = pipReader.Get_EChit_M3();
    double pipEChit_M4 = pipReader.Get_EChit_M4();
    
    myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(pipMom, pipBeta);
    myHistManager.GetBeta_Recalc()->Fill(pipBeta,2);
    myHistManager.GetDBeta_VS_Momentum(2)->Fill(pipMom, pipdBeta);
    
    myChPionID.SetCut_PosPionDiffBeta(pipdBeta); // cut on pi+ beta difference
    if(myChPionID.GetCut_PosPionDiffBeta()) myCounter.Increment("Pos. Pion ID (dBeta)");
    
    myChPionID.SetCut_ChargedPionPair();
    if(myChPionID.GetCut_ChargedPionPair()){
        myCounter.Increment("Charged Pion ID (All)");

        myHistManager.GetEChit_M2_VS_scMsq(1)->Fill(pimEChit_M2,pimSCMassSq);
        myHistManager.GetEChit_M3_VS_scMsq(1)->Fill(pimEChit_M3,pimSCMassSq);
        myHistManager.GetEChit_M4_VS_scMsq(1)->Fill(pimEChit_M4,pimSCMassSq);

        myHistManager.GetEChit_M2_VS_scMsq(2)->Fill(pipEChit_M2,pipSCMassSq);
        myHistManager.GetEChit_M3_VS_scMsq(2)->Fill(pipEChit_M3,pipSCMassSq);
        myHistManager.GetEChit_M4_VS_scMsq(2)->Fill(pipEChit_M4,pipSCMassSq);
    }
    
    return myChPionID.GetCut_ChargedPionPair();
}

// Photon pair ID
bool Analyze_Photons(PartReader photon1Reader, PartReader photon2Reader)
{
    PhotonID myPhotID;
    EC_geometry myECgeom;
    
    double eventStartTime = 0.0; // event start time from HEAD bank
    
    double LIGHTSPEED = myHistManager.GetLIGHTSPEED();
    double MASS_PHOTON = myHistManager.GetMASS_PHOTON();
    
    double mom_phot1 = photon1Reader.Get_Mom();
    double mom_phot2 = photon2Reader.Get_Mom();
    
    double beta_phot1 = photon1Reader.Get_Beta();
    double beta_phot2 = photon2Reader.Get_Beta(); // use beta from CLAS

    double dBeta_phot1 = photon1Reader.Get_BetaDifference(MASS_PHOTON); // difference in beta for measured and ideal beta
    double dBeta_phot2 = photon2Reader.Get_BetaDifference(MASS_PHOTON); // difference in beta for measured and ideal beta

    double scMassSq_phot1 = photon1Reader.Get_TOF_MassSquared(); // calculate the TOF mass-squared for photon 1
    double scMassSq_phot2 = photon2Reader.Get_TOF_MassSquared(); // calculate the TOF mass-squared for photon 2

    double echit_M2_phot1 = photon1Reader.Get_EChit_M2();
    double echit_M3_phot1 = photon1Reader.Get_EChit_M3();
    double echit_M4_phot1 = photon1Reader.Get_EChit_M4();

    double echit_M2_phot2 = photon2Reader.Get_EChit_M2();
    double echit_M3_phot2 = photon2Reader.Get_EChit_M3();
    double echit_M4_phot2 = photon2Reader.Get_EChit_M4();
    
    double ectime_phot1 = photon1Reader.Get_ECtime();
    double ectime_phot2 = photon2Reader.Get_ECtime();
    double ecpath_phot1 = photon1Reader.Get_ECpath();
    double ecpath_phot2 = photon2Reader.Get_ECpath();
    
    double timing_phot1 = ectime_phot1 - ecpath_phot1/LIGHTSPEED;
    double timing_phot2 = ectime_phot2 - ecpath_phot2/LIGHTSPEED;

    double ecu_phot1 = photon1Reader.Get_ECu();
    double ecu_phot2 = photon2Reader.Get_ECu();
    double ecv_phot1 = photon1Reader.Get_ECv();
    double ecv_phot2 = photon2Reader.Get_ECv();
    double ecw_phot1 = photon1Reader.Get_ECw();
    double ecw_phot2 = photon2Reader.Get_ECw();
    
    double ecin_phot1 = photon1Reader.Get_ECin();
    double ecin_phot2 = photon2Reader.Get_ECin();
    double ecout_phot1 = photon1Reader.Get_ECout();
    double ecout_phot2 = photon2Reader.Get_ECout();
    double ectot_phot1 = photon1Reader.Get_ECtot();
    double ectot_phot2 = photon2Reader.Get_ECtot();

    TLorentzVector photon1 = photon1Reader.GetLorentzVector(MASS_PHOTON);
    int Sector_index_phot1 = GetSectorByPhi(photon1.Phi());

    TLorentzVector photon2 = photon1Reader.GetLorentzVector(MASS_PHOTON);
    int Sector_index_phot2 = GetSectorByPhi(photon2.Phi());
    
    myPhotID.SetCut_PhotonMom(mom_phot1,1);
    myPhotID.SetCut_PhotonMom(mom_phot2,2);
    myPhotID.SetCut_PhotonMom_All();
    if(myPhotID.GetCut_PhotonMom_All()) myCounter.Increment("Photon ID (Mom)");

    myHistManager.GetMomentumPhoton1()->Fill(mom_phot1);
    myHistManager.GetMomentumPhoton2()->Fill(mom_phot2);
    if(myPhotID.GetCut_PhotonMom(1)) myHistManager.GetMomentumPhoton1_cut()->Fill(mom_phot1);
    if(myPhotID.GetCut_PhotonMom(2)) myHistManager.GetMomentumPhoton2_cut()->Fill(mom_phot2);

    myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(mom_phot1, beta_phot1);
    myHistManager.GetBeta_Recalc()->Fill(beta_phot1,3);
    myHistManager.GetDBeta_VS_Momentum(3)->Fill(mom_phot1, dBeta_phot1);
    myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(mom_phot2, beta_phot2);
    myHistManager.GetBeta_Recalc()->Fill(beta_phot2,4);
    myHistManager.GetDBeta_VS_Momentum(4)->Fill(mom_phot2, dBeta_phot2);

    myPhotID.SetCut_PhotonBeta(beta_phot1,1);
    myPhotID.SetCut_PhotonBeta(beta_phot2,2);
    myPhotID.SetCut_PhotonBeta_All();
    if(myPhotID.GetCut_PhotonBeta_All()) myCounter.Increment("Photon ID (Beta)");

    myHistManager.GetBetaPhoton1()->Fill(beta_phot1);
    myHistManager.GetBetaPhoton2()->Fill(beta_phot2);
    if(myPhotID.GetCut_PhotonBeta(1)) myHistManager.GetBetaPhoton1_cut()->Fill(beta_phot1);
    if(myPhotID.GetCut_PhotonBeta(2)) myHistManager.GetBetaPhoton2_cut()->Fill(beta_phot2);

    myPhotID.SetCut_PhotonECu(ecu_phot1,1);
    myPhotID.SetCut_PhotonECu(ecu_phot2,2);
    myPhotID.SetCut_PhotonECv(ecv_phot1,1);
    myPhotID.SetCut_PhotonECv(ecv_phot2,2);
    myPhotID.SetCut_PhotonECw(ecw_phot1,1);
    myPhotID.SetCut_PhotonECw(ecw_phot2,2);
    myPhotID.SetCut_PhotonECfid(1);
    myPhotID.SetCut_PhotonECfid(2);
    myPhotID.SetCut_PhotonECfid_All();
    if(myPhotID.GetCut_PhotonECfid_All()) myCounter.Increment("Photon ID (ECfid)");

    myHistManager.GetECuPhoton1()->Fill(ecu_phot1);
    myHistManager.GetECuPhoton2()->Fill(ecu_phot2);
    myHistManager.GetECvPhoton1()->Fill(ecv_phot1);
    myHistManager.GetECvPhoton2()->Fill(ecv_phot2);
    myHistManager.GetECwPhoton1()->Fill(ecw_phot1);
    myHistManager.GetECwPhoton2()->Fill(ecw_phot2);
    if(myPhotID.GetCut_PhotonECu(1)) myHistManager.GetECuPhoton1_cut()->Fill(ecu_phot1);
    if(myPhotID.GetCut_PhotonECu(2)) myHistManager.GetECuPhoton2_cut()->Fill(ecu_phot2);
    if(myPhotID.GetCut_PhotonECv(1)) myHistManager.GetECvPhoton1_cut()->Fill(ecv_phot1);
    if(myPhotID.GetCut_PhotonECv(2)) myHistManager.GetECvPhoton2_cut()->Fill(ecv_phot2);
    if(myPhotID.GetCut_PhotonECw(1)) myHistManager.GetECwPhoton1_cut()->Fill(ecw_phot1);
    if(myPhotID.GetCut_PhotonECw(2)) myHistManager.GetECwPhoton2_cut()->Fill(ecw_phot2);

    myECgeom.Put_UVW(ecu_phot1,ecv_phot1,ecw_phot1);
    myHistManager.GetEC_XvsY_local_Sector_Photon1(Sector_index_phot1 - 1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    if(myPhotID.GetCut_PhotonECfid(1)){
        myHistManager.GetEC_XvsY_local_FidCut_Photon1(Sector_index_phot1 - 1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    }else{
        myHistManager.GetEC_XvsY_local_AntiFidCut_Photon1(Sector_index_phot1 - 1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    }

    myECgeom.Put_UVW(ecu_phot2,ecv_phot2,ecw_phot2);
    myHistManager.GetEC_XvsY_local_Sector_Photon2(Sector_index_phot2 - 1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    if(myPhotID.GetCut_PhotonECfid(2)){
        myHistManager.GetEC_XvsY_local_FidCut_Photon2(Sector_index_phot2 - 1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    }else{
        myHistManager.GetEC_XvsY_local_AntiFidCut_Photon2(Sector_index_phot2 - 1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    }

    myHistManager.GetECtimePhoton1()->Fill(ectime_phot1);
    myHistManager.GetECtimePhoton2()->Fill(ectime_phot2);
    myHistManager.GetECpathPhoton1()->Fill(ecpath_phot1);
    myHistManager.GetECpathPhoton2()->Fill(ecpath_phot2);
    myHistManager.GetECpathtimePhoton1()->Fill(ecpath_phot1/LIGHTSPEED);
    myHistManager.GetECpathtimePhoton2()->Fill(ecpath_phot2/LIGHTSPEED);

    myPhotID.SetCut_PhotonTiming(timing_phot1,1);
    myPhotID.SetCut_PhotonTiming(timing_phot2,2);
    myPhotID.SetCut_PhotonTiming_All();
    if(myPhotID.GetCut_PhotonTiming_All()) myCounter.Increment("Photon ID (Timing)");

    myHistManager.GetECtime_ECl_Photon1()->Fill(timing_phot1);
    myHistManager.GetECtime_ECl_Photon2()->Fill(timing_phot2);

    myHistManager.GetECtime_ECl_Start_Photon1()->Fill(timing_phot1 - eventStartTime);
    myHistManager.GetECtime_ECl_Start_Photon2()->Fill(timing_phot2 - eventStartTime);

    if(myPhotID.GetCut_PhotonTiming(1)) myHistManager.GetECtime_ECl_Photon1_cut()->Fill(ectime_phot1 - ecpath_phot1/LIGHTSPEED);
    if(myPhotID.GetCut_PhotonTiming(2)) myHistManager.GetECtime_ECl_Photon2_cut()->Fill(ectime_phot2 - ecpath_phot2/LIGHTSPEED);

    myPhotID.SetCut_PhotonECinTimesECout(ecin_phot1,ecout_phot1,1);
    myPhotID.SetCut_PhotonECinTimesECout(ecin_phot2,ecout_phot2,2);
    myPhotID.SetCut_PhotonECinTimesECout_All();

    myHistManager.GetECtotP_vs_P_Photon1()->Fill(mom_phot1,ectot_phot1/mom_phot1);
    myHistManager.GetECtotP_vs_P_Photon2()->Fill(mom_phot2,ectot_phot2/mom_phot2);
    myHistManager.GetECin_vs_ECout_Photon1()->Fill(ecin_phot1,ecout_phot1);
    myHistManager.GetECin_vs_ECout_Photon2()->Fill(ecin_phot2,ecout_phot2);

    if(myPhotID.GetCut_PhotonECinTimesECout(1)){
        myHistManager.GetECtotP_vs_P_InOutZeroCut_Photon1()->Fill(mom_phot1,ectot_phot1/mom_phot1);
        myHistManager.GetECin_vs_ECout_InOutZeroCut_Photon1()->Fill(ecin_phot1,ecout_phot1);
    }

    if(myPhotID.GetCut_PhotonECinTimesECout(2)){
        myHistManager.GetECtotP_vs_P_InOutZeroCut_Photon2()->Fill(mom_phot2,ectot_phot2/mom_phot2);
        myHistManager.GetECin_vs_ECout_InOutZeroCut_Photon2()->Fill(ecin_phot2,ecout_phot2);
    }

    myPhotID.SetCut_PhotonSCMsq(scMassSq_phot1,1);
    myPhotID.SetCut_PhotonSCMsq(scMassSq_phot2,2);
    myPhotID.SetCut_PhotonSCMsq_All();
    if(myPhotID.GetCut_PhotonSCMsq_All()) myCounter.Increment("Photon ID (TOF Msq)");
    
    myPhotID.SetCut_Photon_All();
    if(myPhotID.GetCut_Photon_All()){
        myCounter.Increment("Photon ID (All)");
        myHistManager.GetEChit_M2_VS_scMsq(3)->Fill(echit_M2_phot1,scMassSq_phot1);
        myHistManager.GetEChit_M3_VS_scMsq(3)->Fill(echit_M3_phot1,scMassSq_phot1);
        myHistManager.GetEChit_M4_VS_scMsq(3)->Fill(echit_M4_phot1,scMassSq_phot1);

        myHistManager.GetEChit_M2_VS_scMsq(4)->Fill(echit_M2_phot2,scMassSq_phot2);
        myHistManager.GetEChit_M3_VS_scMsq(4)->Fill(echit_M3_phot2,scMassSq_phot2);
        myHistManager.GetEChit_M4_VS_scMsq(4)->Fill(echit_M4_phot2,scMassSq_phot2);
    }
    
    return myPhotID.GetCut_Photon_All();
}

bool Analyze_Electron(PartReader elecReader, double Qsq, double targMass)
{
    ElectronID myElecID;
    EC_geometry myECgeom;
    
    bool cuts_ElecID;
    
    double MASS_ELECTRON = myHistManager.GetMASS_ELECTRON();
    
    double emMom = elecReader.Get_Mom();
    double emECtot = elecReader.Get_ECtot();
    double emECin = elecReader.Get_ECin();
    double emECout = elecReader.Get_ECout();
    double emECu = elecReader.Get_ECu();
    double emECv = elecReader.Get_ECv();
    double emECw = elecReader.Get_ECw();
    double emEChit_M2 = elecReader.Get_EChit_M2();
    double emEChit_M3 = elecReader.Get_EChit_M3();
    double emEChit_M4 = elecReader.Get_EChit_M4();
    double emCCnphe = elecReader.Get_CCnphe();
    double emSCMassSq = elecReader.Get_TOF_MassSquared();
    double emECtime = elecReader.Get_ECtime();
    double emSCtime = elecReader.Get_SCtime();
    double emECpath = elecReader.Get_ECpath();
    double emSCpath = elecReader.Get_SCpath();
    double emdt = emECtime - emSCtime - ECSC_TimeOffset;
    
    TLorentzVector electron = elecReader.GetLorentzVector(MASS_ELECTRON);
    TVector3 elec_vert = elecReader.GetVertex();
    int Sector_index = GetSectorByPhi(electron.Phi());
    
    double emBeta = elecReader.Get_Beta();
    myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(emMom, emBeta);
    myHistManager.GetBeta_Recalc()->Fill(emBeta,0);
    
    double emdBeta = elecReader.Get_BetaDifference(MASS_ELECTRON); // difference in beta for measured and ideal beta
    myHistManager.GetDBeta_VS_Momentum(0)->Fill(emMom, emdBeta);
    
    myElecID.SetCut_ElecMom(emMom); // e- momentum cut
    myElecID.SetCut_ElecCCnphe(emCCnphe); // e- CC nphe cut
    myElecID.SetCut_ElecECoverP(emMom,emECtot,Sector_index,targMass); // e- EC total energy vs momentum cut
    myElecID.SetCut_ElecECin(emECin); // e- EC inner cut
    myElecID.SetCut_Elec_dtECSC(emdt); // e- timing cut
    myElecID.SetCut_ElecECfid(emECu, emECv, emECw); // e- fiducial cuts
    myElecID.SetCut_Elec_All(); // all e- ID cuts
    
    myECgeom.Put_UVW(emECu,emECv,emECw);
    myHistManager.GetEC_XvsY_local_Sector(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    myHistManager.GetECtotP_VS_P_Sector(Sector_index-1)->Fill(emMom,emECtot/emMom);
    myHistManager.GetECinP_VS_ECoutP(Sector_index-1)->Fill(emECin/emMom, emECout/emMom);
    
    myElecID.SetCut_ElecECinP_VS_ECoutP(emMom, emECin, emECout, Sector_index);
    if(myElecID.GetCut_ElecECinP_VS_ECoutP()) {
        myHistManager.GetECinP_VS_ECoutP_cut(Sector_index-1)->Fill(emECin/emMom, emECout/emMom);
    }
    
    if(emMom > 0.5 && emMom <= 1.0) {
        myHistManager.GetECinP_VS_ECoutP_Range(0)->Fill(emECin/emMom, emECout/emMom);
    } else if(emMom > 1.0 && emMom <= 1.5) {
        myHistManager.GetECinP_VS_ECoutP_Range(1)->Fill(emECin/emMom, emECout/emMom);
    } else if(emMom > 1.5 && emMom <= 2.0) {
        myHistManager.GetECinP_VS_ECoutP_Range(2)->Fill(emECin/emMom, emECout/emMom);
    } else if(emMom > 2.0 && emMom <= 2.5) {
        myHistManager.GetECinP_VS_ECoutP_Range(3)->Fill(emECin/emMom, emECout/emMom);
    } else if(emMom > 2.5 && emMom <= 3.0) {
        myHistManager.GetECinP_VS_ECoutP_Range(4)->Fill(emECin/emMom, emECout/emMom);
    }
    
    if(myElecID.GetCut_ElecMom()) myCounter.Increment("Electron ID (Mom)");
    if(myElecID.GetCut_ElecCCnphe()) myCounter.Increment("Electron ID (CCnphe)");
    if(myElecID.GetCut_Elec_dtECSC()) myCounter.Increment("Electron ID (dtECSC)");
    if(myElecID.GetCut_ElecECin()) myCounter.Increment("Electron ID (ECin)");
    if(myElecID.GetCut_ElecECoverP()){
        myCounter.Increment("Electron ID (ECPvsP)");
        myHistManager.GetECtotP_VS_P_ECPCut(Sector_index-1)->Fill(emMom,emECtot/emMom);
    }
    if(myElecID.GetCut_ElecECfid()){
        myCounter.Increment("Electron ID (ECfid)");
        myHistManager.GetECin_VS_ECout_ECfid()->Fill(emECin,emECout);
        myHistManager.GetEC_XvsY_local_FidCut(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    }else{
        myHistManager.GetEC_XvsY_local_AntiFidCut(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    }
    
    if(myElecID.GetCut_Elec_All()){
        myCounter.Increment("Electron ID (All)");
        myHistManager.GetECin_VS_ECout_elecID_All()->Fill(emECin,emECout);
        myHistManager.GetEChit_M2_VS_scMsq(0)->Fill(emEChit_M2,emSCMassSq);
        myHistManager.GetEChit_M3_VS_scMsq(0)->Fill(emEChit_M3,emSCMassSq);
        myHistManager.GetEChit_M4_VS_scMsq(0)->Fill(emEChit_M4,emSCMassSq);
    }
    if (emECout < 0.01){
        myHistManager.GetBeta_VS_Momentum_ECoutCut()->Fill(emMom, emBeta);
        myHistManager.GetTheta_VS_Phi_ECoutCut()->Fill(electron.Theta() * TMath::RadToDeg(), electron.Phi() * TMath::RadToDeg());
        myHistManager.GetElecZVert_ECoutCut()->Fill(elec_vert.Z());
        myHistManager.GetQ2_ECoutCut()->Fill(Qsq);
        
        myHistManager.GetECtot_VS_P_ECoutCut()->Fill(emMom,emECtot);
        myHistManager.GetECtotP_VS_P_ECoutCut()->Fill(emMom,emECtot/emMom);
        myHistManager.GetECtotMinusECin_ECoutCut()->Fill(emECtot-emECin);
        
        myHistManager.GetEC_XvsY_local_ECoutCut(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
    }
    else {
        myHistManager.GetBeta_VS_Momentum_AntiECoutCut()->Fill(emMom, emBeta);
        myHistManager.GetTheta_VS_Phi_AntiECoutCut()->Fill(electron.Theta() * TMath::RadToDeg(), electron.Phi() * TMath::RadToDeg());
        myHistManager.GetElecZVert_AntiECoutCut()->Fill(elec_vert.Z());
        myHistManager.GetQ2_AntiECoutCut()->Fill(Qsq);
        
        myHistManager.GetECtot_VS_P_AntiECoutCut()->Fill(emMom,emECtot);
        myHistManager.GetECtotP_VS_P_AntiECoutCut()->Fill(emMom,emECtot/emMom);
        myHistManager.GetECtotMinusECin_AntiECoutCut()->Fill(emECtot-emECin);
    }
    
    // Testing the electron ID
    for(int ii=0; ii<myElecID.Get_nElecID(); ii++){
        cuts_ElecID = false; // intialize the cuts
        
        if (myElecID.Get_elecIDLabel(ii).compare("No Cuts")==0) {
            cuts_ElecID = true;
        }else if (myElecID.Get_elecIDLabel(ii).compare("Momentum")==0) {
            cuts_ElecID = myElecID.Check_ElecMom(emMom);
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
            cuts_ElecID = myElecID.Check_ElecECoverP(emMom,emECtot,Sector_index,targMass);
        }else{
            cuts_ElecID = false;
        }
        
        if(cuts_ElecID){
            myHistManager.GetCCnphe_elecID()->Fill(emCCnphe,ii);
            myHistManager.GetMom_elecID()->Fill(emMom,ii);
            myHistManager.GetECu_elecID()->Fill(emECu,ii);
            myHistManager.GetECv_elecID()->Fill(emECv,ii);
            myHistManager.GetECw_elecID()->Fill(emECw,ii);
            myHistManager.GetDtime_ECSC_elecID()->Fill(emdt,ii);
            myHistManager.GetECtot_VS_P_elecID(ii)->Fill(emMom,emECtot);
            myHistManager.GetECtotP_VS_P_elecID(ii)->Fill(emMom,emECtot/emMom);
            myHistManager.GetECin_VS_ECout_elecID(ii)->Fill(emECin,emECout);
            myHistManager.GetMom_VS_ECout_elecID(ii)->Fill(emMom,emECout);
            myHistManager.GetECu_VS_ECout_elecID(ii)->Fill(emECu,emECout);
            myHistManager.GetECv_VS_ECout_elecID(ii)->Fill(emECv,emECout);
            myHistManager.GetECw_VS_ECout_elecID(ii)->Fill(emECw,emECout);
        }
    }
    return myElecID.GetCut_Elec_All();
}

void PrintUsage(char *processName)
{
    cerr << processName << " <options> <filename>\n";
    cerr << "\toptions are:\n";
    cerr << "\t-o<filename>\tROOT output file (def. = ctProcess_omega.root).\n";
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
    bool printCuts = true; // print the cut parameters
    
    string inFile;
    string outFile = "ctProcess_omega.root";

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
  
    myCounter.Init(); // zero out the counters
    
    myHistManager.BookHist(); // declare histograms
    
    for (i = optind; i < argc; ++i) {
        inFile = argv[i]; // process all arguments on command line.
        if (inFile != '-') { // we have a file to process
            cout << "Analyzing file " << inFile << endl; // let user know which file is being processed
            // process the root file and return number of processed events
            TotEvents = process(inFile,MaxEvents,dEvents,targMass,printCuts);
            cout<<TotEvents<<" events processed"<<endl; // print out stats
            printCuts = false;
        }
    }
    
    myHistManager.WriteHist(outFile); // write histograms to a file
    
    myCounter.Print(); // Print the counter statistics
    
    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);
}
#endif
