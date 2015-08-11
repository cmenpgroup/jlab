/************************************************************************/
/*  dmsProcess_omega.cc                                                 */
/*                                                                      */
/*  Created by Angelo Licastro and Andy Beiter, Canisius College        */
/*  July 2014 - Modified by M. H. Wood, Canisius College                */
/*                                                                      */
/************************************************************************/

#include "dmsProcess_omega.h"

int process (string inFile, int MaxEvents, int dEvents, int targMass) {
    int i, ii, j, k, kk;
    
    int Sector_index;
    int Vz_index;
    int BankIndex_part[5];
    
    bool cutPi0Mass;
    bool cutPipPimMass;
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
    
    // electron id cuts
    bool cuts_ElecID;
    bool ElecID_All;
    bool ElecID_Mom;
    bool ElecID_ECvsP;
    bool ElecID_ECin;
    bool ElecID_ECfid;
    bool ElecID_dtECSC;
    
    //pi0 id cuts
    bool cuts_photID1_mom;
    bool cuts_photID2_mom;
    bool cuts_photID_mom;
    bool cuts_photID1_beta;
    bool cuts_photID2_beta;
    bool cuts_photID_beta;
    bool cuts_photID1_fidu;
    bool cuts_photID2_fidu;
    bool cuts_photID1_fidv;
    bool cuts_photID2_fidv;
    bool cuts_photID1_fidw;
    bool cuts_photID2_fidw;
    bool cuts_photID1_fid;
    bool cuts_photID2_fid;
    bool cuts_photID_fid;
    bool cuts_photID1_time;
    bool cuts_photID2_time;
    bool cuts_photID_time;
    bool cuts_photID1_ECinTimesECout;
    bool cuts_photID2_ECinTimesECout;
    bool cuts_photID_ECinTimesECout;
    bool cuts_photID;
    bool cuts_pi0fit_mass;
    
    // charged pion id cuts
    bool cuts_nPion_scmsq;
    bool cuts_nPion_dbeta;
    bool cuts_pPion_scmsq;
    bool cuts_pPion_dbeta;
    bool cuts_chPion;
    
	double TwoPhotonAngle, elecPhoton1Angle, elecPhoton2Angle;
    double Qsq, nu, Mx, z_fracEnergy, W;
    double sinHalfTheta;
    double partMom;
    double timeEC, timeSC, pathEC, pathSC;
    double dt_ECminusSC[5];
    
    double pimBeta, pimSCpath, pimSCtime, pimSCMassSq, pimBetaMass, pimdBeta; // variables for pi- id cuts
    double pipBeta, pipSCpath, pipSCtime, pipSCMassSq, pipBetaMass, pipdBeta; // variables for pi- id cuts
    
    double emECu, emECv, emECw, emECin, emECout, emECtot, emCCnphe, emdt; // variables for electron id cuts
    double emECtime, emECpath, emSCtime, emSCpath; // more variables for electron id cuts
    double emBeta, emSCMassSq, emBetaMass, emdBeta; // more variables for electron id cuts
    double eventStartTime; // event start time from HEAD bank
    
    double ectime_phot1, ecpath_phot1, ecu_phot1, ecv_phot1, ecw_phot1, phot1Beta, scMassSq_phot1; //variables for photon 1 cuts
    double phot1BetaMass, phot1dBeta;
    double ectime_phot2, ecpath_phot2, ecu_phot2, ecv_phot2, ecw_phot2, phot2Beta, scMassSq_phot2; //variables for photon 2 cuts
    double phot2BetaMass, phot2dBeta;
    double timing_phot1, timing_phot2; // time difference between photon ECtime and ECpath/c
    
    int ctr_elecID_Mom = 0;
    int ctr_elecID_ECPvsP = 0;
    int ctr_elecID_ECin = 0;
    int ctr_elecID_dtECSC = 0;
    int ctr_elecID_ECfid = 0;
    int ctr_elecID = 0;
    
    int ctr_photID_Mom = 0;
    int ctr_photID_Beta = 0;
    int ctr_photID_timing = 0;
    int ctr_photID_ECfid = 0;
    int ctr_photID = 0;

    int ctr_omegaID_Mpi0 = 0;
    int ctr_omegaID_Q2 = 0;
    int ctr_omegaID_W = 0;
    int ctr_omegaID_Zmatch = 0;
    int ctr_omegaID_OpAngElecPhoton = 0;
    int ctr_omegaID = 0;
    
    EG2Target myTgt;
    EG2Cuts myCuts;
    OmegaMixedEvent myMixEvt;
    ElectronID myElecID;
    PhotonID myPhotID;
    ChargedPionID myChPionID;
    EC_geometry myECgeom;
    
    myMixEvt.Put_NumberOfEventsToMix(1); // add number of mixed event iterations
    myMixEvt.Put_OffsetOfEventsToMix(5); // add offset of the entry number for mixed events
    
    int NUM_MIXING_METHODS = myMixEvt.Get_nLabel(); // number of methods for mixing events
    int NUM_ENTRIES_OFFSET = myMixEvt.Get_NumberOfEventsToMix(); // retreive number of mixed event iterations
    int ENTRIES_OFFSET = myMixEvt.Get_OffsetOfEventsToMix(); // retrieve offset of the entry number for mixed events
    
    myCuts.Print_Cuts();
    myMixEvt.Print_Info();
    myElecID.Print_ElectronID();
    myPhotID.Print_PhotonID();
    myChPionID.Print_ChargedPionID();
    
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
        cutPipPimMass = false;
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
        
        // initialize electron id cuts
        ElecID_All = false;
        ElecID_Mom =false;
        ElecID_ECvsP = false;
        ElecID_ECin = false;
        ElecID_ECfid = false;
        ElecID_dtECSC = false;

        //initialize photon id cuts
        cuts_photID1_mom = false;
        cuts_photID2_mom = false;
        cuts_photID_mom = false;
        cuts_photID1_beta = false;
        cuts_photID2_beta = false;
        cuts_photID_beta = false;
        cuts_photID1_fidu = false;
        cuts_photID2_fidu = false;
        cuts_photID1_fidv = false;
        cuts_photID2_fidv = false;
        cuts_photID1_fidw = false;
        cuts_photID2_fidw = false;
        cuts_photID1_fid = false;
        cuts_photID2_fid = false;
        cuts_photID_fid = false;
        cuts_photID1_time = false;
        cuts_photID2_time = false;
        cuts_photID_time = false;
        cuts_photID1_ECinTimesECout = false;
        cuts_photID2_ECinTimesECout = false;
        cuts_photID_ECinTimesECout = false;
        cuts_photID = false;
        cuts_pi0fit_mass = false;
        
        // initialize charged pion cuts
        cuts_nPion_dbeta = false;
        cuts_nPion_scmsq = false;
        cuts_pPion_dbeta = false;
        cuts_pPion_scmsq = false;
        cuts_chPion = false;
        
        if (!(processed % dEvents)) cout << "Processed Entries: " << processed << endl;
        if (DEBUG) reader.printEvent();
        
        reader.readEntry(processed);
        
        // HEAD bank info
        eventStartTime = reader.getStartTime(); // evetn start time
        myHistManager.GetStartTime()->Fill(eventStartTime);
        // get the first electron lorentz vector and vertex
		TLorentzVector elec = reader.getLorentzVector(ID_ELECTRON, 0, MASS_ELECTRON);
		TVector3 elec_vert = reader.getVertex(ID_ELECTRON, 0);
        BankIndex_part[0] = reader.getIndexByPid(ID_ELECTRON, 0);
        
		//TLorentzVector prot = reader.getLorentzVector(ID_PROTON, 0, MASS_PROTON);
		//TVector3 prot_vert = reader.getVertex(ID_PROTON, 0);
        
        // get the pi- lorentz vector and vertex
        TLorentzVector nPion = reader.getLorentzVector(ID_PION_NEG, 0,MASS_PION_CHARGED);
		TVector3 nPion_vert = reader.getVertex(ID_PION_NEG, 0);
        BankIndex_part[1] = reader.getIndexByPid(ID_PION_NEG, 0);
        
        // get the pi+ lorentz vector and vertex
		TLorentzVector pPion = reader.getLorentzVector(ID_PION_POS, 0, MASS_PION_CHARGED);
		TVector3 pPion_vert = reader.getVertex(ID_PION_POS, 0);
        BankIndex_part[2] = reader.getIndexByPid(ID_PION_POS, 0);

        // get the first photon lorentz vector and vertex
		TLorentzVector photon1 = reader.getLorentzVector(ID_PHOTON, 0, MASS_PHOTON);
		TVector3 photon1_vert = reader.getVertex(ID_PHOTON, 0);
        BankIndex_part[3] = reader.getIndexByPid(ID_PHOTON, 0);

        // get the second photon lorentz vector and vertex
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
        
        // plots of z vertex difference between scattered electron and other decay particle
		myHistManager.GetZVertDiff()->Fill(elecNPionZVertDiff,1);
		myHistManager.GetZVertDiff()->Fill(elecPPionZVertDiff,2);
		myHistManager.GetZVertDiff()->Fill(elecPhoton1ZVertDiff,3);
		myHistManager.GetZVertDiff()->Fill(elecPhoton2ZVertDiff,4);
        
        // plots of x vs y vertices
        myHistManager.GetXvert()->Fill(elec_vert.X(), 0);
        myHistManager.GetXvert()->Fill(nPion_vert.X(), 1);
        myHistManager.GetXvert()->Fill(pPion_vert.X(), 2);
        myHistManager.GetXvert()->Fill(photon1_vert.X(), 3);
        myHistManager.GetXvert()->Fill(photon2_vert.X(), 4);

        myHistManager.GetYvert()->Fill(elec_vert.Y(), 0);
        myHistManager.GetYvert()->Fill(nPion_vert.Y(), 1);
        myHistManager.GetYvert()->Fill(pPion_vert.Y(), 2);
        myHistManager.GetYvert()->Fill(photon1_vert.Y(), 3);
        myHistManager.GetYvert()->Fill(photon2_vert.Y(), 4);
        
        myHistManager.GetXvert_VS_Yvert(0)->Fill(elec_vert.X(), elec_vert.Y());
		myHistManager.GetXvert_VS_Yvert(1)->Fill(nPion_vert.X(), nPion_vert.Y());
		myHistManager.GetXvert_VS_Yvert(2)->Fill(pPion_vert.X(), pPion_vert.Y());
		myHistManager.GetXvert_VS_Yvert(3)->Fill(photon1_vert.X(), photon1_vert.Y());
		myHistManager.GetXvert_VS_Yvert(4)->Fill(photon2_vert.X(), photon2_vert.Y());
        
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
     
        // Use EVNT to use the automatic eg2 cuts
        if(reader.getBankRows("EVNT")){
            cuts_photID = true;
            ElecID_All = true;
            cuts_chPion = true;
        }
        
        // Use EXPB to select the PID cuts
        if (reader.getBankRows("EXPB")) {
            //
            // Start of charged pion ID
            //
            pimSCtime = reader.getProperty("sctime",BankIndex_part[1]);
            pimSCpath = reader.getProperty("scpath",BankIndex_part[1]);
            pimBeta = (pimSCpath/pimSCtime)/LIGHTSPEED; // re-calculate beta
            myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(nPion.P(), pimBeta);
            myHistManager.GetBeta_Recalc()->Fill(pimBeta,1);
            
            pimBetaMass = Get_BetaFromMass(nPion.P(),MASS_PION_CHARGED); // calculate beta from ideal pi- mass
            pimdBeta = pimBeta - pimBetaMass; // difference in beta for measured and ideal beta
            myHistManager.GetDBeta_VS_Momentum(1)->Fill(nPion.P(), pimdBeta);
            
            cuts_nPion_dbeta = myChPionID.Check_ChargedPionDiffBeta(pimdBeta); // cut on pi- beta difference
            
            pimSCMassSq = Get_scMassSquared(nPion.P(),pimBeta); // calculate the TOF mass-squared
            myHistManager.GetScMassSquared()->Fill(pimSCMassSq,1);
            
            pipSCtime = reader.getProperty("sctime",BankIndex_part[2]);
            pipSCpath = reader.getProperty("scpath",BankIndex_part[2]);
            pipBeta = (pipSCpath/pipSCtime)/LIGHTSPEED; // re-calculate beta
            myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(pPion.P(), pipBeta);
            myHistManager.GetBeta_Recalc()->Fill(pipBeta,2);

            pipBetaMass = Get_BetaFromMass(pPion.P(),MASS_PION_CHARGED); // calculate beta from ideal pi+ mass
            pipdBeta = pipBeta - pipBetaMass; // difference in beta for measured and ideal beta
            myHistManager.GetDBeta_VS_Momentum(2)->Fill(pPion.P(), pipdBeta);
            
            cuts_pPion_dbeta = myChPionID.Check_ChargedPionDiffBeta(pipdBeta); // cut on pi+ beta difference
            
            pipSCMassSq = Get_scMassSquared(pPion.P(),pipBeta); // calculate the TOF mass-squared
            myHistManager.GetScMassSquared()->Fill(pipSCMassSq,2);
            
            cuts_chPion = cuts_nPion_dbeta && cuts_pPion_dbeta;
            //
            // End of charged pion ID
            //
            
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
        
                myHistManager.GetCCnphe()->Fill(reader.getProperty("ccnphe",BankIndex_part[ii]),ii);
                myHistManager.GetECu()->Fill(reader.getProperty("ecu",BankIndex_part[ii]),ii);
                myHistManager.GetECv()->Fill(reader.getProperty("ecv",BankIndex_part[ii]),ii);
                myHistManager.GetECw()->Fill(reader.getProperty("ecw",BankIndex_part[ii]),ii);
                myHistManager.GetECtot_VS_P(ii)->Fill(partMom,reader.getProperty("ectot",BankIndex_part[ii]));
                myHistManager.GetECtotP_VS_P(ii)->Fill(partMom,reader.getProperty("ectot",BankIndex_part[ii])/partMom);
                myHistManager.GetECin_VS_ECout(ii)->Fill(reader.getProperty("ecin",BankIndex_part[ii]),reader.getProperty("ecout",BankIndex_part[ii]));
            
                timeEC = reader.getProperty("ectime",BankIndex_part[ii]);
                timeSC = reader.getProperty("sctime",BankIndex_part[ii]);
                pathEC = reader.getProperty("ecpath",BankIndex_part[ii]);
                pathSC = reader.getProperty("scpath",BankIndex_part[ii]);
                dt_ECminusSC[ii] = timeEC - timeSC - 0.7;
                myHistManager.GetDtime_ECSC()->Fill(dt_ECminusSC[ii],ii);
            }

            // Electron ID cuts
            emECtot = reader.getProperty("ectot",BankIndex_part[0]);
            emECin = reader.getProperty("ecin",BankIndex_part[0]);
            emECout = reader.getProperty("ecout",BankIndex_part[0]);
            emECu = reader.getProperty("ecu",BankIndex_part[0]);
            emECv = reader.getProperty("ecv",BankIndex_part[0]);
            emECw = reader.getProperty("ecw",BankIndex_part[0]);
            emCCnphe = reader.getProperty("ccnphe",BankIndex_part[0]);
            emECtime = reader.getProperty("ectime",BankIndex_part[0]);
            emSCtime = reader.getProperty("sctime",BankIndex_part[0]);
            emECpath = reader.getProperty("ecpath",BankIndex_part[0]);
            emSCpath = reader.getProperty("scpath",BankIndex_part[0]);
            emdt = emECtime - emSCtime - 0.7;
        
            emBeta = (emSCpath/emSCtime)/LIGHTSPEED; // re-calculate beta
            myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(elec.P(), emBeta);
            myHistManager.GetBeta_Recalc()->Fill(emBeta,0);
            
            emBetaMass = Get_BetaFromMass(elec.P(),MASS_ELECTRON); // calculate beta from ideal pi- mass
            emdBeta = emBeta - emBetaMass; // difference in beta for measured and ideal beta
            myHistManager.GetDBeta_VS_Momentum(0)->Fill(elec.P(), emdBeta);
            
            emSCMassSq = Get_scMassSquared(elec.P(),emBeta);
            myHistManager.GetScMassSquared()->Fill(emSCMassSq,0);
            
            ElecID_Mom = myElecID.Check_ElecMom(elec.P()); // e- momentum cut
            ElecID_ECvsP = myElecID.Check_ElecECoverP(elec.P(),emECtot,Sector_index,targMass); // e- EC total energy vs momentum cut
            ElecID_ECin = myElecID.Check_ElecECin(emECin); // e- EC inner cut
            ElecID_dtECSC = myElecID.Check_Elec_dtECSC(emdt); // e- timing cut
            ElecID_ECfid = myElecID.Check_ElecECu(emECu) && myElecID.Check_ElecECv(emECv) && myElecID.Check_ElecECw(emECw); // e- fiducial cuts
            ElecID_All = (ElecID_Mom && ElecID_ECvsP && ElecID_ECin && ElecID_dtECSC && ElecID_ECfid);  // all e- ID cuts
        
            myECgeom.Put_UVW(emECu,emECv,emECw);
            myHistManager.GetEC_XvsY_local_Sector(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());

            myHistManager.GetECtotP_VS_P_Sector(Sector_index-1)->Fill(elec.P(),emECtot/elec.P());

            myHistManager.GetECinP_VS_ECoutP(Sector_index-1)->Fill(emECin/elec.P(), emECout/elec.P());
            /*
             1	0.192702	-0.586568	0.0420221	-0.104084
             2	0.2118	    -0.631934	0.0478416	-0.120601
             3	0.179616	-0.515158	0.052625	-0.149614
             4	0.198436	-0.621269	0.0387745	-0.0961029
             5	0.231849	-0.734145	0.0349574	-0.0709758
             6	0.177784	-0.451267	0.0620507	-0.17344
             */
            double aMean = 0.0, bMean = 0.0, aSigma = 0.0, bSigma = 0.0;
            double aAbove = 0.0, bAbove = 0.0, aBelow = 0.0, bBelow = 0.0;
            switch(Sector_index) {
            case 1:
                aMean = 0.192702;
                bMean = -0.586568;
                aSigma = 0.0420221;
                bSigma = -0.104084;
                break;
            case 2:
                aMean = 0.2118;
                bMean = -0.631934;
                aSigma = 0.0478416;
                bSigma = -0.120601;
                break;
            case 3:
                aMean = 0.179616;
                bMean = -0.515158;
                aSigma = 0.052625;
                bSigma = -0.149614;
                break;
            case 4:
                aMean = 0.198436;
                bMean = -0.621269;
                aSigma = 0.0387745;
                bSigma = -0.0961029;
                break;
            case 5:
                aMean = 0.231849;
                bMean = -0.734145;
                aSigma = 0.0349574;
                bSigma = -0.0709758;
                break;
            case 6:
                aMean = 0.177784;
                bMean = -0.451267;
                aSigma = 0.0620507;
                bSigma = -0.17344;
                break;
            }
            aAbove = aMean + 2 * aSigma;
            bAbove = bMean + 2 * bSigma;
            aBelow = aMean - 2 * aSigma;
            bBelow = bMean - 2 * bSigma;

            if((emECout / elec.P() - (bBelow * emECin / elec.P() + aBelow) > 0) && (emECout / elec.P() - (bAbove * emECin / elec.P() + aAbove) < 0)) {
                myHistManager.GetECinP_VS_ECoutP_cut(Sector_index-1)->Fill(emECin/elec.P(), emECout/elec.P());
            }

            if(elec.P() > 0.5 && elec.P() <= 1.0) {
                myHistManager.GetECinP_VS_ECoutP_Range(0)->Fill(emECin/elec.P(), emECout/elec.P());
            } else if(elec.P() > 1.0 && elec.P() <= 1.5) {
                myHistManager.GetECinP_VS_ECoutP_Range(1)->Fill(emECin/elec.P(), emECout/elec.P());
            } else if(elec.P() > 1.5 && elec.P() <= 2.0) {
                myHistManager.GetECinP_VS_ECoutP_Range(2)->Fill(emECin/elec.P(), emECout/elec.P());
            } else if(elec.P() > 2.0 && elec.P() <= 2.5) {
                myHistManager.GetECinP_VS_ECoutP_Range(3)->Fill(emECin/elec.P(), emECout/elec.P());
            } else if(elec.P() > 2.5 && elec.P() <= 3.0) {
                myHistManager.GetECinP_VS_ECoutP_Range(4)->Fill(emECin/elec.P(), emECout/elec.P());
            }

            if(ElecID_Mom){
                ctr_elecID_Mom++;
            }
            if(ElecID_dtECSC){
                ctr_elecID_dtECSC++;
            }
            if(ElecID_ECin){
                ctr_elecID_ECin++;
            }
            
            if(ElecID_ECvsP){
                ctr_elecID_ECPvsP++;
                myHistManager.GetECtotP_VS_P_ECPCut(Sector_index-1)->Fill(elec.P(),emECtot/elec.P());
            }
        
            if(ElecID_ECfid){
                ctr_elecID_ECfid++;
                myHistManager.GetECin_VS_ECout_ECfid()->Fill(emECin,emECout);
                myHistManager.GetEC_XvsY_local_FidCut(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }else{
                myHistManager.GetEC_XvsY_local_AntiFidCut(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }
        
            if(ElecID_All){
                ctr_elecID++;
                myHistManager.GetECin_VS_ECout_elecID_All()->Fill(emECin,emECout);
            }

            if (emECout < 0.01){
                myHistManager.GetBeta_VS_Momentum_ECoutCut()->Fill(elec.P(), elec.Beta());
                myHistManager.GetTheta_VS_Phi_ECoutCut()->Fill(elec.Theta() * TMath::RadToDeg(), elec.Phi() * TMath::RadToDeg());
                myHistManager.GetElecZVert_ECoutCut()->Fill(elec_vert.Z());
                myHistManager.GetQ2_ECoutCut()->Fill(Qsq);
            
                myHistManager.GetECtot_VS_P_ECoutCut()->Fill(elec.P(),emECtot);
                myHistManager.GetECtotP_VS_P_ECoutCut()->Fill(elec.P(),emECtot/elec.P());
                myHistManager.GetECtotMinusECin_ECoutCut()->Fill(emECtot-emECin);
            
                myHistManager.GetEC_XvsY_local_ECoutCut(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }
            else {
                myHistManager.GetBeta_VS_Momentum_AntiECoutCut()->Fill(elec.P(), elec.Beta());
                myHistManager.GetTheta_VS_Phi_AntiECoutCut()->Fill(elec.Theta() * TMath::RadToDeg(), elec.Phi() * TMath::RadToDeg());
                myHistManager.GetElecZVert_AntiECoutCut()->Fill(elec_vert.Z());
                myHistManager.GetQ2_AntiECoutCut()->Fill(Qsq);
            
                myHistManager.GetECtot_VS_P_AntiECoutCut()->Fill(elec.P(),emECtot);
                myHistManager.GetECtotP_VS_P_AntiECoutCut()->Fill(elec.P(),emECtot/elec.P());
                myHistManager.GetECtotMinusECin_AntiECoutCut()->Fill(emECtot-emECin);
            }
        
            // Testing the electron ID
            for(ii=0; ii<myElecID.Get_nElecID(); ii++){
                cuts_ElecID = false; // intialize the cuts
            
                if (myElecID.Get_elecIDLabel(ii).compare("No Cuts")==0) {
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
                    myHistManager.GetCCnphe_elecID()->Fill(emCCnphe,ii);
                    myHistManager.GetMom_elecID()->Fill(elec.P(),ii);
                    myHistManager.GetECu_elecID()->Fill(emECu,ii);
                    myHistManager.GetECv_elecID()->Fill(emECv,ii);
                    myHistManager.GetECw_elecID()->Fill(emECw,ii);
                    myHistManager.GetDtime_ECSC_elecID()->Fill(emdt,ii);
                    myHistManager.GetECtot_VS_P_elecID(ii)->Fill(elec.P(),emECtot);
                    myHistManager.GetECtotP_VS_P_elecID(ii)->Fill(elec.P(),emECtot/elec.P());
                    myHistManager.GetECin_VS_ECout_elecID(ii)->Fill(emECin,emECout);
                    myHistManager.GetMom_VS_ECout_elecID(ii)->Fill(elec.P(),emECout);
                    myHistManager.GetECu_VS_ECout_elecID(ii)->Fill(emECu,emECout);
                    myHistManager.GetECv_VS_ECout_elecID(ii)->Fill(emECv,emECout);
                    myHistManager.GetECw_VS_ECout_elecID(ii)->Fill(emECw,emECout);
                }
            }
            //
            // End of  Electron ID
            //

            //
            // Start of Photon ID
            //
            ectime_phot1 = reader.getProperty("ectime",BankIndex_part[3]);
            ectime_phot2 = reader.getProperty("ectime",BankIndex_part[4]);
            ecpath_phot1 = reader.getProperty("ecpath",BankIndex_part[3]);
            ecpath_phot2 = reader.getProperty("ecpath",BankIndex_part[4]);
            
            ecu_phot1 = reader.getProperty("ecu",BankIndex_part[3]);
            ecu_phot2 = reader.getProperty("ecu",BankIndex_part[4]);
            ecv_phot1 = reader.getProperty("ecv",BankIndex_part[3]);
            ecv_phot2 = reader.getProperty("ecv",BankIndex_part[4]);
            ecw_phot1 = reader.getProperty("ecw",BankIndex_part[3]);
            ecw_phot2 = reader.getProperty("ecw",BankIndex_part[4]);
            
            cuts_photID1_mom = myPhotID.Check_PhotonMom(photon1.P());
            cuts_photID2_mom = myPhotID.Check_PhotonMom(photon2.P());
            cuts_photID_mom = cuts_photID1_mom && cuts_photID2_mom;
            if(cuts_photID_mom){
                ctr_photID_Mom++;
            }
            
            myHistManager.GetMomentumPhoton1()->Fill(photon1.P());
            myHistManager.GetMomentumPhoton2()->Fill(photon2.P());
            if(cuts_photID1_mom) myHistManager.GetMomentumPhoton1_cut()->Fill(photon1.P());
            if(cuts_photID2_mom) myHistManager.GetMomentumPhoton2_cut()->Fill(photon2.P());
            
            phot1Beta = (ecpath_phot1/ectime_phot1)/LIGHTSPEED; // re-calculate beta
            myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(photon1.P(), phot1Beta);
            myHistManager.GetBeta_Recalc()->Fill(phot1Beta,3);
            
            phot1BetaMass = Get_BetaFromMass(photon1.P(),MASS_PHOTON); // calculate beta from ideal photon mass
            phot1dBeta = phot1Beta - phot1BetaMass; // difference in beta for measured and ideal beta
            myHistManager.GetDBeta_VS_Momentum(3)->Fill(photon1.P(), phot1dBeta);
            
            scMassSq_phot1 = Get_scMassSquared(photon1.P(),phot1Beta);
            myHistManager.GetScMassSquared()->Fill(scMassSq_phot1,3);
            
            phot2Beta = (ecpath_phot2/ectime_phot2)/LIGHTSPEED; // re-calculate beta
            myHistManager.GetBeta_VS_Momentum_Recalc()->Fill(photon2.P(), phot2Beta);
            myHistManager.GetBeta_Recalc()->Fill(phot2Beta,4);

            phot2BetaMass = Get_BetaFromMass(photon2.P(),MASS_PHOTON); // calculate beta from ideal photon mass
            phot2dBeta = phot2Beta - phot2BetaMass; // difference in beta for measured and ideal beta
            myHistManager.GetDBeta_VS_Momentum(4)->Fill(photon2.P(), phot2dBeta);
            
            scMassSq_phot2 = Get_scMassSquared(photon2.P(),phot2Beta);
            myHistManager.GetScMassSquared()->Fill(scMassSq_phot2,4);
            
            cuts_photID1_beta = myPhotID.Check_PhotonBeta(photon1.Beta());
            cuts_photID2_beta = myPhotID.Check_PhotonBeta(photon2.Beta());
            cuts_photID_beta = cuts_photID1_beta && cuts_photID2_beta;
            if(cuts_photID_beta){
                ctr_photID_Beta++;
            }
            
            myHistManager.GetBetaPhoton1()->Fill(photon1.Beta());
            myHistManager.GetBetaPhoton2()->Fill(photon2.Beta());
            if(cuts_photID1_beta) myHistManager.GetBetaPhoton1_cut()->Fill(photon1.Beta());
            if(cuts_photID2_beta) myHistManager.GetBetaPhoton2_cut()->Fill(photon2.Beta());

            cuts_photID1_fidu = myPhotID.Check_PhotonECu(ecu_phot1);
            cuts_photID2_fidu = myPhotID.Check_PhotonECu(ecu_phot2);
            cuts_photID1_fidv = myPhotID.Check_PhotonECv(ecv_phot1);
            cuts_photID2_fidv = myPhotID.Check_PhotonECv(ecv_phot2);
            cuts_photID1_fidw = myPhotID.Check_PhotonECw(ecw_phot1);
            cuts_photID2_fidw = myPhotID.Check_PhotonECw(ecw_phot2);
            cuts_photID1_fid = cuts_photID1_fidu && cuts_photID1_fidv && cuts_photID1_fidw;
            cuts_photID2_fid = cuts_photID2_fidu && cuts_photID2_fidv && cuts_photID2_fidw;
            cuts_photID_fid = cuts_photID1_fid && cuts_photID2_fid;
            if(cuts_photID_fid){
                ctr_photID_ECfid++;
            }
            
            myHistManager.GetECuPhoton1()->Fill(ecu_phot1);
            myHistManager.GetECuPhoton2()->Fill(ecu_phot2);
            myHistManager.GetECvPhoton1()->Fill(ecv_phot1);
            myHistManager.GetECvPhoton2()->Fill(ecv_phot2);
            myHistManager.GetECwPhoton1()->Fill(ecw_phot1);
            myHistManager.GetECwPhoton2()->Fill(ecw_phot2);
            if(cuts_photID1_fidu) myHistManager.GetECuPhoton1_cut()->Fill(ecu_phot1);
            if(cuts_photID2_fidu) myHistManager.GetECuPhoton2_cut()->Fill(ecu_phot2);
            if(cuts_photID1_fidv) myHistManager.GetECvPhoton1_cut()->Fill(ecv_phot1);
            if(cuts_photID2_fidv) myHistManager.GetECvPhoton2_cut()->Fill(ecv_phot2);
            if(cuts_photID1_fidw) myHistManager.GetECwPhoton1_cut()->Fill(ecw_phot1);
            if(cuts_photID2_fidw) myHistManager.GetECwPhoton2_cut()->Fill(ecw_phot2);

            myECgeom.Put_UVW(ecu_phot1,ecv_phot1,ecw_phot1);
            myHistManager.GetEC_XvsY_local_Sector_Photon1(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            if(cuts_photID1_fid){
                myHistManager.GetEC_XvsY_local_FidCut_Photon1(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }else{
                myHistManager.GetEC_XvsY_local_AntiFidCut_Photon1(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }
            
            myECgeom.Put_UVW(ecu_phot2,ecv_phot2,ecw_phot2);
            myHistManager.GetEC_XvsY_local_Sector_Photon2(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            if(cuts_photID2_fid){
                myHistManager.GetEC_XvsY_local_FidCut_Photon2(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }else{
                myHistManager.GetEC_XvsY_local_AntiFidCut_Photon2(Sector_index-1)->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }

            myHistManager.GetECtimePhoton1()->Fill(ectime_phot1);
            myHistManager.GetECtimePhoton2()->Fill(ectime_phot2);
            myHistManager.GetECpathPhoton1()->Fill(ecpath_phot1);
            myHistManager.GetECpathPhoton2()->Fill(ecpath_phot2);
            myHistManager.GetECpathtimePhoton1()->Fill(ecpath_phot1/LIGHTSPEED);
            myHistManager.GetECpathtimePhoton2()->Fill(ecpath_phot2/LIGHTSPEED);
            
            timing_phot1 = ectime_phot1 - ecpath_phot1/LIGHTSPEED;
            timing_phot2 = ectime_phot2 - ecpath_phot2/LIGHTSPEED;
            
            cuts_photID1_time = myPhotID.Check_PhotonTiming(timing_phot1,1);
            cuts_photID2_time = myPhotID.Check_PhotonTiming(timing_phot2,2);
            cuts_photID_time = cuts_photID1_time && cuts_photID2_time;
            if(cuts_photID_time){
                ctr_photID_timing++;
            }
            
            myHistManager.GetECtime_ECl_Photon1()->Fill(timing_phot1);
            myHistManager.GetECtime_ECl_Photon2()->Fill(timing_phot2);
            
            myHistManager.GetECtime_ECl_Start_Photon1()->Fill(timing_phot1 - eventStartTime);
            myHistManager.GetECtime_ECl_Start_Photon2()->Fill(timing_phot2 - eventStartTime);
            
            if(cuts_photID1_time) myHistManager.GetECtime_ECl_Photon1_cut()->Fill(ectime_phot1 - ecpath_phot1/LIGHTSPEED);
            if(cuts_photID2_time) myHistManager.GetECtime_ECl_Photon2_cut()->Fill(ectime_phot2 - ecpath_phot2/LIGHTSPEED);

            Double_t ecin_phot1 = reader.getProperty("ecin",BankIndex_part[3]);
            Double_t ecin_phot2 = reader.getProperty("ecin",BankIndex_part[4]);
            Double_t ecout_phot1 = reader.getProperty("ecout",BankIndex_part[3]);
            Double_t ecout_phot2 = reader.getProperty("ecout",BankIndex_part[4]);
            Double_t ectot_phot1 = reader.getProperty("ectot",BankIndex_part[3]);
            Double_t ectot_phot2 = reader.getProperty("ectot",BankIndex_part[4]);

            cuts_photID1_ECinTimesECout = myPhotID.Check_PhotonECinTimesECout(ecin_phot1,ecout_phot1);
            cuts_photID2_ECinTimesECout = myPhotID.Check_PhotonECinTimesECout(ecin_phot2,ecout_phot2);
            cuts_photID_ECinTimesECout = cuts_photID1_ECinTimesECout && cuts_photID2_ECinTimesECout;
            
            myHistManager.GetECtotP_vs_P_Photon1()->Fill(photon1.P(),ectot_phot1/photon1.P());
            myHistManager.GetECtotP_vs_P_Photon2()->Fill(photon2.P(),ectot_phot2/photon2.P());
            myHistManager.GetECin_vs_ECout_Photon1()->Fill(ecin_phot1,ecout_phot1);
            myHistManager.GetECin_vs_ECout_Photon2()->Fill(ecin_phot2,ecout_phot2);

            if(cuts_photID1_ECinTimesECout){
                myHistManager.GetECtotP_vs_P_InOutZeroCut_Photon1()->Fill(photon1.P(),ectot_phot1/photon1.P());
                myHistManager.GetECin_vs_ECout_InOutZeroCut_Photon1()->Fill(ecin_phot1,ecout_phot1);
            }

            if(cuts_photID2_ECinTimesECout){
                myHistManager.GetECtotP_vs_P_InOutZeroCut_Photon2()->Fill(photon2.P(),ectot_phot2/photon2.P());
                myHistManager.GetECin_vs_ECout_InOutZeroCut_Photon2()->Fill(ecin_phot2,ecout_phot2);
            }
            
            cuts_photID = cuts_photID_mom && cuts_photID_beta && cuts_photID_fid && cuts_photID_time;
//            cuts_photID = cuts_photID_mom && cuts_photID_beta && cuts_photID_fid && cuts_photID_time && cuts_photID_ECinTimesECout;
            if(cuts_photID){
                ctr_photID++;
            }
            
            //
            // End of Photon ID
            //
        } // end of if(EXPB)
        
        // SumLo = 0.0866593 and SumHi = 0.193608
        if(0.0866593 < TwoPhoton.M() && TwoPhoton.M() < 0.193608) {
            cuts_pi0fit_mass = true;
        }

        if(cuts_photID && ElecID_All) {
            myHistManager.GetPi0Mass_PhotIDcuts()->Fill(TwoPhoton.M());
        }
        
        if(cuts_photID && ElecID_All && cuts_pi0fit_mass) {
            myHistManager.GetOmegaMass_AllCuts()->Fill(Omega.M());
        }
        
        if(ElecID_All){
            myHistManager.GetScMassSquared_elecID()->Fill(emSCMassSq,0);
            myHistManager.GetScMassSquared_elecID()->Fill(pimSCMassSq,1);
            myHistManager.GetScMassSquared_elecID()->Fill(pipSCMassSq,2);
            myHistManager.GetScMassSquared_elecID()->Fill(scMassSq_phot1,3);
            myHistManager.GetScMassSquared_elecID()->Fill(scMassSq_phot2,4);
        }
        
        if(cuts_photID){
            myHistManager.GetScMassSquared_photID()->Fill(emSCMassSq,0);
            myHistManager.GetScMassSquared_photID()->Fill(pimSCMassSq,1);
            myHistManager.GetScMassSquared_photID()->Fill(pipSCMassSq,2);
            myHistManager.GetScMassSquared_photID()->Fill(scMassSq_phot1,3);
            myHistManager.GetScMassSquared_photID()->Fill(scMassSq_phot2,4);
        }

        /*
            1 - pi+ pi-
            2 - pi+ pi0
            3 - pi- pi0
        */

        TLorentzVector pipPi0 = pPion + TwoPhoton;
        TLorentzVector pimPi0 = nPion + TwoPhoton;
        if(TwoPion.M() < 0.48 || TwoPion.M() > 0.51) {
            myHistManager.GetMass2Pions_VS_massOmega_NC(0)->Fill(TwoPion.M(), Omega.M()); // no cuts
            myHistManager.GetMass2Pions_VS_massOmega_NC(1)->Fill(pipPi0.M(), Omega.M()); // no cuts
            myHistManager.GetMass2Pions_VS_massOmega_NC(2)->Fill(pimPi0.M(), Omega.M()); // no cuts
        }
        //
        // Start omega ID
        //
        if(ElecID_All && cuts_photID && cuts_chPion){
            myHistManager.GetDBeta_VS_Momentum_EPC(0)->Fill(elec.P(), emdBeta);
            myHistManager.GetDBeta_VS_Momentum_EPC(1)->Fill(nPion.P(), pimdBeta);
            myHistManager.GetDBeta_VS_Momentum_EPC(2)->Fill(pPion.P(), pipdBeta);
            myHistManager.GetDBeta_VS_Momentum_EPC(3)->Fill(photon1.P(), phot1dBeta);
            myHistManager.GetDBeta_VS_Momentum_EPC(4)->Fill(photon2.P(), phot2dBeta);

            myHistManager.GetScMassSquared_elecIDphotID()->Fill(emSCMassSq,0);
            myHistManager.GetScMassSquared_elecIDphotID()->Fill(pimSCMassSq,1);
            myHistManager.GetScMassSquared_elecIDphotID()->Fill(pipSCMassSq,2);
            myHistManager.GetScMassSquared_elecIDphotID()->Fill(scMassSq_phot1,3);
            myHistManager.GetScMassSquared_elecIDphotID()->Fill(scMassSq_phot2,4);
            
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
            myHistManager.GetOpAng_VS_IMOmega(Vz_index)->Fill(Omega.M(), TwoPhotonAngle); // variable = 2 photon opening angle

            // set the cuts
            cutPi0Mass = myCuts.Check_MassPi0(TwoPhoton.M()); // pi0 mass cut
            cutQSquared = myCuts.Check_QSquared(Qsq); // Q^2 cut
            cutW = myCuts.Check_Wcut(W); // W cut
            cutPipPimMass = myCuts.Check_MassPipPim(TwoPion.M()); // pi+ pi- inv. mass cut
            
            /*
             1 - pi+ pi-
             2 - pi+ pi0
             3 - pi- pi0
             */
            
            if(cutPipPimMass) {
                myHistManager.GetMass2Pions_VS_massOmega_EPC(0)->Fill(TwoPion.M(), Omega.M()); // electron and photon cuts
                myHistManager.GetMass2Pions_VS_massOmega_EPC(1)->Fill(pipPi0.M(), Omega.M()); // electron and photon cuts
                myHistManager.GetMass2Pions_VS_massOmega_EPC(2)->Fill(pimPi0.M(), Omega.M()); // electron and photon cuts
            }
            
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
            myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),0); // inv. mass of 2 photons
            myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons
            myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons
            myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),0); // inv. mass of pi+ pi- 2 photons

            if(cutPi0Mass) { // applying the pi0 mass cut
                ctr_omegaID_Mpi0++;
                myHistManager.GetOpAng_VS_E_MassPi0Cut(Vz_index)->Fill(TwoPhoton.E(),TwoPhotonAngle);
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),1);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),1);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),1);
            }
            cuts_woPi0Mass = (cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutBetaPhoton && cutW);
            if(cuts_woPi0Mass){ // applying all cuts except the pi0 mass
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),1);
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),8);
            }
            
            if(cutQSquared) { // applying the Q^2 cut
                ctr_omegaID_Q2++;
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),2);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),2);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),2);
            }
            cuts_woQsquared = (cutPi0Mass && cutZDiff && cutOpAng_ElecPhoton && cutBetaPhoton && cutW);
            if(cuts_woQsquared){ // applying all cuts except Q^2
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),2);
            }
            
            if(cutW){ // applying the W cut
                ctr_omegaID_W++;
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),3);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),3);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),3);
            }
            cuts_woW = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutBetaPhoton);
            if(cuts_woW){ // applying all cuts except W
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),3);
            }
            
            if(cutZDiff){ // applying the e-,pion z-vertex matching cut
                ctr_omegaID_Zmatch++;
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),4);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),4);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),4);
            }
            cuts_woZDiff = (cutPi0Mass && cutQSquared && cutOpAng_ElecPhoton && cutBetaPhoton && cutW);
            if(cuts_woZDiff){ // applying all cuts except the e-,pion z-vertex matching
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),4);
            }
            
            if(cutOpAng_ElecPhoton) { // applying the e-,photon opening angle cut
                ctr_omegaID_OpAngElecPhoton++;
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),5);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),5);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),5);
            }
            cuts_woOpAng_ElecPhoton = (cutPi0Mass && cutZDiff && cutQSquared && cutBetaPhoton && cutW);
            if(cuts_woOpAng_ElecPhoton){ // applying all cuts except the e-,photon opening angle
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),5);
            }
            
            if(cutBetaPhoton){ // applying the photon beta cut
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),6);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),6);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),6);
            }
            cuts_woBetaPhoton = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutW);
            if(cuts_woBetaPhoton){ // applying all cuts except the photon beta
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),6);
            }

            if(cutElecR){ // applying the photon beta cut
                myHistManager.GetIM2Photons(Vz_index)->Fill(TwoPhoton.M(),7);
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),7);
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),7);
            }
            cuts_woElecR = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutW && cutBetaPhoton);
            if(cuts_woElecR){ // applying all cuts except the photon beta
                myHistManager.GetIMOmega_woCut(Vz_index)->Fill(Omega.M(),7);
            }
            
             // applying all cuts
            cutsAll = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutW && cutPipPimMass);
            if(cutsAll){
                ctr_omegaID++;
                myHistManager.GetIMOmega(Vz_index)->Fill(Omega.M(),8);
                myHistManager.GetW_VS_IMOmega_AllCuts(Vz_index)->Fill(W, Omega.M()); // variable = W
                myHistManager.GetIM2Pions_VS_IMOmega_AllCuts(Vz_index)->Fill(TwoPion.M(), Omega.M()); // variable = pion pair inv. mass
                myHistManager.GetPtSq_Omega_AllCuts(Vz_index)->Fill(Omega.Perp2());

                myHistManager.GetXvert_VS_Yvert_AllCuts(Vz_index)->Fill(elec_vert.X(), elec_vert.Y());

                /*
                    1 - pi+ pi-
                    2 - pi+ pi0
                    3 - pi- pi0
                */
                if(TwoPion.M() < 0.48 || TwoPion.M() > 0.51) {
                    myHistManager.GetMass2Pions_VS_massOmega_EPOC(0)->Fill(TwoPion.M(), Omega.M()); // electron, photon, and omega cuts
                    myHistManager.GetMass2Pions_VS_massOmega_EPOC(1)->Fill(pipPi0.M(), Omega.M()); // electron, photon, and omega cuts
                    myHistManager.GetMass2Pions_VS_massOmega_EPOC(2)->Fill(pimPi0.M(), Omega.M()); // electron, photon, and omega cuts
                }
                if(cutOmegaMass){
                    myHistManager.GetXvert_VS_Yvert_Omega(Vz_index)->Fill(elec_vert.X(), elec_vert.Y());
                    myHistManager.GetPtSq_Omega_AllCuts_IMOmegaCut(Vz_index)->Fill(Omega.Perp2());
                }
                if(cutOmegaMass_sb){
                    myHistManager.GetPtSq_Omega_AllCuts_IMOmegaSBCut(Vz_index)->Fill(Omega.Perp2());
                }
            }else{
                myHistManager.GetIMOmega_antiCut(Vz_index)->Fill(Omega.M(),8);
            }
        
            for(k=0; k<NUM_ENTRIES_OFFSET; k++){ // loop over number of mixed event iterations
                for(j=0; j<NUM_MIXING_METHODS; j++){ // loop over number of mixed event methods
                    myHistManager.GetIM2Photons_ME(Vz_index)->Fill(Mass_TwoPhoton_ME[j][k],j); // inv. mass of 2 photons
                    myHistManager.GetIMOmega_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j); // inv. mass of pi+ pi- 2 photons
                    if(cutOpAng_ElecPhoton1 && cutOpAng_ElecPhoton2) {
                        myHistManager.GetIM2Photons_OpAng_ElecPhoton_Cut_ME(Vz_index)->Fill(Mass_TwoPhoton_ME[j][k],j);
                        myHistManager.GetIMOmega_OpAng_ElecPhoton_Cut_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                    if(cutPi0Mass) {
                        myHistManager.GetIMOmega_MassPi0Cut_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                    if(cutZDiff_ElectronNPion && cutZDiff_ElectronPPion){
                        myHistManager.GetIMOmega_ZVertCut_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                    if(cutQSquared) {
                        myHistManager.GetIMOmega_QsqCut_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                    if(cutsAll){
                        myHistManager.GetIMOmega_AllCuts_ME(Vz_index)->Fill(Mass_Omega_ME[j][k],j);
                    }
                }
            }
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
    myHistManager.SetMass2Pions_VS_massOmega_NCX(0);
    myHistManager.SetMass2Pions_VS_massOmega_NCX(1);
    myHistManager.SetMass2Pions_VS_massOmega_NCX(2);
    myHistManager.SetMass2Pions_VS_massOmega_EPCX(0);
    myHistManager.SetMass2Pions_VS_massOmega_EPCX(1);
    myHistManager.SetMass2Pions_VS_massOmega_EPCX(2);
    myHistManager.SetMass2Pions_VS_massOmega_EPOCX(0);
    myHistManager.SetMass2Pions_VS_massOmega_EPOCX(1);
    myHistManager.SetMass2Pions_VS_massOmega_EPOCX(2);
    cout<<endl;
    cout<<"Statistics on cuts"<<endl;
    cout<<"Electron ID "<<ctr_elecID<<" ("<<float(ctr_elecID)/float(entries)<<")"<<endl;
    cout<<"Electron ID (Mom.) "<<ctr_elecID_Mom<<" ("<<float(ctr_elecID_Mom)/float(entries)<<")"<<endl;
    cout<<"Electron ID (ECPvsP) "<<ctr_elecID_ECPvsP<<" ("<<float(ctr_elecID_ECPvsP)/float(entries)<<")"<<endl;
    cout<<"Electron ID (ECin) "<<ctr_elecID_ECin<<" ("<<float(ctr_elecID_ECin)/float(entries)<<")"<<endl;
    cout<<"Electron ID (dtECSC) "<<ctr_elecID_dtECSC<<" ("<<float(ctr_elecID_dtECSC)/float(entries)<<")"<<endl;
    cout<<"Electron ID (ECfid) "<<ctr_elecID_ECfid<<" ("<<float(ctr_elecID_ECfid)/float(entries)<<")"<<endl<<endl;
    cout<<"Photon ID "<<ctr_photID<<" ("<<float(ctr_photID)/float(entries)<<")"<<endl;
    cout<<"Photon ID (Mom.) "<<ctr_photID_Mom<<" ("<<float(ctr_photID_Mom)/float(entries)<<")"<<endl;
    cout<<"Photon ID (beta) "<<ctr_photID_Beta<<" ("<<float(ctr_photID_Beta)/float(entries)<<")"<<endl;
    cout<<"Photon ID (timing) "<<ctr_photID_timing<<" ("<<float(ctr_photID_timing)/float(entries)<<")"<<endl;
    cout<<"Photon ID (ECfid) "<<ctr_photID_ECfid<<" ("<<float(ctr_photID_ECfid)/float(entries)<<")"<<endl<<endl;
    cout<<"Omega ID "<<ctr_omegaID<<" ("<<float(ctr_omegaID)/float(entries)<<")"<<endl;
    cout<<"Omega ID (Mpi0) "<<ctr_omegaID_Mpi0<<" ("<<float(ctr_omegaID_Mpi0)/float(entries)<<")"<<endl;
    cout<<"Omega ID (OpAngElecPhoton) "<<ctr_omegaID_OpAngElecPhoton<<" ("<<float(ctr_omegaID_OpAngElecPhoton)/float(entries)<<")"<<endl;
    cout<<"Omega ID (Zmatch) "<<ctr_omegaID_Zmatch<<" ("<<float(ctr_omegaID_Zmatch)/float(entries)<<")"<<endl;
    cout<<"Omega ID (W) "<<ctr_omegaID_W<<" ("<<float(ctr_omegaID_W)/float(entries)<<")"<<endl;
    cout<<"Omega ID (Q2) "<<ctr_omegaID_Q2<<" ("<<float(ctr_omegaID_Q2)/float(entries)<<")"<<endl;
    
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

// Return the Time-of-Flight mass squared
//
//          fMom = particle momentum
//          fBeta = particle beta
//
double Get_scMassSquared(double fMom, double fBeta){

    double ret;
    double fBetaSq = fBeta*fBeta;
    
    if(fBetaSq){
        ret = fMom*fMom*(1.0-fBetaSq)/fBetaSq;
    }else{
        ret = -99.0;
    }
    return ret;
}


// Return the expected particle beta given the mass and momentum
//
//          fMom = particle momentum
//          fMass = particle mass
//
double Get_BetaFromMass(double fMom, double fMass){
    
    double ret;
    
    if(fMom){
        ret = 1.0/sqrt((fMass*fMass)/(fMom*fMom) + 1.0);
    }else{
        ret = -99.0;
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
  
 
    myHistManager.BookHist(); // declare histograms
    
    for (i = optind; i < argc; ++i) {
        inFile = argv[i]; // process all arguments on command line.
        if (inFile != '-') { // we have a file to process
            
            cout << "Analyzing file " << inFile << endl; // let user know which file is being processed
            
            // process the root file and return number of processed events
            TotEvents = process(inFile,MaxEvents,dEvents,targMass);
        }
    }
    cout<<TotEvents<<" events processed"<<endl; // print out stats
    
    myHistManager.WriteHist(outFile); // write histograms to a file
    
    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);
}
#endif
