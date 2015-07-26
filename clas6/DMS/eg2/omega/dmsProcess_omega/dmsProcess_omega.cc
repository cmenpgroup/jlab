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
    
	double TwoPhotonAngle, elecPhoton1Angle, elecPhoton2Angle;
    double Qsq, nu, Mx, z_fracEnergy, W;
    double sinHalfTheta;
    double partMom;
    double timeEC, timeSC, pathEC, pathSC;
    double dt_ECminusSC[5];
    
    double emECu, emECv, emECw, emECin, emECout, emECtot, emCCnphe, emdt; // variables for electron id cuts
    double eventStartTime; // event start time from HEAD bank
    
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
        
        if (!(processed % dEvents)) cout << "Processed Entries: " << processed << endl;
        if (DEBUG) reader.printEvent();
        
        reader.readEntry(processed);
        
        // HEAD bank info
        eventStartTime = reader.getStartTime(); // evetn start time
        StartTime->Fill(eventStartTime);
        
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
            elecZVertSector->Fill(elec_vert.Z(),Sector_index);
            elecZVertSector_Corr->Fill(elec_vert_corr.Z(),Sector_index);
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
        elecZVert_VS_Phi_Corr->Fill(elec.Phi() * TMath::RadToDeg(),elec_vert_corr.Z()); // fill electron z vertex vs phi histogram

        
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
     
        // Use EVNT to use the automatic eg2 cuts
        if(reader.getBankRows("EVNT")){
            cuts_photID = true;
            ElecID_All = true;
        }
        
        // Use EXPB to select the PID cuts
        if (reader.getBankRows("EXPB")) {
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
        
            ElecID_Mom = myElecID.Check_ElecMom(elec.P()); // e- momentum cut
            ElecID_ECvsP = myElecID.Check_ElecECoverP(elec.P(),emECtot,Sector_index,targMass); // e- EC total energy vs momentum cut
            ElecID_ECin = myElecID.Check_ElecECin(emECin); // e- EC inner cut
            ElecID_dtECSC = myElecID.Check_Elec_dtECSC(emdt); // e- timing cut
            ElecID_ECfid = myElecID.Check_ElecECu(emECu) && myElecID.Check_ElecECv(emECv) && myElecID.Check_ElecECw(emECw); // e- fiducial cuts
            ElecID_All = (ElecID_Mom && ElecID_ECvsP && ElecID_ECin && ElecID_dtECSC && ElecID_ECfid);  // all e- ID cuts
        
            myECgeom.Put_UVW(emECu,emECv,emECw);
            EC_XvsY_local_Sector[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());

            ECtotP_VS_P_Sector[Sector_index-1]->Fill(elec.P(),emECtot/elec.P());

            ECinP_VS_ECoutP[Sector_index-1]->Fill(emECin/elec.P(), emECout/elec.P());
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
                ECinP_VS_ECoutP_cut[Sector_index-1]->Fill(emECin/elec.P(), emECout/elec.P());
            }

            if(elec.P() > 0.5 && elec.P() <= 1.0) {
                ECinP_VS_ECoutP_Range[0]->Fill(emECin/elec.P(), emECout/elec.P());
            } else if(elec.P() > 1.0 && elec.P() <= 1.5) {
                ECinP_VS_ECoutP_Range[1]->Fill(emECin/elec.P(), emECout/elec.P());
            } else if(elec.P() > 1.5 && elec.P() <= 2.0) {
                ECinP_VS_ECoutP_Range[2]->Fill(emECin/elec.P(), emECout/elec.P());
            } else if(elec.P() > 2.0 && elec.P() <= 2.5) {
                ECinP_VS_ECoutP_Range[3]->Fill(emECin/elec.P(), emECout/elec.P());
            } else if(elec.P() > 2.5 && elec.P() <= 3.0) {
                ECinP_VS_ECoutP_Range[4]->Fill(emECin/elec.P(), emECout/elec.P());
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
                ECtotP_VS_P_ECPCut[Sector_index-1]->Fill(elec.P(),emECtot/elec.P());
            }
        
            if(ElecID_ECfid){
                ctr_elecID_ECfid++;
                ECin_VS_ECout_ECfid->Fill(emECin,emECout);
                EC_XvsY_local_FidCut[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }else{
                EC_XvsY_local_AntiFidCut[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }
        
            if(ElecID_All){
                ctr_elecID++;
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
            
                EC_XvsY_local_ECoutCut[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
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

            //
            // Start of Photon ID
            //
            cuts_photID1_mom = myPhotID.Check_PhotonMom(photon1.P());
            cuts_photID2_mom = myPhotID.Check_PhotonMom(photon2.P());
            cuts_photID_mom = cuts_photID1_mom && cuts_photID2_mom;
            if(cuts_photID_mom){
                ctr_photID_Mom++;
            }
            
            MomentumPhoton1->Fill(photon1.P());
            MomentumPhoton2->Fill(photon2.P());
            if(cuts_photID1_mom) MomentumPhoton1_cut->Fill(photon1.P());
            if(cuts_photID2_mom) MomentumPhoton2_cut->Fill(photon2.P());
            
            cuts_photID1_beta = myPhotID.Check_PhotonBeta(photon1.Beta());
            cuts_photID2_beta = myPhotID.Check_PhotonBeta(photon2.Beta());
            cuts_photID_beta = cuts_photID1_beta && cuts_photID2_beta;
            if(cuts_photID_beta){
                ctr_photID_Beta++;
            }
            
            BetaPhoton1->Fill(photon1.Beta());
            BetaPhoton2->Fill(photon2.Beta());
            if(cuts_photID1_beta) BetaPhoton1_cut->Fill(photon1.Beta());
            if(cuts_photID2_beta) BetaPhoton2_cut->Fill(photon2.Beta());

            Double_t ecu_phot1 = reader.getProperty("ecu",BankIndex_part[3]);
            Double_t ecu_phot2 = reader.getProperty("ecu",BankIndex_part[4]);
            Double_t ecv_phot1 = reader.getProperty("ecv",BankIndex_part[3]);
            Double_t ecv_phot2 = reader.getProperty("ecv",BankIndex_part[4]);
            Double_t ecw_phot1 = reader.getProperty("ecw",BankIndex_part[3]);
            Double_t ecw_phot2 = reader.getProperty("ecw",BankIndex_part[4]);

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
            
            ECuPhoton1->Fill(ecu_phot1);
            ECuPhoton2->Fill(ecu_phot2);
            ECvPhoton1->Fill(ecv_phot1);
            ECvPhoton2->Fill(ecv_phot2);
            ECwPhoton1->Fill(ecw_phot1);
            ECwPhoton2->Fill(ecw_phot2);
            if(cuts_photID1_fidu) ECuPhoton1_cut->Fill(ecu_phot1);
            if(cuts_photID2_fidu) ECuPhoton2_cut->Fill(ecu_phot2);
            if(cuts_photID1_fidv) ECvPhoton1_cut->Fill(ecv_phot1);
            if(cuts_photID2_fidv) ECvPhoton2_cut->Fill(ecv_phot2);
            if(cuts_photID1_fidw) ECwPhoton1_cut->Fill(ecw_phot1);
            if(cuts_photID2_fidw) ECwPhoton2_cut->Fill(ecw_phot2);

            myECgeom.Put_UVW(ecu_phot1,ecv_phot1,ecw_phot1);
            EC_XvsY_local_Sector_Photon1[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            if(cuts_photID1_fid){
                EC_XvsY_local_FidCut_Photon1[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }else{
                EC_XvsY_local_AntiFidCut_Photon1[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }
            
            myECgeom.Put_UVW(ecu_phot2,ecv_phot2,ecw_phot2);
            EC_XvsY_local_Sector_Photon2[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            if(cuts_photID2_fid){
                EC_XvsY_local_FidCut_Photon2[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }else{
                EC_XvsY_local_AntiFidCut_Photon2[Sector_index-1]->Fill(myECgeom.Get_Xlocal(),myECgeom.Get_Ylocal());
            }
            
            Double_t ectime_phot1 = reader.getProperty("ectime",BankIndex_part[3]);
            Double_t ectime_phot2 = reader.getProperty("ectime",BankIndex_part[4]);
            Double_t ecpath_phot1 = reader.getProperty("ecpath",BankIndex_part[3]);
            Double_t ecpath_phot2 = reader.getProperty("ecpath",BankIndex_part[4]);

            ECtimePhoton1->Fill(ectime_phot1);
            ECtimePhoton2->Fill(ectime_phot2);
            ECpathPhoton1->Fill(ecpath_phot1);
            ECpathPhoton2->Fill(ecpath_phot2);
            ECpathtimePhoton1->Fill(ecpath_phot1/LIGHTSPEED);
            ECpathtimePhoton2->Fill(ecpath_phot2/LIGHTSPEED);
            
            timing_phot1 = ectime_phot1 - ecpath_phot1/LIGHTSPEED;
            timing_phot2 = ectime_phot2 - ecpath_phot2/LIGHTSPEED;
            
            cuts_photID1_time = myPhotID.Check_PhotonTiming(timing_phot1,1);
            cuts_photID2_time = myPhotID.Check_PhotonTiming(timing_phot2,2);
            cuts_photID_time = cuts_photID1_time && cuts_photID2_time;
            if(cuts_photID_time){
                ctr_photID_timing++;
            }
            
            ECtime_ECl_Photon1->Fill(timing_phot1);
            ECtime_ECl_Photon2->Fill(timing_phot2);
            
            ECtime_ECl_Start_Photon1->Fill(timing_phot1 - eventStartTime);
            ECtime_ECl_Start_Photon2->Fill(timing_phot2 - eventStartTime);
            
            if(cuts_photID1_time) ECtime_ECl_Photon1_cut->Fill(ectime_phot1 - ecpath_phot1/LIGHTSPEED);
            if(cuts_photID2_time) ECtime_ECl_Photon2_cut->Fill(ectime_phot2 - ecpath_phot2/LIGHTSPEED);

            Double_t ecin_phot1 = reader.getProperty("ecin",BankIndex_part[3]);
            Double_t ecin_phot2 = reader.getProperty("ecin",BankIndex_part[4]);
            Double_t ecout_phot1 = reader.getProperty("ecout",BankIndex_part[3]);
            Double_t ecout_phot2 = reader.getProperty("ecout",BankIndex_part[4]);
            Double_t ectot_phot1 = reader.getProperty("ectot",BankIndex_part[3]);
            Double_t ectot_phot2 = reader.getProperty("ectot",BankIndex_part[4]);

            cuts_photID1_ECinTimesECout = myPhotID.Check_PhotonECinTimesECout(ecin_phot1,ecout_phot1);
            cuts_photID2_ECinTimesECout = myPhotID.Check_PhotonECinTimesECout(ecin_phot2,ecout_phot2);
            cuts_photID_ECinTimesECout = cuts_photID1_ECinTimesECout && cuts_photID2_ECinTimesECout;
            
            ECtotP_vs_P_Photon1->Fill(photon1.P(),ectot_phot1/photon1.P());
            ECtotP_vs_P_Photon2->Fill(photon2.P(),ectot_phot2/photon2.P());
            ECin_vs_ECout_Photon1->Fill(ecin_phot1,ecout_phot1);
            ECin_vs_ECout_Photon2->Fill(ecin_phot2,ecout_phot2);

            if(cuts_photID1_ECinTimesECout){
                ECtotP_vs_P_InOutZeroCut_Photon1->Fill(photon1.P(),ectot_phot1/photon1.P());
                ECin_vs_ECout_InOutZeroCut_Photon1->Fill(ecin_phot1,ecout_phot1);
            }

            if(cuts_photID2_ECinTimesECout){
                ECtotP_vs_P_InOutZeroCut_Photon2->Fill(photon2.P(),ectot_phot2/photon2.P());
                ECin_vs_ECout_InOutZeroCut_Photon2->Fill(ecin_phot2,ecout_phot2);
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
            Pi0Mass_PhotIDcuts->Fill(TwoPhoton.M());
        }
        
        if(cuts_photID && ElecID_All && cuts_pi0fit_mass) {
            OmegaMass_AllCuts->Fill(Omega.M());
        }
        
        //
        // Start omega ID
        //
        if(ElecID_All && cuts_photID){
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
                ctr_omegaID_Mpi0++;
                OpAng_VS_E_MassPi0Cut[Vz_index]->Fill(TwoPhoton.E(),TwoPhotonAngle);
                IM2Photons[Vz_index]->Fill(TwoPhoton.M(),1);
                IMOmega[Vz_index]->Fill(Omega.M(),1);
            }else{
                IMOmega_antiCut[Vz_index]->Fill(Omega.M(),1);
            }
            cuts_woPi0Mass = (cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutBetaPhoton && cutW);
            if(cuts_woPi0Mass){ // applying all cuts except the pi0 mass
                IMOmega_woCut[Vz_index]->Fill(Omega.M(),1);
                IM2Photons[Vz_index]->Fill(TwoPhoton.M(),8);
            }
            
            if(cutQSquared) { // applying the Q^2 cut
                ctr_omegaID_Q2++;
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
                ctr_omegaID_W++;
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
                ctr_omegaID_Zmatch++;
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
                ctr_omegaID_OpAngElecPhoton++;
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
            cutsAll = (cutPi0Mass && cutZDiff && cutQSquared && cutOpAng_ElecPhoton && cutW);
            if(cutsAll){
                ctr_omegaID++;
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
        }
        //
        // End omega ID
        //
        
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

    sprintf(hname,"StartTime");
    sprintf(htitle,"Event Start Time");
    StartTime = new TH1D(hname,htitle, 420, -40., 100.);
    
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

    sprintf(hname,"elecZVertSector");
    sprintf(htitle,"Z Vertex of Electron vs Sector");
    elecZVertSector = new TH2D(hname, htitle, 300, -40, -10, 6, 0.5, 6.5);

    sprintf(hname,"elecZVertSector_Corr");
    sprintf(htitle,"Z Vertex of Electron vs Sector, Corrected");
    elecZVertSector_Corr = new TH2D(hname, htitle, 300, -40, -10, 6, 0.5, 6.5);
    
    sprintf(hname,"elecZVert_VS_Phi");
    sprintf(htitle,"Z Vertex  vs. #phi, Electrons");
    elecZVert_VS_Phi = new TH2D(hname,htitle, 360, -180., 180., 300, -35., -20.);

    sprintf(hname,"elecZVert_VS_Phi_Corr");
    sprintf(htitle,"Z Vertex (corrected) vs. #phi, Electrons");
    elecZVert_VS_Phi_Corr = new TH2D(hname,htitle, 360, -180., 180., 300, -35., -20.);

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
    
    for(i=0; i < 5; i++) {
        sprintf(hname, "ECinP_VS_ECout/P_%.1f_P_%.1f",i*0.5+0.5, i*0.5+1.0);
        sprintf(htitle,"ECin/P vs ECout/P, %.1f < P <= %.1f",i*0.5+0.5, i*0.5+1.0);
        ECinP_VS_ECoutP_Range[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.5);
    }

    for(i=0; i<MAX_SECTORS; i++){
        sprintf(hname,"ECtotP_VS_P_Sector%i",i+1);
        sprintf(htitle,"ECtot/P vs P, Sector %i",i+1);
        ECtotP_VS_P_Sector[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);

        sprintf(hname,"ECinP_VS_ECoutP_Sector%i",i+1);
        sprintf(htitle,"ECin/P vs ECout/P, Sector %i",i+1);
        ECinP_VS_ECoutP[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.5);

        sprintf(hname,"ECinP_VS_ECoutP_cut_Sector%i",i+1);
        sprintf(htitle,"ECin/P vs ECout/P cut, Sector %i",i+1);
        ECinP_VS_ECoutP_cut[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.5);

        sprintf(hname,"ECtotP_VS_P_ECPCut%i",i+1);
        sprintf(htitle,"ECtot/P vs P, Sector %i",i+1);
        ECtotP_VS_P_ECPCut[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);
        
        sprintf(hname,"EC_XvsY_local_Sector%i",i+1);
        sprintf(htitle,"EC local X vs local Y, Sector %i",i+1);
        EC_XvsY_local_Sector[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_ECoutCut%i",i+1);
        sprintf(htitle,"EC local X vs local Y, EC_{out} cut, Sector %i",i+1);
        EC_XvsY_local_ECoutCut[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);

        sprintf(hname,"EC_XvsY_local_FidCut%i",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut, Sector %i",i+1);
        EC_XvsY_local_FidCut[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);

        sprintf(hname,"EC_XvsY_local_AntiFidCut%i",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut (anti), Sector %i",i+1);
        EC_XvsY_local_AntiFidCut[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
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

    sprintf(hname,"MomentumPhoton1");
    sprintf(htitle,"Total Momentum of Photon 1");
	MomentumPhoton1 = new TH1D(hname,htitle, 500, 0, 5);

    sprintf(hname,"MomentumPhoton2");
    sprintf(htitle,"Total Momentum of Photon 2");
	MomentumPhoton2 = new TH1D(hname,htitle, 500, 0, 5);

    sprintf(hname,"MomentumPhoton1_cut");
    sprintf(htitle,"Total Momentum of Photon 1 Cut");
	MomentumPhoton1_cut = new TH1D(hname,htitle, 500, 0, 5);

    sprintf(hname,"MomentumPhoton2_cut");
    sprintf(htitle,"Total Momentum of Photon 2 Cut");
	MomentumPhoton2_cut = new TH1D(hname,htitle, 500, 0, 5);

    sprintf(hname,"BetaPhoton1");
    sprintf(htitle,"Beta of Photon 1");
	BetaPhoton1 = new TH1D(hname,htitle, 100, 0.8, 2.1);

    sprintf(hname,"BetaPhoton2");
    sprintf(htitle,"Beta of Photon 2");
	BetaPhoton2 = new TH1D(hname,htitle, 100, 0.8, 2.1);

    sprintf(hname,"BetaPhoton1_cut");
    sprintf(htitle,"Beta of Photon 1 Cut");
	BetaPhoton1_cut = new TH1D(hname,htitle, 100, 0.8, 2.1);

    sprintf(hname,"BetaPhoton2_cut");
    sprintf(htitle,"Beta of Photon 2 Cut");
	BetaPhoton2_cut = new TH1D(hname,htitle, 100, 0.8, 2.1);

    sprintf(hname,"ECuPhoton1");
    sprintf(htitle,"ECu of Photon 1");
	ECuPhoton1 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECuPhoton2");
    sprintf(htitle,"ECu of Photon 2");
	ECuPhoton2 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECvPhoton1");
    sprintf(htitle,"ECv of Photon 1");
	ECvPhoton1 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECvPhoton2");
    sprintf(htitle,"ECv of Photon 2");
	ECvPhoton2 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECwPhoton1");
    sprintf(htitle,"ECw of Photon 1");
	ECwPhoton1 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECwPhoton2");
    sprintf(htitle,"ECw of Photon 2");
	ECwPhoton2 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECuPhoton1_cut");
    sprintf(htitle,"ECu of Photon 1 Cut");
	ECuPhoton1_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECuPhoton2_cut");
    sprintf(htitle,"ECu of Photon 2 Cut");
	ECuPhoton2_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECvPhoton1_cut");
    sprintf(htitle,"ECv of Photon 1 Cut");
	ECvPhoton1_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECvPhoton2_cut");
    sprintf(htitle,"ECv of Photon 2 Cut");
	ECvPhoton2_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECwPhoton1_cut");
    sprintf(htitle,"ECw of Photon 1 Cut");
	ECwPhoton1_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECwPhoton2_cut");
    sprintf(htitle,"ECw of Photon 2 Cut");
	ECwPhoton2_cut = new TH1D(hname,htitle, 450, 0, 450);

    for(i=0; i<MAX_SECTORS; i++){
        sprintf(hname,"EC_XvsY_local_Sector_%i_Photon1",i+1);
        sprintf(htitle,"EC local X vs local Y, Sector %i, Photon 1",i+1);
        EC_XvsY_local_Sector_Photon1[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_FidCut%i_Photon1",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut, Sector %i, Photon 1",i+1);
        EC_XvsY_local_FidCut_Photon1[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_AntiFidCut%i_Photon1",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut (anti), Sector %i, Photon1",i+1);
        EC_XvsY_local_AntiFidCut_Photon1[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);

        sprintf(hname,"EC_XvsY_local_Sector_%i_Photon2",i+1);
        sprintf(htitle,"EC local X vs local Y, Sector %i, Photon 2",i+1);
        EC_XvsY_local_Sector_Photon2[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_FidCut%i_Photon2",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut, Sector %i, Photon 2",i+1);
        EC_XvsY_local_FidCut_Photon2[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_AntiFidCut%i_Photon2",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut (anti), Sector %i, Photon2",i+1);
        EC_XvsY_local_AntiFidCut_Photon2[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
    }
    
    sprintf(hname,"ECtime_ECl_Start_Photon1");
    sprintf(htitle,"ECtime - StartTime - EClength/c of Photon 1");
    ECtime_ECl_Start_Photon1 = new TH1D(hname,htitle, 400, -80., 20.);
    
    sprintf(hname,"ECtime_ECl_Start_Photon2");
    sprintf(htitle,"ECtime - StartTime - EClength/c of Photon 2");
    ECtime_ECl_Start_Photon2 = new TH1D(hname,htitle, 400, -80., 20.);
    
    sprintf(hname,"ECtime_ECl_Photon1");
    sprintf(htitle,"ECtime - EClength/c of Photon 1");
    ECtime_ECl_Photon1 = new TH1D(hname,htitle, 200, -3, 2);

    sprintf(hname,"ECtime_ECl_Photon2");
    sprintf(htitle,"ECtime - EClength/c of Photon 2");
    ECtime_ECl_Photon2 = new TH1D(hname,htitle, 200, -3, 2);

    sprintf(hname,"ECtime_ECl_Photon1_cut");
    sprintf(htitle,"ECtime - EClength/c of Photon 1 Cut");
    ECtime_ECl_Photon1_cut = new TH1D(hname,htitle, 200, -3, 2);

    sprintf(hname,"ECtime_ECl_Photon2_cut");
    sprintf(htitle,"ECtime - EClength/c of Photon 2 Cut");
    ECtime_ECl_Photon2_cut = new TH1D(hname,htitle, 200, -3, 2);

    sprintf(hname,"ECtimePhoton1");
    sprintf(htitle,"ECtime of Photon 1");
    ECtimePhoton1 = new TH1D(hname,htitle, 300, 10, 25);

    sprintf(hname,"ECtimePhoton2");
    sprintf(htitle,"ECtime of Photon 2");
    ECtimePhoton2 = new TH1D(hname,htitle, 300, 10, 25);

    sprintf(hname,"ECpathPhoton1");
    sprintf(htitle,"ECpath of Photon 1");
    ECpathPhoton1 = new TH1D(hname,htitle, 200, 450, 650);

    sprintf(hname,"ECpathPhoton2");
    sprintf(htitle,"ECpath of Photon 2");
    ECpathPhoton2 = new TH1D(hname,htitle, 200, 450, 650);

    sprintf(hname,"ECpathtimePhoton1");
    sprintf(htitle,"ECpath/c of Photon 1");
    ECpathtimePhoton1 = new TH1D(hname,htitle, 300, 17, 20);
    
    sprintf(hname,"ECpathtimePhoton2");
    sprintf(htitle,"ECpath/c of Photon 2");
    ECpathtimePhoton2 = new TH1D(hname,htitle, 300, 17, 20);
    
    sprintf(hname,"ECtotP_vs_P_Photon1");
    sprintf(htitle,"ECtot / P vs P of Photon 1");
    ECtotP_vs_P_Photon1 = new TH2D(hname,htitle, 500, 0, 5, 500, 0, 0.5);

    sprintf(hname,"ECtotP_vs_P_Photon2");
    sprintf(htitle,"ECtot / P vs P of Photon 2");
    ECtotP_vs_P_Photon2 = new TH2D(hname,htitle, 500, 0, 5, 500, 0, 0.5);

    sprintf(hname,"ECin_vs_ECout_Photon1");
    sprintf(htitle,"EC_{in} vs EC_{out} of Photon 1");
    ECin_vs_ECout_Photon1 = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.35);

    sprintf(hname,"ECin_vs_ECout_Photon2");
    sprintf(htitle,"EC_{in} vs EC_{out} of Photon 2");
    ECin_vs_ECout_Photon2 = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.35);

    sprintf(hname,"ECtotP_vs_P_InOutZeroCut_Photon1");
    sprintf(htitle,"ECtot / P vs P of Photon 1");
    ECtotP_vs_P_InOutZeroCut_Photon1 = new TH2D(hname,htitle, 500, 0, 5, 500, 0, 0.5);
    
    sprintf(hname,"ECtotP_vs_P_InOutZeroCut_Photon2");
    sprintf(htitle,"ECtot / P vs P of Photon 2");
    ECtotP_vs_P_InOutZeroCut_Photon2 = new TH2D(hname,htitle, 500, 0, 5, 500, 0, 0.5);
    
    sprintf(hname,"ECin_vs_ECout_InOutZeroCut_Photon1");
    sprintf(htitle,"EC_{in} vs EC_{out} of Photon 1");
    ECin_vs_ECout_InOutZeroCut_Photon1 = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.35);
    
    sprintf(hname,"ECin_vs_ECout_InOutZeroCut_Photon2");
    sprintf(htitle,"EC_{in} vs EC_{out} of Photon 2");
    ECin_vs_ECout_InOutZeroCut_Photon2 = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.35);
    
    sprintf(hname,"Recon_Pi0_Mass_PhotID_Cuts");
    sprintf(htitle,"Reconstructed #pi^{0} Mass - Photon ID Cuts");
    Pi0Mass_PhotIDcuts = new TH1D(hname,htitle, 100, 0, 1);

    sprintf(hname,"Recon_Omega_Mass_All_Cuts");
    sprintf(htitle,"Reconstructed #omega Mass - All Cuts");
    OmegaMass_AllCuts = new TH1D(hname,htitle,100,0,3);
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

    StartTime->GetXaxis()->SetTitle("Event Start Time (ns)");
    StartTime->GetYaxis()->SetTitle("Counts");
    StartTime->Write();
    
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

    elecZVertSector->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVertSector->GetYaxis()->SetTitle("Sector");
    elecZVertSector->Write();

    elecZVertSector_Corr->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVertSector_Corr->GetYaxis()->SetTitle("Sector");
    elecZVertSector_Corr->Write();
    
    elecZVert_VS_Phi->GetXaxis()->SetTitle("#phi (deg.)");
    elecZVert_VS_Phi->GetYaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_VS_Phi->Write();

    elecZVert_VS_Phi_Corr->GetXaxis()->SetTitle("#phi (deg.)");
    elecZVert_VS_Phi_Corr->GetYaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_VS_Phi_Corr->Write();
    
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
    }
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        ECtotP_VS_P[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P[i]->Write();
    }

    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
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

    for(i=0; i<MAX_SECTORS; i++){
        ECtotP_VS_P_Sector[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P_Sector[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P_Sector[i]->Write();
    }

    for(i=0; i<MAX_SECTORS; i++){
        ECinP_VS_ECoutP[i]->GetXaxis()->SetTitle("EC_{in}/Mom.");
        ECinP_VS_ECoutP[i]->GetYaxis()->SetTitle("EC_{out}/Mom.");
        ECinP_VS_ECoutP[i]->Write();
        ECinP_VS_ECoutP_cut[i]->GetXaxis()->SetTitle("EC_{in}/Mom.");
        ECinP_VS_ECoutP_cut[i]->GetYaxis()->SetTitle("EC_{out}/Mom.");
        ECinP_VS_ECoutP_cut[i]->Write();
    }

    for(i=0; i<5; i++) {
        ECinP_VS_ECoutP_Range[i]->GetXaxis()->SetTitle("EC_{in}/Mom.");
        ECinP_VS_ECoutP_Range[i]->GetYaxis()->SetTitle("EC_{out}/Mom.");
        ECinP_VS_ECoutP_Range[i]->Write();
    }

    for(i=0; i<MAX_SECTORS; i++){
        ECtotP_VS_P_ECPCut[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P_ECPCut[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P_ECPCut[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_Sector[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_Sector[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_Sector[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_ECoutCut[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_ECoutCut[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_ECoutCut[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_FidCut[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_FidCut[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_FidCut[i]->Write();
    }

    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_AntiFidCut[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_AntiFidCut[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_AntiFidCut[i]->Write();
    }

    // create a directory for photon id
    TDirectory *cdPhotID = out->mkdir("PhotonID");
    cdPhotID->cd();
    
    MomentumPhoton1->GetXaxis()->SetTitle("P (GeV/c)");
    MomentumPhoton1->GetYaxis()->SetTitle("Counts");
	MomentumPhoton1->Write();

    MomentumPhoton2->GetXaxis()->SetTitle("P (GeV/c)");
    MomentumPhoton2->GetYaxis()->SetTitle("Counts");
	MomentumPhoton2->Write();

    MomentumPhoton1_cut->GetXaxis()->SetTitle("P (GeV/c)");
    MomentumPhoton1_cut->GetYaxis()->SetTitle("Counts");
	MomentumPhoton1_cut->Write();

    MomentumPhoton2_cut->GetXaxis()->SetTitle("P (GeV/c)");
    MomentumPhoton2_cut->GetYaxis()->SetTitle("Counts");
	MomentumPhoton2_cut->Write();

    BetaPhoton1->GetXaxis()->SetTitle("#beta");
    BetaPhoton1->GetYaxis()->SetTitle("Counts");
	BetaPhoton1->Write();

    BetaPhoton2->GetXaxis()->SetTitle("#beta");
    BetaPhoton2->GetYaxis()->SetTitle("Counts");
	BetaPhoton2->Write();

    BetaPhoton1_cut->GetXaxis()->SetTitle("#beta");
    BetaPhoton1_cut->GetYaxis()->SetTitle("Counts");
	BetaPhoton1_cut->Write();

    BetaPhoton2_cut->GetXaxis()->SetTitle("#beta");
    BetaPhoton2_cut->GetYaxis()->SetTitle("Counts");
	BetaPhoton2_cut->Write();

    ECuPhoton1->GetXaxis()->SetTitle("EC U (cm)");
    ECuPhoton1->GetYaxis()->SetTitle("Counts");
	ECuPhoton1->Write();

    ECuPhoton2->GetXaxis()->SetTitle("EC U (cm)");
    ECuPhoton2->GetYaxis()->SetTitle("Counts");
	ECuPhoton2->Write();

    ECvPhoton1->GetXaxis()->SetTitle("EC V (cm)");
    ECvPhoton1->GetYaxis()->SetTitle("Counts");
	ECvPhoton1->Write();

    ECvPhoton2->GetXaxis()->SetTitle("EC V (cm)");
    ECvPhoton2->GetYaxis()->SetTitle("Counts");
	ECvPhoton2->Write();

    ECwPhoton1->GetXaxis()->SetTitle("EC W (cm)");
    ECwPhoton1->GetYaxis()->SetTitle("Counts");
	ECwPhoton1->Write();

    ECwPhoton2->GetXaxis()->SetTitle("EC W (cm)");
    ECwPhoton2->GetYaxis()->SetTitle("Counts");
	ECwPhoton2->Write();

    ECuPhoton1_cut->GetXaxis()->SetTitle("EC U (cm)");
    ECuPhoton1_cut->GetYaxis()->SetTitle("Counts");
	ECuPhoton1_cut->Write();

    ECuPhoton2_cut->GetXaxis()->SetTitle("EC U (cm)");
    ECuPhoton2_cut->GetYaxis()->SetTitle("Counts");
	ECuPhoton2_cut->Write();

    ECvPhoton1_cut->GetXaxis()->SetTitle("EC V (cm)");
    ECvPhoton1_cut->GetYaxis()->SetTitle("Counts");
	ECvPhoton1_cut->Write();

    ECvPhoton2_cut->GetXaxis()->SetTitle("EC V (cm)");
    ECvPhoton2_cut->GetYaxis()->SetTitle("Counts");
	ECvPhoton2_cut->Write();

    ECwPhoton1_cut->GetXaxis()->SetTitle("EC W (cm)");
    ECwPhoton1_cut->GetYaxis()->SetTitle("Counts");
	ECwPhoton1_cut->Write();

    ECwPhoton2_cut->GetXaxis()->SetTitle("EC W (cm)");
    ECwPhoton2_cut->GetYaxis()->SetTitle("Counts");
	ECwPhoton2_cut->Write();

    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_Sector_Photon1[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_Sector_Photon1[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_Sector_Photon1[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_FidCut_Photon1[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_FidCut_Photon1[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_FidCut_Photon1[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_AntiFidCut_Photon1[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_AntiFidCut_Photon1[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_AntiFidCut_Photon1[i]->Write();
    }

    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_Sector_Photon2[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_Sector_Photon2[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_Sector_Photon2[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_FidCut_Photon2[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_FidCut_Photon2[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_FidCut_Photon2[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_AntiFidCut_Photon2[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_AntiFidCut_Photon2[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_AntiFidCut_Photon2[i]->Write();
    }
    
    ECtime_ECl_Start_Photon1->GetXaxis()->SetTitle("t_{EC} - t_{start} - l_{EC}/c (cm)");
    ECtime_ECl_Start_Photon1->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Start_Photon1->Write();
    
    ECtime_ECl_Start_Photon2->GetXaxis()->SetTitle("t_{EC} - t_{start} - l_{EC}/c (cm)");
    ECtime_ECl_Start_Photon2->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Start_Photon2->Write();
    
    ECtime_ECl_Photon1->GetXaxis()->SetTitle("t_{EC} - l_{EC}/c  (ns)");
    ECtime_ECl_Photon1->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Photon1->Write();

    ECtime_ECl_Photon2->GetXaxis()->SetTitle("t_{EC} - l_{EC}/c  (ns)");
    ECtime_ECl_Photon2->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Photon2->Write();

    ECtime_ECl_Photon1_cut->GetXaxis()->SetTitle("t_{EC} - l_{EC}/c (ns)");
    ECtime_ECl_Photon1_cut->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Photon1_cut->Write();

    ECtime_ECl_Photon2_cut->GetXaxis()->SetTitle("t_{EC} - l_{EC}/c  (ns)");
    ECtime_ECl_Photon2_cut->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Photon2_cut->Write();

    ECtimePhoton1->GetXaxis()->SetTitle("EC_{time} (ns)");
    ECtimePhoton1->GetYaxis()->SetTitle("Counts");
    ECtimePhoton1->Write();

    ECtimePhoton2->GetXaxis()->SetTitle("EC_{time} (ns)");
    ECtimePhoton2->GetYaxis()->SetTitle("Counts");
    ECtimePhoton2->Write();

    ECpathPhoton1->GetXaxis()->SetTitle("EC_{path} (cm)");
    ECpathPhoton1->GetYaxis()->SetTitle("Counts");
    ECpathPhoton1->Write();

    ECpathPhoton2->GetXaxis()->SetTitle("EC_{path} (cm)");
    ECpathPhoton2->GetYaxis()->SetTitle("Counts");
    ECpathPhoton2->Write();

    ECpathtimePhoton1->GetXaxis()->SetTitle("EC_{path}/c (ns)");
    ECpathtimePhoton1->GetYaxis()->SetTitle("Counts");
    ECpathtimePhoton1->Write();
    
    ECpathtimePhoton2->GetXaxis()->SetTitle("EC_{path}/c (ns)");
    ECpathtimePhoton2->GetYaxis()->SetTitle("Counts");
    ECpathtimePhoton2->Write();
    
    ECtotP_vs_P_Photon1->GetXaxis()->SetTitle("P (GeV/c)");
    ECtotP_vs_P_Photon1->GetYaxis()->SetTitle("EC_{tot} / P");
    ECtotP_vs_P_Photon1->Write();

    ECtotP_vs_P_Photon2->GetXaxis()->SetTitle("P (GeV/c)");
    ECtotP_vs_P_Photon2->GetYaxis()->SetTitle("EC_{tot} / P");
    ECtotP_vs_P_Photon2->Write();

    ECin_vs_ECout_Photon1->GetXaxis()->SetTitle("EC_{in}");
    ECin_vs_ECout_Photon1->GetYaxis()->SetTitle("EC_{out}");
    ECin_vs_ECout_Photon1->Write();

    ECin_vs_ECout_Photon2->GetXaxis()->SetTitle("EC_{in}");
    ECin_vs_ECout_Photon2->GetYaxis()->SetTitle("EC_{out}");
    ECin_vs_ECout_Photon2->Write();

    ECtotP_vs_P_InOutZeroCut_Photon1->GetXaxis()->SetTitle("P (GeV/c)");
    ECtotP_vs_P_InOutZeroCut_Photon1->GetYaxis()->SetTitle("EC_{tot} / P");
    ECtotP_vs_P_InOutZeroCut_Photon1->Write();
    
    ECtotP_vs_P_InOutZeroCut_Photon2->GetXaxis()->SetTitle("P (GeV/c)");
    ECtotP_vs_P_InOutZeroCut_Photon2->GetYaxis()->SetTitle("EC_{tot} / P");
    ECtotP_vs_P_InOutZeroCut_Photon2->Write();
    
    ECin_vs_ECout_InOutZeroCut_Photon1->GetXaxis()->SetTitle("EC_{in}");
    ECin_vs_ECout_InOutZeroCut_Photon1->GetYaxis()->SetTitle("EC_{out}");
    ECin_vs_ECout_InOutZeroCut_Photon1->Write();
    
    ECin_vs_ECout_InOutZeroCut_Photon2->GetXaxis()->SetTitle("EC_{in}");
    ECin_vs_ECout_InOutZeroCut_Photon2->GetYaxis()->SetTitle("EC_{out}");
    ECin_vs_ECout_InOutZeroCut_Photon2->Write();
    
    // create a directory for reconstructed particle cuts
    TDirectory *cdReconC = out->mkdir("ReconCuts");
    cdReconC->cd();

    Pi0Mass_PhotIDcuts->GetXaxis()->SetTitle("Reconstructed #pi^{0} mass (GeV/c^{2})");
    Pi0Mass_PhotIDcuts->GetYaxis()->SetTitle("Counts");
    Pi0Mass_PhotIDcuts->Write();

    OmegaMass_AllCuts->GetXaxis()->SetTitle("Reconstructed #omega mass (GeV/c^{2})");
    OmegaMass_AllCuts->GetYaxis()->SetTitle("Counts");
    OmegaMass_AllCuts->Write();
    
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
