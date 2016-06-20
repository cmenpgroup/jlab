#include "Ana_pythiaCPP_Omega.h"

void PrintUsage(char *processName)
{
    cerr << processName << " <options> <filename>\n";
    cerr << "\toptions are:\n";
    cerr << "\t-o<filename>\tROOT output file (def. = ctProcess_omega.root).\n";
    cerr << "\t-M#\t\tprocess maximum # of events.\n";
    cerr << "\t-D#\t\tinform user when # of events have been processed (def. = 1000).\n";
    cerr << "\t-A#\t\tAnalysis type (def. = 0)\n";
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

int main (int argc, char **argv) {
    
    extern char *optarg;
    int c;
    extern int optind;
    
    int i;
    
    int dEvents = 1000; // increment of events for processing print statement
    int MaxEvents = 0; // max. number of events to process
    int iAna = 0; // analysis type
    int totOmegas = 0;
    
    string inFile;
    string outFile = "Ana_pythiaCPP_Omega.root";
    
    float timeStart = clock(); // start time
    
    for (i = 0; i < argc; ++i) cerr << argv[i] << " "; cerr << endl;
    while ((c = getopt(argc,argv, "o:M:D:T:ih")) != -1 ) {
        switch (c) {
            case 'o': outFile = optarg; break;
            case 'M': MaxEvents = atoi(optarg); break;
            case 'D': dEvents = atoi(optarg); break;
            case 'A': iAna = atoi(optarg); break;
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
    
    myHistManager.BookHist();
    
    for (i = optind; i < argc; ++i) {
        inFile = argv[i]; // process all arguments on command line.
        if (inFile != '-') { // we have a file to process
            cout << "Analyzing file " << inFile << endl; // let user know which file is being processed
            // process the root file and return number of processed events
            totOmegas = process(inFile,iAna,MaxEvents,dEvents);
            cout<<"=========================="<<endl;
            cout<<totOmegas<<" omegas found"<<endl; // print out stats
            cout<<"=========================="<<endl;
        }
    }
    
    myHistManager.WriteHist(outFile); // write the histograms to file
    
    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);
}

int process (string inFile, int iAna, int MaxEvents, int dEvents) {
    
    int i, j;
    int ii, jj;
    int iPC; // particle combination index
    int iPhoton;
    
    int totOmegas = 0;
    int nPip, nPim, nPi0, nPhoton, nOmega, nEta, nEtaPrime, nPhi; // particle counters per event
    
    TLorentzVector beam;
    TLorentzVector target;
    TLorentzVector electron;
    TLorentzVector recoil;
    TLorentzVector BeamMinusElectron;
    TLorentzVector missing;
    TLorentzVector W_TLV;
    TLorentzVector pip;
    TLorentzVector pim;
    TLorentzVector pi0;
    TLorentzVector photon[2];
    TLorentzVector twoPhotons;
    TLorentzVector omega;
    TLorentzVector PipPim;
    TLorentzVector PipPimPi0;
    TLorentzVector PipPim2Photons;
    TLorentzVector X;
    TLorentzVector Xrecoil;
    TLorentzVector Xpip;
    TLorentzVector Xpim;
    TLorentzVector Xpi0;
    
    Float_t Qsq, W;
    
    int ntNpart, ntTargetA, ntTargetZ, ntProcess;
    int ntks[MAX_TRACKS], ntPID[MAX_TRACKS], ntParent[MAX_TRACKS];
    double ntEbeam, ntNu;
    double ntPx[MAX_TRACKS], ntPy[MAX_TRACKS], ntPz[MAX_TRACKS], ntP[MAX_TRACKS], ntE[MAX_TRACKS];
    
    double ChargedPionAngle;
    double ThreePionAngle;
    double TwoPhotonAngle;
    
    vector<int> tempPid;
    vector<int> tempParent;
    vector<int> tempKs;
    
    bool cuts_Qsq;
    bool cuts_W;
    bool cuts_beam;
    bool cuts_target;
    bool cuts_electron;
    bool cuts_recoil;
    bool no_other_mesons;
    
    // open text files for writing events with certain particle combinations
    char OutFile[100];
    sprintf(OutFile,"PartList09.dat");
    std::ofstream fout09(OutFile);
    
    sprintf(OutFile,"PartList10.dat");
    std::ofstream fout10(OutFile);
    
    sprintf(OutFile,"PartList11.dat");
    std::ofstream fout11(OutFile);
    
    sprintf(OutFile,"PartList12.dat");
    std::ofstream fout12(OutFile);
    
    PDG_pythiaCPP_Omega myPDG; // declare the PDG object
    
    TFile *fm = new TFile(inFile.c_str(),"READ"); // data file containing the trees
    
    TTree *myTree = (TTree*)fm->Get("MC");
    
    int nEntries = (int)myTree->GetEntries();
    
    cout<<"Number of entries = "<<nEntries<<endl;
    
    if(MaxEvents==0) MaxEvents = nEntries;
    
    myTree->SetBranchAddress("Ntracks",&ntNpart);
    myTree->SetBranchAddress("Eb",&ntEbeam);
    myTree->SetBranchAddress("tarA",&ntTargetA);
    myTree->SetBranchAddress("tarZ",&ntTargetZ);
    myTree->SetBranchAddress("process",&ntProcess);
    myTree->SetBranchAddress("nu",&ntNu);
    
    myTree->SetBranchAddress("ks",&ntks);
    myTree->SetBranchAddress("type",&ntPID);
    myTree->SetBranchAddress("parent",&ntParent);
    myTree->SetBranchAddress("px",&ntPx);
    myTree->SetBranchAddress("py",&ntPy);
    myTree->SetBranchAddress("pz",&ntPz);
    myTree->SetBranchAddress("p",&ntP);
    myTree->SetBranchAddress("E",&ntE);
    
    for(ii=0; ii<MaxEvents; ii++){

        myTree->GetEntry(ii); // retrieve the event from the ntuple

        if(!(ii % dEvents)) cout<<ii<<endl;
        
        cuts_beam = false;
        cuts_target = false;
        cuts_electron = false;
        cuts_recoil = false;
        cuts_Qsq = false;
        cuts_W =false;
        no_other_mesons=false;
        
        iPhoton = 0;
        nPip = 0;
        nPim = 0;
        nPi0 = 0;
        nPhoton = 0;
        nOmega = 0;
        nEta = 0;
        nEtaPrime = 0;
        nPhi = 0;
        tempKs.clear();
        tempPid.clear();
        tempParent.clear();
        
        for (i=0; i<ntNpart; i++) {
            tempKs.push_back(ntks[i]);
            tempPid.push_back(ntPID[i]);
            tempParent.push_back(ntParent[i]);
            
            if(ntPID[i]==ID_ELECTRON && ntParent[i]==0){
                beam.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                cuts_beam = true;
            }
            
            if(ntPID[i]==ID_PROTON && ntParent[i]==0){
                target.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                cuts_target =true;
            }
            
            if(ntPID[i]==ID_PROTON && ntks[i]==1){
                recoil.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                cuts_recoil = true;
            }
            
            if(ntPID[i]==ID_ELECTRON && ntks[i]==1){
                electron.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                if(Cut_CLAS6_Theta_Ana(electron, iAna, "EC")) cuts_electron = true;
            }
            
            if(ntPID[i]==ID_PION_POS && ntks[i]==1){
                pip.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                if(Cut_CLAS6_Theta_Ana(pip, iAna, "TOF")) nPip++;
            }
            if(ntPID[i]==ID_PION_NEG && ntks[i]==1){
                pim.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                if(Cut_CLAS6_Theta_Ana(pim, iAna, "TOF")) nPim++;
            }
            if(ntPID[i]==ID_PION_ZERO && ntks[i]==11){
                pi0.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                if(Cut_CLAS6_Theta_Ana(pi0, iAna, "TOF")) nPi0++;
            }
            if(ntPID[i]==ID_PHOTON && ntks[i]==1){
                if(iPhoton<2){
                    photon[iPhoton].SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                    if(Cut_CLAS6_Theta_Ana(photon[iPhoton], iAna, "EC")) nPhoton++;
                }
                iPhoton++;
            }
            if(ntPID[i]==ID_OMEGA_MESON && ntks[i]==11){
                totOmegas++;
                nOmega++;
                omega.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
            }

            if(ntPID[i]==ID_ETA_MESON && ntks[i]==11){
                nEta++;
            }

            if(ntPID[i]==ID_PHI_MESON && ntks[i]==11){
                nPhi++;
            }

            if(ntPID[i]==ID_ETA_PRIME_MESON && ntks[i]==11){
                nEtaPrime++;
            }
            
            if(ntks[i]==11) myHistManager.Get_hIntermedPart()->Fill(ntPID[i]);
        }
        
        if(nEta==0 && nEtaPrime==0 && nPhi==0) no_other_mesons = true; // set if these mesons are not in the event
        
        if(cuts_beam && cuts_target && cuts_electron){
        
            BeamMinusElectron = beam - electron;
            
            Qsq = -1.0*BeamMinusElectron.Mag2();
            myHistManager.Get_hQsq_NoCuts()->Fill(Qsq);
        
            myHistManager.Get_hNu_NoCuts()->Fill(ntNu);
        
            W_TLV = BeamMinusElectron + target;
            W = W_TLV.M();
            myHistManager.Get_hW_NoCuts()->Fill(W);
            
            cuts_Qsq = (Qsq >= LIMIT_QSQ);
            cuts_W = (W >= LIMIT_W);
        
            if(cuts_Qsq && cuts_W){
                
                myHistManager.Get_hQsq()->Fill(Qsq);
                myHistManager.Get_hNu()->Fill(ntNu);
                myHistManager.Get_hW()->Fill(W);
                
                myHistManager.Get_hPartPerEvt()->Fill(nPip,1);
                myHistManager.Get_hPartPerEvt()->Fill(nPim,2);
                myHistManager.Get_hPartPerEvt()->Fill(nPi0,3);
                myHistManager.Get_hPartPerEvt()->Fill(nPhoton,4);
                myHistManager.Get_hPartPerEvt()->Fill(nOmega,5);
        
                if(nPip>=1 && nPim>=1 && nPi0>=1) myHistManager.Get_hPartPerEvt()->Fill(1,6);
                if(nPip==1 && nPim==1 && nPi0==1) myHistManager.Get_hPartPerEvt()->Fill(2,6);
                if(nPip>=1 && nPim>=1 && nPhoton>=2) myHistManager.Get_hPartPerEvt()->Fill(3,6);
                if(nPip==1 && nPim==1 && nPhoton==2) myHistManager.Get_hPartPerEvt()->Fill(4,6);
                if(nPip>=1 && nPim>=1 && nPi0>=1 && nOmega>=1) myHistManager.Get_hPartPerEvt()->Fill(5,6);
                if(nPip==1 && nPim==1 && nPi0==1 && nOmega==1) myHistManager.Get_hPartPerEvt()->Fill(6,6);
                if(nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega>=1) myHistManager.Get_hPartPerEvt()->Fill(7,6);
                if(nPip==1 && nPim==1 && nPhoton==2 && nOmega==1) myHistManager.Get_hPartPerEvt()->Fill(8,6);
                if(nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0) myHistManager.Get_hPartPerEvt()->Fill(9,6);
                if(nPip==1 && nPim==1 && nPi0==1 && nOmega==0) myHistManager.Get_hPartPerEvt()->Fill(10,6);
                if(nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0) myHistManager.Get_hPartPerEvt()->Fill(11,6);
                if(nPip==1 && nPim==1 && nPhoton==2 && nOmega==0) myHistManager.Get_hPartPerEvt()->Fill(12,6);
                if(nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0 && nEta==0 && nPhi==0) myHistManager.Get_hPartPerEvt()->Fill(13,6);
                if(nPip==1 && nPim==1 && nPi0==1 && nOmega==0 && nEta==0 && nPhi==0) myHistManager.Get_hPartPerEvt()->Fill(14,6);
                if(nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0 && nEta==0 && nPhi==0) myHistManager.Get_hPartPerEvt()->Fill(15,6);
                if(nPip==1 && nPim==1 && nPhoton==2 && nOmega==0 && nEta==0 && nPhi==0) myHistManager.Get_hPartPerEvt()->Fill(16,6);
            
                if(cuts_recoil){
                    X = BeamMinusElectron + target;
                    Xrecoil = X - recoil;
                    Xpip = X - pip;
                    Xpim = X - pim;
                    Xpi0 = X - pi0;
                }
      
                if(nPip>=1 && nPim>=1 && nPi0>=1){
                    PipPim = pip + pim;
                    PipPimPi0 = pip + pim + pi0;
                    FillHists_omega(pip,pim,pi0,ntNu,1);
                    if(cuts_recoil){
                        myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),1);
                        myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),1);
                        myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),1);
                        myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),1);
                        myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),1);
                    }
                    if(nOmega){
                        FillHists_omega(pip,pim,pi0,ntNu,5);
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),5);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),5);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),5);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),5);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),5);
                        }
                    }else{
                        FillHists_omega(pip,pim,pi0,ntNu,9);
                        fout09<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout09<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_PDGname(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),9);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),9);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),9);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),9);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),9);
                        }
                        
                        if(no_other_mesons) FillHists_omega(pip,pim,pi0,ntNu,13);
                    }
                }

                if(nPip==1 && nPim==1 && nPi0==1){
                    PipPim = pip + pim;
                    PipPimPi0 = pip + pim + pi0;
            
                    FillHists_omega(pip,pim,pi0,ntNu,2);
                    if(cuts_recoil){
                        myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),2);
                        myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),2);
                        myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),2);
                        myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),2);
                        myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),2);
                    }
                
                    if(nOmega){
                        FillHists_omega(pip,pim,pi0,ntNu,6);
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),6);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),6);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),6);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),6);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),6);
                        }
                    }else{
                        FillHists_omega(pip,pim,pi0,ntNu,10);
                        fout10<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout10<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_PDGname(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),10);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),10);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),10);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),10);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),10);
                        }
                        
                        if(no_other_mesons) FillHists_omega(pip,pim,pi0,ntNu,14);
                    }
                }
            
                if(nPip>=1 && nPim>=1 && nPhoton>=2){
                    PipPim = pip + pim;
                    twoPhotons = photon[0] + photon[1];
                    PipPim2Photons = pip + pim + twoPhotons;
            
                    TwoPhotonAngle = TMath::RadToDeg()*photon[0].Angle(photon[1].Vect());
            
                    FillHists_omega(pip,pim,twoPhotons,ntNu,3);
                    myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,0);
                    if(cuts_recoil){
                        myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),3);
                        myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),3);
                        myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),3);
                        myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),3);
                        myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),3);
                    }
                    
                    if(nOmega){
                        FillHists_omega(pip,pim,twoPhotons,ntNu,7);
                        myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,2);
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),7);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),7);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),7);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),7);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),7);
                        }
                    }else{
                        FillHists_omega(pip,pim,twoPhotons,ntNu,11);
                        myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,4);
                        fout11<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout11<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_PDGname(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),11);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),11);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),11);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),11);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),11);
                        }
                        
                        if(no_other_mesons) FillHists_omega(pip,pim,twoPhotons,ntNu,15);
                    }
                }
            
                if(nPip==1 && nPim==1 && nPhoton==2){
                    PipPim = pip + pim;
                    twoPhotons = photon[0] + photon[1];
                    PipPim2Photons = pip + pim + twoPhotons;
            
                    TwoPhotonAngle = TMath::RadToDeg()*photon[0].Angle(photon[1].Vect());
            
                    FillHists_omega(pip,pim,twoPhotons,ntNu,4);
                    myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,1);
                    if(cuts_recoil){
                        myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),4);
                        myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),4);
                        myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),4);
                        myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),4);
                        myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),4);
                    }
                    
                    if(nOmega){
                        FillHists_omega(pip,pim,twoPhotons,ntNu,8);
                        myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,3);
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),8);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),8);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),8);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),8);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),8);
                        }
                    }else{
                        FillHists_omega(pip,pim,twoPhotons,ntNu,12);
                        myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,5);
                        fout12<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout12<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_PDGname(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),12);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),12);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),12);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),12);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),12);
                        }
                        
                        if(no_other_mesons){
                            FillHists_omega(pip,pim,twoPhotons,ntNu,16);
/*                            if(PipPim2Photons.M()>=1.0 && PipPim2Photons.M()<1.02){
                                cout<<"Check B"<<endl;
                                for(jj=0; jj<ntNpart; jj++){
                                    cout<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_PDGname(tempPid[jj])<<endl;
                                }
                            } */
                        }
                    }
                }
            }
        }
    }
    
    fout09.close();
    fout10.close();
    fout11.close();
    fout12.close();
    
    return totOmegas;
}

void FillHists_omega(TLorentzVector pip, TLorentzVector pim, TLorentzVector pi0, double nu, int iPC){
    TLorentzVector PipPim = pip + pim;
    TLorentzVector PipPimPi0 = pip + pim + pi0;
    
    double ChargedPionAngle = TMath::RadToDeg()*pip.Angle(pim.Vect());
    double ThreePionAngle = TMath::RadToDeg()*PipPim.Angle(pi0.Vect());
    
    myHistManager.Get_hIMomega()->Fill(PipPimPi0.M(),iPC);
    myHistManager.Get_hIMomega_VS_IMPipPim(iPC-1)->Fill(PipPim.M(),PipPimPi0.M());
    myHistManager.Get_hOpAng_PipPim()->Fill(ChargedPionAngle,iPC);
    myHistManager.Get_hOpAng_PipPimPi0()->Fill(ThreePionAngle,iPC);
    myHistManager.Get_hz_fracEnergy()->Fill(PipPimPi0.E()/nu,iPC);
    myHistManager.Get_hDalitz_pip(iPC-1)->Fill((pi0+pip).M2(),(pim+pip).M2());
}

// Check that polar angle is within the geometrical acceptance of the EC.  Input angle in degrees.
bool Cut_CLAS6_Theta_EC(double theta){
    double theta_min = 10.0; // lower limit on polar angle, in degrees
    double theta_max = 60.0; // lower limit on polar angle, in degrees

    return (theta >= theta_min && theta < theta_max) ? true : false;
}

// Check that polar angle is within the geometrical acceptance of the TOF.  Input angle in degrees.
bool Cut_CLAS6_Theta_TOF(double theta){
    double theta_min = 10.0; // lower limit on polar angle, in degrees
    double theta_max = 160.0; // lower limit on polar angle, in degrees
    
    return (theta >= theta_min && theta < theta_max) ? true : false;
}

bool Cut_CLAS6_Theta_Ana(TLorentzVector V4, int iAna, string detName){
    bool ret = false;
    
    if(iAna==1){
        if(detName.compare("EC")==0){
            if(Cut_CLAS6_Theta_EC(V4.Theta()* TMath::RadToDeg())) ret = true;
        }
        if(detName.compare("TOF")==0){
            if(Cut_CLAS6_Theta_TOF(V4.Theta()* TMath::RadToDeg())) ret = true;
        }
    }else{
        ret = true;
    }
    return ret;
}

void PrintTLorentzVector(TLorentzVector TLV){
    cout <<"Px "<<TLV.Px()<<"\t";
    cout <<"Py "<<TLV.Py()<<"\t";
    cout <<"Pz "<<TLV.Pz()<<"\t";
    cout <<"E " <<TLV.E() <<"\t";
    cout <<"M " <<TLV.M() <<"\t";
    cout <<"M^2 " <<TLV.M2() <<endl;
}
