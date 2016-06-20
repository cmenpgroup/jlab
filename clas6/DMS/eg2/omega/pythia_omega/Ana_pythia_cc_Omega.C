#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define MAX_ANA 2
#define MAX_TRACKS 50
#define LIMIT_QSQ 1.0
#define LIMIT_W 2.0

#define nCtr 6
// 1 : Number of pi+
// 2 : Number of pi-
// 3 : Number of pi0
// 4 : Number of 2 photons
// 5 : Number of omega mesons
// 6 : Number of particle combinations listed below

class PartComb
{
    vector<string> LabelPartComb;
    
public:
    PartComb();
    Int_t Get_nPartComb() {return LabelPartComb.size();};
    string Get_PartComb(int num) {return LabelPartComb[num];};
};

PartComb::PartComb()
{
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPi0>=1");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPi0==1");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPhoton>=2");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPhoton==2");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPi0>=1 && nOmega>=1");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPi0==1 && nOmega==1");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega>=1");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPhoton==2 && nOmega==1");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPi0==1 && nOmega==0");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPhoton==2 && nOmega==0");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPi0==1 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPhoton==2 && nOmega==0 && nEta==0");
}

#define nCombPhoton 8 // number of particle combination with 2 photons

//gROOT->Reset();   // start from scratch

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

char hname[100];
char htitle[100];

void BookHists();
void WriteHists(string outRoot);
string Get_PDGname(int PDGcode);
bool Cut_CLAS6_Theta_EC(Double_t theta);
bool Cut_CLAS6_Theta_TOF(Double_t theta);
bool Cut_CLAS6_Theta_Ana(TLorentzVector V4, Int_t iAna, string detName);

void FillHists_omega(TLorentzVector pip, TLorentzVector pim, TLorentzVector pi0, Double_t nu, Int_t iPC){
    TLorentzVector PipPim = pip + pim;
    TLorentzVector PipPimPi0 = pip + pim + pi0;
    
    Double_t ChargedPionAngle = TMath::RadToDeg()*pip.Angle(pim.Vect());
    Double_t ThreePionAngle = TMath::RadToDeg()*PipPim.Angle(pi0.Vect());
    
    hIMomega->Fill(PipPimPi0.M(),iPC);
    hIMomega_VS_IMPipPim[iPC-1]->Fill(PipPim.M(),PipPimPi0.M());
    hOpAng_PipPim->Fill(ChargedPionAngle,iPC);
    hOpAng_PipPimPi0->Fill(ThreePionAngle,iPC);
    hz_fracEnergy->Fill(PipPimPi0.E()/nu,iPC);
}

void Ana_pythia_cc_Omega(string inRoot, string outRoot="Ana_pythia_cc_Omega.root", Int_t iAna=0, Int_t MaxEvt=0, Int_t ModEvt=1000){

    Int_t i, j;
    Int_t ii, jj;
    Int_t iPC; // particle combination index
    Int_t iPhoton;
    
    Int_t nPip, nPim, nPi0, nPhoton, nOmega, nEta, nEtaPrime, nPhi; // particle counters per event
    
    Int_t ID_ELECTRON = 11; // PDG electron id
    Int_t ID_PHOTON = 22;  // PDG photon id
    Int_t ID_PION_POS = 211;  // PDG pi+ id
    Int_t ID_PION_NEG = -211;  // PDG pi- id
    Int_t ID_PION_ZERO = 111;  // PDG pi0 id
    Int_t ID_ETA_MESON = 221;  // PDG eta id
    Int_t ID_OMEGA_MESON = 223;  // PDG omega id
    Int_t ID_ETA_PRIME_MESON = 331;  // PDG eta' (958) id
    Int_t ID_PHI_MESON = 333;  // PDG phi id
    Int_t ID_PROTON = 2212; // PDG proton id
    
    Double_t MASS_PHOTON = 0.0; // mass of photon in GeV/c^2
    Double_t MASS_ELECTRON = 0.000511; // mass of charged pion in GeV/c^2
    Double_t MASS_PION_CHARGED = 0.138; // mass of charged pion in GeV/c^2
    Double_t MASS_PION_NEUTRAL = 0.135; // mass of neutral pion in GeV/c^2
    Double_t MASS_PROTON = 0.938; // mass of proton in GeV/c^2
    
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
    
    Int_t ntNpart, ntTargetA, ntTargetZ, ntProcess;
    Int_t ntks[MAX_TRACKS], ntPID[MAX_TRACKS], ntParent[MAX_TRACKS];
    Double_t ntEbeam, ntNu;
    Double_t ntPx[MAX_TRACKS], ntPy[MAX_TRACKS], ntPz[MAX_TRACKS], ntP[MAX_TRACKS], ntE[MAX_TRACKS];

    Double_t ChargedPionAngle;
    Double_t ThreePionAngle;
    Double_t TwoPhotonAngle;
    
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
    ofstream fout09(OutFile);

    sprintf(OutFile,"PartList10.dat");
    ofstream fout10(OutFile);
    
    sprintf(OutFile,"PartList11.dat");
    ofstream fout11(OutFile);
    
    sprintf(OutFile,"PartList12.dat");
    ofstream fout12(OutFile);
    
    BookHists();
    
    // data files contain the trees
    printf("Analyzing file %s\n",inRoot.c_str());
    TFile *fm = new TFile(inRoot.c_str(),"READ");
    
    TTree *myTree = (TTree*)fm->Get("MC");
    
    Int_t nEntries = (Int_t)myTree->GetEntries();
    
    cout<<"Number of entries = "<<nEntries<<endl;
    
    if(MaxEvt==0) MaxEvt = nEntries;
    
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
    
    for(ii=0; ii<MaxEvt; ii++){

        myTree->GetEntry(ii); // retrieve the event from the ntuple

        if(!(ii % ModEvt)) cout<<ii<<endl;
        
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
            
            if(ntks[i]==11) hIntermedPart->Fill(ntPID[i]);
        }
        
        if(nEta==0 && nEtaPrime==0 && nPhi==0) no_other_mesons = true; // set if these mesons are not in the event
        
        if(cuts_beam && cuts_target && cuts_electron){
        
            BeamMinusElectron = beam - electron;
            
            Qsq = -1.0*BeamMinusElectron.Mag2();
            hQsq_NoCuts->Fill(Qsq);
        
            hNu_NoCuts->Fill(ntNu);
        
            W_TLV = BeamMinusElectron + target;
            W = W_TLV.M();
            hW_NoCuts->Fill(W);
            
            cuts_Qsq = (Qsq >= LIMIT_QSQ);
            cuts_W = (W >= LIMIT_W);
        
            if(cuts_Qsq && cuts_W){
                
                hQsq->Fill(Qsq);
                hNu->Fill(ntNu);
                hW->Fill(W);
                
                hPartPerEvt->Fill(nPip,1);
                hPartPerEvt->Fill(nPim,2);
                hPartPerEvt->Fill(nPi0,3);
                hPartPerEvt->Fill(nPhoton,4);
                hPartPerEvt->Fill(nOmega,5);
        
                if(nPip>=1 && nPim>=1 && nPi0>=1) hPartPerEvt->Fill(1,6);
                if(nPip==1 && nPim==1 && nPi0==1) hPartPerEvt->Fill(2,6);
                if(nPip>=1 && nPim>=1 && nPhoton>=2) hPartPerEvt->Fill(3,6);
                if(nPip==1 && nPim==1 && nPhoton==2) hPartPerEvt->Fill(4,6);
                if(nPip>=1 && nPim>=1 && nPi0>=1 && nOmega>=1) hPartPerEvt->Fill(5,6);
                if(nPip==1 && nPim==1 && nPi0==1 && nOmega==1) hPartPerEvt->Fill(6,6);
                if(nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega>=1) hPartPerEvt->Fill(7,6);
                if(nPip==1 && nPim==1 && nPhoton==2 && nOmega==1) hPartPerEvt->Fill(8,6);
                if(nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0) hPartPerEvt->Fill(9,6);
                if(nPip==1 && nPim==1 && nPi0==1 && nOmega==0) hPartPerEvt->Fill(10,6);
                if(nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0) hPartPerEvt->Fill(11,6);
                if(nPip==1 && nPim==1 && nPhoton==2 && nOmega==0) hPartPerEvt->Fill(12,6);
                if(nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0 && nEta==0 && nPhi==0) hPartPerEvt->Fill(13,6);
                if(nPip==1 && nPim==1 && nPi0==1 && nOmega==0 && nEta==0 && nPhi==0) hPartPerEvt->Fill(14,6);
                if(nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0 && nEta==0 && nPhi==0) hPartPerEvt->Fill(15,6);
                if(nPip==1 && nPim==1 && nPhoton==2 && nOmega==0 && nEta==0 && nPhi==0) hPartPerEvt->Fill(16,6);
            
                if(cuts_recoil){
                    X = BeamMinusElectron + target;
                    Xrecoil = X - recoil;
                    Xpip = X - pip;
                    Xpim = X - pim;
                    Xpi0 = X - pi0;
                }
      
                
                
/*                if(nEtaPrime!=0){
                cout<<"Eta' "<<nEtaPrime<<endl;
                for(jj=0; jj<ntNpart; jj++){
                    cout<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                }
                }
*/
                if(nPip>=1 && nPim>=1 && nPi0>=1){
                    PipPim = pip + pim;
                    PipPimPi0 = pip + pim + pi0;
                    FillHists_omega(pip,pim,pi0,ntNu,1);
                    if(cuts_recoil){
                        hMMsq_X->Fill(X.Mag2(),1);
                        hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),1);
                        hMMsq_Xpip->Fill(Xpip.Mag2(),1);
                        hMMsq_Xpim->Fill(Xpim.Mag2(),1);
                        hMMsq_Xpi0->Fill(Xpi0.Mag2(),1);
                    }
                    if(nOmega){
                        FillHists_omega(pip,pim,pi0,ntNu,5);
                        if(cuts_recoil){
                            hMMsq_X->Fill(X.Mag2(),5);
                            hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),5);
                            hMMsq_Xpip->Fill(Xpip.Mag2(),5);
                            hMMsq_Xpim->Fill(Xpim.Mag2(),5);
                            hMMsq_Xpi0->Fill(Xpi0.Mag2(),5);
                        }
                    }else{
                        FillHists_omega(pip,pim,pi0,ntNu,9);
                        fout09<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout09<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            hMMsq_X->Fill(X.Mag2(),9);
                            hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),9);
                            hMMsq_Xpip->Fill(Xpip.Mag2(),9);
                            hMMsq_Xpim->Fill(Xpim.Mag2(),9);
                            hMMsq_Xpi0->Fill(Xpi0.Mag2(),9);
                        }
                        
                        if(no_other_mesons) FillHists_omega(pip,pim,pi0,ntNu,13);
                    }
                }

                if(nPip==1 && nPim==1 && nPi0==1){
                    PipPim = pip + pim;
                    PipPimPi0 = pip + pim + pi0;
            
                    FillHists_omega(pip,pim,pi0,ntNu,2);
                    if(cuts_recoil){
                        hMMsq_X->Fill(X.Mag2(),2);
                        hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),2);
                        hMMsq_Xpip->Fill(Xpip.Mag2(),2);
                        hMMsq_Xpim->Fill(Xpim.Mag2(),2);
                        hMMsq_Xpi0->Fill(Xpi0.Mag2(),2);
                    }
                
                    if(nOmega){
                        FillHists_omega(pip,pim,pi0,ntNu,6);
                        if(cuts_recoil){
                            hMMsq_X->Fill(X.Mag2(),6);
                            hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),6);
                            hMMsq_Xpip->Fill(Xpip.Mag2(),6);
                            hMMsq_Xpim->Fill(Xpim.Mag2(),6);
                            hMMsq_Xpi0->Fill(Xpi0.Mag2(),6);
                        }
                    }else{
                        FillHists_omega(pip,pim,pi0,ntNu,10);
                        fout10<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout10<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            hMMsq_X->Fill(X.Mag2(),10);
                            hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),10);
                            hMMsq_Xpip->Fill(Xpip.Mag2(),10);
                            hMMsq_Xpim->Fill(Xpim.Mag2(),10);
                            hMMsq_Xpi0->Fill(Xpi0.Mag2(),10);
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
                    hOpAng_TwoPhoton->Fill(TwoPhotonAngle,0);
                    if(cuts_recoil){
                        hMMsq_X->Fill(X.Mag2(),3);
                        hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),3);
                        hMMsq_Xpip->Fill(Xpip.Mag2(),3);
                        hMMsq_Xpim->Fill(Xpim.Mag2(),3);
                        hMMsq_Xpi0->Fill(Xpi0.Mag2(),3);
                    }
                    
                    if(nOmega){
                        FillHists_omega(pip,pim,twoPhotons,ntNu,7);
                        hOpAng_TwoPhoton->Fill(TwoPhotonAngle,2);
                        if(cuts_recoil){
                            hMMsq_X->Fill(X.Mag2(),7);
                            hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),7);
                            hMMsq_Xpip->Fill(Xpip.Mag2(),7);
                            hMMsq_Xpim->Fill(Xpim.Mag2(),7);
                            hMMsq_Xpi0->Fill(Xpi0.Mag2(),7);
                        }
                    }else{
                        FillHists_omega(pip,pim,twoPhotons,ntNu,11);
                        hOpAng_TwoPhoton->Fill(TwoPhotonAngle,4);
                        fout11<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout11<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            hMMsq_X->Fill(X.Mag2(),11);
                            hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),11);
                            hMMsq_Xpip->Fill(Xpip.Mag2(),11);
                            hMMsq_Xpim->Fill(Xpim.Mag2(),11);
                            hMMsq_Xpi0->Fill(Xpi0.Mag2(),11);
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
                    hOpAng_TwoPhoton->Fill(TwoPhotonAngle,1);
                    if(cuts_recoil){
                        hMMsq_X->Fill(X.Mag2(),4);
                        hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),4);
                        hMMsq_Xpip->Fill(Xpip.Mag2(),4);
                        hMMsq_Xpim->Fill(Xpim.Mag2(),4);
                        hMMsq_Xpi0->Fill(Xpi0.Mag2(),4);
                    }
                    
                    if(nOmega){
                        FillHists_omega(pip,pim,twoPhotons,ntNu,8);
                        hOpAng_TwoPhoton->Fill(TwoPhotonAngle,3);
                        if(cuts_recoil){
                            hMMsq_X->Fill(X.Mag2(),8);
                            hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),8);
                            hMMsq_Xpip->Fill(Xpip.Mag2(),8);
                            hMMsq_Xpim->Fill(Xpim.Mag2(),8);
                            hMMsq_Xpi0->Fill(Xpi0.Mag2(),8);
                        }
                    }else{
                        FillHists_omega(pip,pim,twoPhotons,ntNu,12);
                        hOpAng_TwoPhoton->Fill(TwoPhotonAngle,5);
                        fout12<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout12<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            hMMsq_X->Fill(X.Mag2(),12);
                            hMMsq_Xrecoil->Fill(Xrecoil.Mag2(),12);
                            hMMsq_Xpip->Fill(Xpip.Mag2(),12);
                            hMMsq_Xpim->Fill(Xpim.Mag2(),12);
                            hMMsq_Xpi0->Fill(Xpi0.Mag2(),12);
                        }
                        
                        if(no_other_mesons){
                            FillHists_omega(pip,pim,twoPhotons,ntNu,16);
/*                            if(PipPim2Photons.M()>=1.0 && PipPim2Photons.M()<1.02){
                                cout<<"Check B"<<endl;
                                for(jj=0; jj<ntNpart; jj++){
                                    cout<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                                }
                            } */
                        }
                    }
                }
            }
        }
    }

    WriteHists(outRoot); // write the histograms to file
    
    fout09.close();
    fout10.close();
    fout11.close();
    fout12.close();
}

void BookHists(){

    Int_t j;

    Int_t nTheta = 180;
    Double_t ThetaLo = 0.0;
    Double_t ThetaHi = 180.0;

    Int_t nIM2pi = 200;
    Double_t IM2piLo = 0.0;
    Double_t IM2piHi = 1.0;

    Int_t nIMomega = 125;
    Double_t IMomegaLo = 0.0;
    Double_t IMomegaHi = 2.5;

    Int_t nMMsq = 200;
    Double_t MMsqLo = 0.0;
    Double_t MMsqHi = 5.0;

    PartComb myComb;
    Int_t nComb = myComb.Get_nPartComb();
    
    cout<<"nComb "<<nComb<<endl;
    
    hIntermedPart = new TH1D("hIntermedPart","Intermediate Particle ID", 3500, 0.5, 3500.5);
    hPartPerEvt = new TH2D("hPartPerEvt", "Part. Comb. per Event",nComb+1,-0.5,nComb+0.5,nCtr,0.5,nCtr+0.5);

    hIMomega = new TH2D("hIMomega", "Mass(#omega)",nIMomega,IMomegaLo,IMomegaHi,nComb,0.5,nComb+0.5);

    hMMsq_X = new TH2D("hMMsq_X", "MM^{2}, total",nMMsq,MMsqLo,MMsqHi+5.0,nComb,0.5,nComb+0.5);
    hMMsq_Xrecoil = new TH2D("hMMsq_Xrecoil", "MM^{2} off proton",nMMsq,MMsqLo,MMsqHi,nComb,0.5,nComb+0.5);
    hMMsq_Xpip = new TH2D("hMMsq_Xpip", "MM^{2} off pi^{+}",nMMsq,MMsqLo,MMsqHi+2.0,nComb,0.5,nComb+0.5);
    hMMsq_Xpim = new TH2D("hMMsq_Xpim", "MM^{2} off pi^{-}",nMMsq,MMsqLo,MMsqHi+2.0,nComb,0.5,nComb+0.5);
    hMMsq_Xpi0 = new TH2D("hMMsq_Xpi0", "MM^{2} off pi^{0}",nMMsq,MMsqLo,MMsqHi+2.0,nComb,0.5,nComb+0.5);
    
    hOpAng_PipPim = new TH2D("hOpAng_PipPim","Opening Angle(#pi^{+},#pi^{-})",nTheta,ThetaLo,ThetaHi,nComb,0.5,nComb+0.5);

    hOpAng_PipPimPi0 = new TH2D("hOpAng_PipPimPi0","Opening Angle(#pi^{+}#pi^{-},#pi^{0})",nTheta,ThetaLo,ThetaHi,nComb,0.5,nComb+0.5);

    hQsq_NoCuts = new TH1D("hQsq_NoCuts","Q^{2}", 100, 0., 4.);
    hNu_NoCuts = new TH1D("hNu_NoCuts","#nu", 100, 0., 5.);
    hW_NoCuts = new TH1D("hW_NoCuts", "W", 250, 0., 5.);

    hQsq = new TH1D("hQsq","Q^{2}", 100, 0., 4.);
    hNu = new TH1D("hNu","#nu", 100, 0., 5.);
    hW = new TH1D("hW", "W", 250, 0., 5.);

    hz_fracEnergy = new TH2D("hz_fracEnergy", "Fractional Energy", 150, 0, 1.5,nComb,0.5,nComb+0.5);
    hOpAng_TwoPhoton = new TH2D("hOpAng_TwoPhoton","Opening Angle(#gamma,#gamma)",nTheta, ThetaLo, ThetaHi,nCombPhoton,0.5,nCombPhoton+0.5);
    
    for(j=0; j<nComb; j++){
        sprintf(hname,"hIMomega_VS_IMPipPim_%i",j);
        sprintf(htitle,"Mass(#omega) vs. Mass(#pi^{+} #pi^{-}), %s",myComb.Get_PartComb(j).c_str());
        hIMomega_VS_IMPipPim[j] = new TH2D(hname,htitle,nIM2pi,IM2piLo,IM2piHi,nIMomega,IMomegaLo,IMomegaHi);
    }
}

void WriteHists(string outRoot){

    PartComb myComb;
    Int_t nComb = myComb.Get_nPartComb();
    
    Int_t j;
    
    TFile *out = new TFile(outRoot.c_str(), "recreate");
    out->cd();

    hIntermedPart->GetXaxis()->SetTitle("Intermediate Particle ID");
    hIntermedPart->GetYaxis()->SetTitle("Counts");
    hIntermedPart->Write();

    hPartPerEvt->GetXaxis()->SetTitle("Number of Particles per Event");
    hPartPerEvt->GetYaxis()->SetTitle("Particle Combination");
    hPartPerEvt->Write();

    hIMomega->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hIMomega->GetYaxis()->SetTitle("Particle Combination");
    hIMomega->Write();

    hMMsq_X->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_X->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_X->Write();

    hMMsq_Xrecoil->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_Xrecoil->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_Xrecoil->Write();
    
    hMMsq_Xpip->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_Xpip->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_Xpip->Write();
    
    hMMsq_Xpim->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_Xpim->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_Xpim->Write();
    
    hMMsq_Xpi0->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_Xpi0->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_Xpi0->Write();
    
    hQsq_NoCuts->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    hQsq_NoCuts->GetYaxis()->SetTitle("Counts");
    hQsq_NoCuts->Write();

    hNu_NoCuts->GetXaxis()->SetTitle("\nu (GeV)");
    hNu_NoCuts->GetYaxis()->SetTitle("Counts");
    hNu_NoCuts->Write();

    hW_NoCuts->GetXaxis()->SetTitle("W (GeV)");
    hW_NoCuts->GetYaxis()->SetTitle("Counts");
    hW_NoCuts->Write();

    hQsq->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    hQsq->GetYaxis()->SetTitle("Counts");
    hQsq->Write();
    
    hNu->GetXaxis()->SetTitle("\nu (GeV)");
    hNu->GetYaxis()->SetTitle("Counts");
    hNu->Write();
    
    hW->GetXaxis()->SetTitle("W (GeV)");
    hW->GetYaxis()->SetTitle("Counts");
    hW->Write();
    
    hz_fracEnergy->GetXaxis()->SetTitle("z");
    hz_fracEnergy->GetYaxis()->SetTitle("Counts");
    hz_fracEnergy->Write();

    hOpAng_PipPim->GetXaxis()->SetTitle("#theta_{#pi^{+},#pi^{-}}(deg.)");
    hOpAng_PipPim->GetYaxis()->SetTitle("Particle Combination");
    hOpAng_PipPim->Write();

    hOpAng_PipPimPi0->GetXaxis()->SetTitle("#theta_{#pi^{+}#pi^{-},#pi^{0}}(deg.)");
    hOpAng_PipPimPi0->GetYaxis()->SetTitle("Particle Combination");
    hOpAng_PipPimPi0->Write();

    hOpAng_TwoPhoton->GetXaxis()->SetTitle("#theta_{#gamma,#gamma}(deg.)");
    hOpAng_TwoPhoton->GetYaxis()->SetTitle("Particle Combination with 2 photons");
    hOpAng_TwoPhoton->Write();
    
    for(j=0; j<nComb; j++){
        hIMomega_VS_IMPipPim[j]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} Mass (GeV/c^{2})");
        hIMomega_VS_IMPipPim[j]->GetYaxis()->SetTitle("#pi^{+} #pi^{-} #pi^{0} Mass (GeV/c^{2})");
        hIMomega_VS_IMPipPim[j]->Write();
    }

    out->Close();
}

string Get_PDGname(int PDGcode){
    string ret;
    switch (PDGcode) {
        case 1: ret = "down"; break;
        case -1: ret = "anti-down"; break;
        case 2: ret = "up"; break;
        case -2: ret = "anti-up"; break;
        case 3: ret = "strange"; break;
        case -3: ret = "anti-strange"; break;
        case 11: ret = "e-"; break;
        case -11: ret = "e+"; break;
        case 21: ret = "gluon"; break;
        case 22: ret = "gamma"; break;
        case 111: ret = "pi0"; break;
        case 211: ret = "pi+"; break;
        case -211: ret = "pi-"; break;
        case 221: ret = "eta"; break;
        case 113: ret = "rho0"; break;
        case 213: ret = "rho+"; break;
        case -213: ret = "rho-"; break;
        case 130: ret = "KL0"; break;
        case 310: ret = "KS0"; break;
        case 311: ret = "K0"; break;
        case -311: ret = "anti-K0"; break;
        case -313: ret = "anti-K*(892)0"; break;
        case 313: ret = "K*(892)0"; break;
        case 323: ret = "K*(892)+"; break;
        case -323: ret = "K*(892)-"; break;
        case 321: ret = "K+"; break;
        case -321: ret = "K-"; break;
        case 331: ret = "eta'(958)"; break;
        case 333: ret = "phi (1020)"; break;
        case 2101: ret = "(ud)0"; break;
        case 2103: ret = "(ud)1"; break;
        case 2203: ret = "(uu)1"; break;
        case 2112: ret = "n"; break;
        case 2212: ret = "p"; break;
        case 1114: ret = "Delta-"; break;
        case 2114: ret = "Delta0"; break;
        case 2214: ret = "Delta+"; break;
        case 2224: ret = "Delta++"; break;
        case 3114: ret = "Sigma*-"; break;
        case 3122: ret = "Lambda"; break;
        case 3112: ret = "Sigma-"; break;
        case 3212: ret = "Sigma0"; break;
        case 3214: ret = "Sigma*0"; break;
        case 3222: ret = "Sigma+"; break;
        case 3224: ret = "Sigma(1385)"; break;
        default: ret = "unknown"; break;
    }
    return ret;
}

// Check that polar angle is within the geometrical acceptance of the EC.  Input angle in degrees.
bool Cut_CLAS6_Theta_EC(Double_t theta){
    Double_t theta_min = 10.0; // lower limit on polar angle, in degrees
    Double_t theta_max = 60.0; // lower limit on polar angle, in degrees

    return (theta >= theta_min && theta < theta_max) ? true : false;
}

// Check that polar angle is within the geometrical acceptance of the TOF.  Input angle in degrees.
bool Cut_CLAS6_Theta_TOF(Double_t theta){
    Double_t theta_min = 10.0; // lower limit on polar angle, in degrees
    Double_t theta_max = 160.0; // lower limit on polar angle, in degrees
    
    return (theta >= theta_min && theta < theta_max) ? true : false;
}

bool Cut_CLAS6_Theta_Ana(TLorentzVector V4, Int_t iAna, string detName){
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
