#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

#define LIMIT_QSQ 1.0
#define LIMIT_W 2.0

#define nCtr 6
// 1 : Number of pi+
// 2 : Number of pi-
// 3 : Number of pi0
// 4 : Number of 2 photons
// 5 : Number of omega mesons
// 6 : Number of particle combinations listed below

#define nComb 12
// Particle combinations:
// 1 : nPip>=1 && nPim>=1 && nPi0>=1
// 2 : nPip==1 && nPim==1 && nPi0==1
// 3 : nPip>=1 && nPim>=1 && nPhoton>=2
// 4 : nPip==1 && nPim==1 && nPhoton==2
// 5 : nPip>=1 && nPim>=1 && nPi0>=1 && nOmega>=1
// 6 : nPip==1 && nPim==1 && nPi0==1 && nOmega==1
// 7 : nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega>=1
// 8 : nPip==1 && nPim==1 && nPhoton==2 && nOmega==1
// 9 : nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0
// 10 : nPip==1 && nPim==1 && nPi0==1 && nOmega==0
// 11 : nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0
// 12 : nPip==1 && nPim==1 && nPhoton==2 && nOmega==0

#define nCombPhoton 6 // number of particle combination with 2 photons

gROOT->Reset();   // start from scratch

TH1D *hIntermedPart;
TH2D *hPartPerEvt;
TH2D *hIMomega;
TH2D *hMMsq_proton;
TH2D *hOpAng_PipPim;
TH2D *hOpAng_PipPimPi0;
TH1D *hQsq;
TH1D *hNu;
TH1D *hW;
TH2D *hz_fracEnergy;
TH2D *hOpAng_TwoPhoton;
TH2D *hIMomega_VS_IMPipPim[nComb];


char hname[100];
char htitle[100];

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

void Ana_pythiaOmega(char* inFile, char* outFile="Ana_pythiaOmega.root", Int_t MaxEvt=1000000, Int_t ModEvt=100){

    Int_t j;
    Int_t ii, jj;
    Int_t iPC; // particle combination index
    Int_t iPhoton;
    Int_t Npart_prev = 1;
    Int_t nEvt = 0;
    
    Int_t nPip, nPim, nPi0, nPhoton, nOmega; // particle counters per event
    
    Int_t ID_ELECTRON = 11; // PDG electron id
    Int_t ID_PHOTON = 22;  // PDG photon id
    Int_t ID_PION_POS = 211;  // PDG pi+ id
    Int_t ID_PION_NEG = -211;  // PDG pi- id
    Int_t ID_PION_ZERO = 111;  // PDG pi0 id
    Int_t ID_ETA_MESON = 221;  // PDG omega id
    Int_t ID_OMEGA_MESON = 223;  // PDG omega id
    Int_t ID_PHI_MESON = 333;  // PDG omega id
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
    
    Float_t Qsq, nu, W;
    Float_t ntNpart, ntks, ntPID, ntParent, ntPx, ntPy, ntPz, ntE, ntM;

    Double_t ChargedPionAngle;
    Double_t ThreePionAngle;
    Double_t TwoPhotonAngle;
    
    vector<int> origNum;
    vector<int> tempNpart;
    vector<int> tempPid;
    vector<int> tempParent;
    vector<int> tempKs;

    bool cuts_Qsq;
    bool cuts_W;
    bool cuts_beam = false;
    bool cuts_target = false;
    bool cuts_electron = false;
    bool cuts_recoil = false;
    
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
    printf("Analyzing file %s\n",inFile);
    TFile *fm = new TFile(inFile,"READ");
    
    TNtuple *myNtuple = (TNtuple*)fm->Get("b");
    
    Int_t nEntries = (Int_t)myNtuple->GetEntries();
    
    cout<<"Number of entries = "<<nEntries<<endl;
    
    if(MaxEvt<0) MaxEvt = nEntries;
    
    myNtuple->SetBranchAddress("i",&ntNpart);
    myNtuple->SetBranchAddress("ks",&ntks);
    myNtuple->SetBranchAddress("kf",&ntPID);
    myNtuple->SetBranchAddress("orig",&ntParent);
    myNtuple->SetBranchAddress("p_x",&ntPx);
    myNtuple->SetBranchAddress("p_y",&ntPy);
    myNtuple->SetBranchAddress("p_z",&ntPz);
    myNtuple->SetBranchAddress("e",&ntE);
    myNtuple->SetBranchAddress("m",&ntM);
    
    for(ii=0; ii<nEntries; ii++){

        myNtuple->GetEntry(ii); // retrieve the event from the ntuple

        if(ntNpart<Npart_prev){
            if(cuts_beam && cuts_target && cuts_electron){
                cuts_Qsq = false;
                cuts_W =false;
                
                BeamMinusElectron = beam - electron;
                
                Qsq = -1.0*BeamMinusElectron.M();
                hQsq->Fill(Qsq);
                
                nu = BeamMinusElectron.E();
                hNu->Fill(nu);
                
                W_TLV = BeamMinusElectron + target;
                W = W_TLV.M();
                hW->Fill(W);
                
                cuts_Qsq = (Qsq >= LIMIT_QSQ);
                cuts_W = (W >= LIMIT_W);
                
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
                
                if(cuts_Qsq && cuts_W){
                    if(cuts_recoil){
                        missing = beam + target - electron - recoil;
                    }
                    
                    if(nPip>=1 && nPim>=1 && nPi0>=1){
                        PipPim = pip + pim;
                        PipPimPi0 = pip + pim + pi0;

                        FillHists_omega(pip,pim,pi0,nu,1);
                        if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),1);

                        if(nOmega){
                            FillHists_omega(pip,pim,pi0,nu,5);
                            if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),5);
                        }else{
                            FillHists_omega(pip,pim,pi0,nu,9);
                            fout09<<endl<<"Event: "<<nEvt<<endl;
                            for(jj=0; jj<Npart_prev; jj++){
                                fout09<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                            }
                            if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),9);
                        }
                    }

                    if(nPip==1 && nPim==1 && nPi0==1){
                        PipPim = pip + pim;
                        PipPimPi0 = pip + pim + pi0;
                    
                        FillHists_omega(pip,pim,pi0,nu,2);
                        if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),2);
                        
                        if(nOmega){
                            FillHists_omega(pip,pim,pi0,nu,6);
                            if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),6);
                        }else{
                            FillHists_omega(pip,pim,pi0,nu,10);
                            fout10<<endl<<"Event: "<<nEvt<<endl;
                            for(jj=0; jj<Npart_prev; jj++){
                                fout10<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                            }
                            if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),10);
                        }
                    }
                
                    if(nPip>=1 && nPim>=1 && nPhoton>=2){
                        PipPim = pip + pim;
                        twoPhotons = photon[0] + photon[1];
                        PipPim2Photons = pip + pim + twoPhotons;
                    
                        TwoPhotonAngle = TMath::RadToDeg()*photon[0].Angle(photon[1].Vect());
                    
                        FillHists_omega(pip,pim,twoPhotons,nu,3);
                        hOpAng_TwoPhoton->Fill(TwoPhotonAngle,0);
                        if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),3);

                        if(nOmega){
                            FillHists_omega(pip,pim,twoPhotons,nu,7);
                            hOpAng_TwoPhoton->Fill(TwoPhotonAngle,2);
                            if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),7);
                        }else{
                            FillHists_omega(pip,pim,twoPhotons,nu,11);
                            hOpAng_TwoPhoton->Fill(TwoPhotonAngle,4);
                            fout11<<endl<<"Event: "<<nEvt<<endl;
                            for(jj=0; jj<Npart_prev; jj++){
                                fout11<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                            }
                            if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),11);
                        }
                    }
                
                    if(nPip==1 && nPim==1 && nPhoton==2){
                        PipPim = pip + pim;
                        twoPhotons = photon[0] + photon[1];
                        PipPim2Photons = pip + pim + twoPhotons;
                    
                        TwoPhotonAngle = TMath::RadToDeg()*photon[0].Angle(photon[1].Vect());
                    
                        FillHists_omega(pip,pim,twoPhotons,nu,4);
                        hOpAng_TwoPhoton->Fill(TwoPhotonAngle,1);
                        if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),4);

                        if(nOmega){
                            FillHists_omega(pip,pim,twoPhotons,nu,8);
                            hOpAng_TwoPhoton->Fill(TwoPhotonAngle,3);
                            if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),8);
                        }else{
                            FillHists_omega(pip,pim,twoPhotons,nu,12);
                            hOpAng_TwoPhoton->Fill(TwoPhotonAngle,5);
                            fout12<<endl<<"Event: "<<nEvt<<endl;
                            for(jj=0; jj<Npart_prev; jj++){
                                fout12<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<Get_PDGname(tempPid[jj])<<endl;
                            }
                            if(cuts_recoil) hMMsq_proton->Fill(missing.Mag2(),12);
                        }
                    }
                }
            
                for (jj=0; jj<origNum.size(); jj++){
                    hIntermedPart->Fill(origNum[jj]);
                }
            }
            
            nEvt++;
            if(!(nEvt % ModEvt)) cout<<nEvt<<endl;
            
            cuts_beam = false;
            cuts_target = false;
            cuts_electron = false;
            cuts_recoil = false;
            iPhoton = 0;
            nPip = 0;
            nPim = 0;
            nPi0 = 0;
            nPhoton = 0;
            nOmega = 0;
            origNum.clear();
            tempNpart.clear();
            tempKs.clear();
            tempPid.clear();
            tempParent.clear();
        }
        
//        if(nEvt>=0){
//            cout<<nEvt<<"\t"<<ii<<"\t"<<ntNpart<<"\t"<<ntks<<"\t"<<ntPID<<"\t"<<ntParent<<"\t";
//            cout<<ntPx<<"\t"<<ntPy<<"\t"<<ntPz<<"\t"<<ntE<<"\t"<<ntM<<endl;
//        }
        
        tempNpart.push_back(ntNpart);
        tempKs.push_back(ntks);
        tempPid.push_back(ntPID);
        tempParent.push_back(ntParent);
        
        if(ntPID==ID_ELECTRON && ntNpart==1){
            beam.SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
            cuts_beam = true;
        }
        
        if(ntPID==ID_PROTON && ntNpart==2){
            target.SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
            cuts_target =true;
        }
        
        if(ntPID==ID_ELECTRON && ntks==1){
            electron.SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
            cuts_electron = true;
        }
        
        if(ntPID==ID_PROTON && ntks==1){
            recoil.SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
            cuts_recoil = true;
        }
        
        if(ntPID==ID_PION_POS && ntks==1){
            nPip++;
            pip.SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
        }
        
        if(ntPID==ID_PION_NEG && ntks==1){
            nPim++;
            pim.SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
        }
        if(ntPID==ID_PION_ZERO && ntks==11){
            nPi0++;
            pi0.SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
        }
        if(ntPID==ID_PHOTON && ntks==1){
            nPhoton++;
            if(iPhoton<2){
                photon[iPhoton].SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
            }
            iPhoton++;
        }
        if(ntPID==ID_OMEGA_MESON && ntks==11){
            nOmega++;
            omega.SetPxPyPzE(ntPx,ntPy,ntPz,ntE);
        }
        
        if(ntks==11){
            origNum.push_back(ntPID);
        }
        
        Npart_prev = ntNpart;

        if(nEvt==MaxEvt) break;
    }
    cout << "Number of nEvt = " << nEvt << endl;

    WriteHists(outFile); // write the histograms to file
    
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

    Int_t nIMomega = 200;
    Double_t IMomegaLo = 0.0;
    Double_t IMomegaHi = 2.5;

    Int_t nMMsq = 200;
    Double_t MMsqLo = 0.0;
    Double_t MMsqHi = 5.0;
    
    char *LabelPartComb[nComb] ={
    "nPip>=1 && nPim>=1 && nPi0>=1",
    "nPip==1 && nPim==1 && nPi0==1",
    "nPip>=1 && nPim>=1 && nPhoton>=2",
    "nPip==1 && nPim==1 && nPhoton==2",
    "nPip>=1 && nPim>=1 && nPi0>=1 && nOmega>=1",
    "nPip==1 && nPim==1 && nPi0==1 && nOmega==1",
    "nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega>=1",
    "nPip==1 && nPim==1 && nPhoton==2 && nOmega==1",
    "nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0",
    "nPip==1 && nPim==1 && nPi0==1 && nOmega==0",
    "nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0",
    "nPip==1 && nPim==1 && nPhoton==2 && nOmega==0"};

    hIntermedPart = new TH1D("hIntermedPart","Intermediate Particle ID", 3500, 0.5, 3500.5);
    hPartPerEvt = new TH2D("hPartPerEvt", "Part. Comb. per Event",nComb+1,-0.5,nComb+0.5,nCtr,0.5,nCtr+0.5);

    hIMomega = new TH2D("hIMomega", "Mass(#omega)",nIMomega,IMomegaLo,IMomegaHi,nComb,0.5,nComb+0.5);

    hMMsq_proton = new TH2D("hMMsq_proton", "MM^{2} off proton",nMMsq,MMsqLo,MMsqHi,nComb,0.5,nComb+0.5);
    
    hOpAng_PipPim = new TH2D("hOpAng_PipPim","Opening Angle(#pi^{+},#pi^{-})",nTheta,ThetaLo,ThetaHi,nComb,0.5,nComb+0.5);

    hOpAng_PipPimPi0 = new TH2D("hOpAng_PipPimPi0","Opening Angle(#pi^{+}#pi^{-},#pi^{0})",nTheta,ThetaLo,ThetaHi,nComb,0.5,nComb+0.5);

    hQsq = new TH1D("hQsq","Q^{2}", 100, 0., 4.);
    hNu = new TH1D("hNu","#nu", 100, 0., 5.);
    hW = new TH1D("hW", "W", 250, 0., 5.);
    hz_fracEnergy = new TH2D("hz_fracEnergy", "Fractional Energy", 150, 0, 1.5,nComb,0.5,nComb+0.5);
    hOpAng_TwoPhoton = new TH2D("hOpAng_TwoPhoton","Opening Angle(#gamma,#gamma)",nTheta, ThetaLo, ThetaHi,nCombPhoton,0.5,nCombPhoton+0.5);
    
    for(j=0; j<nComb; j++){
        sprintf(hname,"hIMomega_VS_IMPipPim_%i",j);
        sprintf(htitle,"Mass(#omega) vs. Mass(#pi^{+} #pi^{-}), %s",LabelPartComb[j]);
        hIMomega_VS_IMPipPim[j] = new TH2D(hname,htitle,nIM2pi,IM2piLo,IM2piHi,nIMomega,IMomegaLo,IMomegaHi);
    }
}

void WriteHists(char *outFile){

    Int_t j;
    
    TFile *out = new TFile(outFile, "recreate");
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

    hMMsq_proton->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_proton->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_proton->Write();
    
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

void PrintTLorentzVector(TLorentzVector TLV){
    cout <<"Px "<<TLV.Px()<<"\t";
    cout <<"Py "<<TLV.Py()<<"\t";
    cout <<"Pz "<<TLV.Pz()<<"\t";
    cout <<"E " <<TLV.E() <<"\t";
    cout <<"M " <<TLV.M() <<"\t";
    cout <<"M^2 " <<TLV.M2() <<endl;
}
