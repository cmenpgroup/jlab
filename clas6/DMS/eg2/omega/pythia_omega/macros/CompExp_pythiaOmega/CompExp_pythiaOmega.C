// CompExp_pythiaOmega.C
//
// macro to plot histograms comapring pythia sim. with
// experimental data
//
//          fSim = output from Ana_pythiaOmega.C
//          fExp = output from data analysis
//
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//gROOT->Reset();   // start from scratch

Float_t Lmar = 0.125;
Float_t Rmar = 0.125;
Float_t yoff = 1.5;

Int_t lcol[10] = {1,2,4,6,7,8,9,13,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};

class Experiment
{
    vector<string> LabelTgt; // list of target types
    vector<string> LabelRun; // list of runs by solid target
    vector<string> LabelHist; // list of histogram prefixes
    vector<string> LabelHistLeg; // list of histogram legends
    vector<string> LabelCuts; // list of cut names
    
public:
    Experiment();
    Int_t Get_nTgt() {return LabelTgt.size();};
    string Get_Tgt(int num) {return LabelTgt[num];};
    Int_t Get_nRun() {return LabelRun.size();};
    string Get_Run(int num) {return LabelRun[num];};
    Int_t Get_nHist() {return LabelHist.size();};
    string Get_Hist(int num) {return LabelHist[num];};
    Int_t Get_nHistLeg() {return LabelHistLeg.size();};
    string Get_HistLeg(int num) {return LabelHistLeg[num];};
    Int_t Get_nCuts() {return LabelCuts.size();};
    string Get_Cuts(int num) {return LabelCuts[num];};
    
    void Check_HistIndex(int index);
    void PrintHistIndex();
    void Check_HistLegIndex(int index);
    void PrintHistLegIndex();
    void Check_TgtIndex(int index);
    void PrintTgtIndex();
    void Check_CutIndex(int index);
    void PrintCutIndex();
};

Experiment::Experiment()
{
    LabelTgt.push_back("NoTarget");
    LabelTgt.push_back("LD2");
    LabelTgt.push_back("Nuc");

    LabelRun.push_back("C12");
    LabelRun.push_back("Fe56");
    LabelRun.push_back("Sn");
    LabelRun.push_back("Pb208");

    LabelHist.push_back("IMOmega_");
    LabelHist.push_back("IMOmega_woCut_");
    LabelHist.push_back("IMOmega_antiCut_");

    LabelHistLeg.push_back("Cuts: ");
    LabelHistLeg.push_back("All Cuts Except:");
    LabelHistLeg.push_back("Anti-Cuts:");

    LabelCuts.push_back("None");
    LabelCuts.push_back("All #omega cuts");
    LabelCuts.push_back("M(#pi^{0})");
    LabelCuts.push_back("Q^{2}");
    LabelCuts.push_back("W");
    LabelCuts.push_back("V_{z} Matching");
    LabelCuts.push_back("Part. Topology");
    LabelCuts.push_back("#theta_{e-,#gamma}");
    LabelCuts.push_back("M(#pi^{+}#pi^{-}");
    LabelCuts.push_back("EC 2nd Moment for #gamma's, Region 1");
    LabelCuts.push_back("EC 2nd Moment for #gamma's, Region 2");
    LabelCuts.push_back("EC 2nd Moment for #gamma's, Region 3");
    LabelCuts.push_back("EC 3rd Moment for #gamma's, Region 1");
    LabelCuts.push_back("EC 3rd Moment for #gamma's, Region 2");
    LabelCuts.push_back("EC 3rd Moment for #gamma's, Region 3");
    LabelCuts.push_back("Photon TOF M^{2}");
    
}

void Experiment::Check_HistIndex(Int_t index)
{
    int MAX_HIST = this->Get_nHist();
    if(index<0 || index>=MAX_HIST){
        cout<<"Histogram index "<<index<<" is out of range [0,"<<MAX_HIST-1<<"]"<<endl;
        this->PrintHistIndex();
        exit(0);
    }
}

void Experiment::PrintHistIndex()
{
    Int_t i;
    cout<<"Histogram Index:"<<endl;
    for(i=0;i<this->Get_nHist();i++){
        cout<<i<<"\t"<<this->Get_Hist(i)<<endl;
    }
}

void Experiment::Check_HistLegIndex(Int_t index)
{
    int MAX_HLEG = this->Get_nHistLeg();
    if(index<0 || index>=MAX_HLEG){
        cout<<"Histogram Legend index "<<index<<" is out of range [0,"<<MAX_HLEG-1<<"]"<<endl;
        this->PrintHistLegIndex();
        exit(0);
    }
}

void Experiment::PrintHistLegIndex()
{
    Int_t i;
    cout<<"Histogram Legend Index:"<<endl;
    for(i=0;i<this->Get_nHistLeg();i++){
        cout<<i<<"\t"<<this->Get_HistLeg(i)<<endl;
    }
}

void Experiment::Check_TgtIndex(Int_t index)
{
    int MAX_TGT = this->Get_nTgt();
    if(index<0 || index>=MAX_TGT){
        cout<<"Target index "<<index<<" is out of range [0,"<<MAX_TGT-1<<"]"<<endl;
        this->PrintTgtIndex();
        exit(0);
    }
}

void Experiment::PrintTgtIndex()
{
    Int_t i;
    cout<<"Target Index:"<<endl;
    for(i=0;i<this->Get_nTgt();i++){
        cout<<i<<"\t"<<this->Get_Tgt(i)<<endl;
    }
}

void Experiment::Check_CutIndex(Int_t index)
{
    int MAX_CUT = this->Get_nCuts();
    if(index<0 || index>=MAX_CUT){
        cout<<"Cut index "<<index<<" is out of range [0,"<<MAX_CUT-1<<"]"<<endl;
        this->PrintCutIndex();
        exit(0);
    }
}

void Experiment::PrintCutIndex()
{
    Int_t i;
    cout<<"Cut Index:"<<endl;
    for(i=0;i<this->Get_nCuts();i++){
        cout<<i<<"\t"<<this->Get_Cuts(i)<<endl;
    }
}

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

void CompExp_pythiaOmega(string fExp, string fSim, Int_t chan = 1, Int_t tgtIndex = 0)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    char hname[50];
    char title[100];
    char legLabel[50];
    
    string fSame;
    
    Experiment myExp;
    Int_t histIndex = 0;
    Int_t cutIndex = 1;
    
    myExp.Check_HistIndex(histIndex);
    myExp.Check_TgtIndex(tgtIndex);
    myExp.Check_CutIndex(cutIndex);
    
    PartComb myComb;
    Int_t nHists = myComb.Get_nPartComb();
    
    if(chan<=0 || chan>nHists){
        cout<<"Incorrect channel "<<chan<<endl;
        exit(0);
    }
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    // data files contain the trees

    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TFile *fm[2];
    TDirectory *dirExp;
    TH2D *h2D[2];
    TH1D *h1D[2];
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.70,0.85,0.875);
    
    for(i=0; i<2; i++){
        // data files contain the trees
        switch(i){
            case 1:
                sprintf(legLabel,"%s","Pythia");
                printf("Analyzing file %s\n",fSim.c_str());
                fm[i] = new TFile(fSim.c_str(),"READ");
                h2D[i] = (TH2D*)fm[i]->Get("hIMomega");
                sprintf(strname,"hIMomega_px_%i_%i",chan,chan);
                h1D[i] = (TH1D*)h2D[i]->ProjectionX(strname,chan,chan,"");
                fSame = "same";
                break;
            case 0:
                sprintf(legLabel,"%s","eg2 Data");
                printf("Analyzing file %s\n",fExp.c_str());
                fm[i] = new TFile(fExp.c_str(),"READ");
                dirExp = fm[i]->GetDirectory(myExp.Get_Tgt(tgtIndex).c_str());
                
                sprintf(hname,"%s%s",myExp.Get_Hist(histIndex).c_str(),myExp.Get_Tgt(tgtIndex).c_str());
                h2D[i] = (TH2D*)dirExp->Get(hname);

                sprintf(strname,"%s_%i",hname,i);
                h1D[i] = (TH1D*)h2D[i]->ProjectionX(strname,cutIndex+1,cutIndex+1,"");
                fSame = "";
                break;
            default: sprintf(legLabel,"Insert File Name"); break;
        }
    
        sprintf(title,"Target: %s, Cuts: %s",myExp.Get_Tgt(tgtIndex).c_str(),myExp.Get_Cuts(cutIndex).c_str());
        h1D[i]->SetTitle(title);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(lcol[i]);
        h1D[i]->Draw(fSame.c_str());
        
        leg->AddEntry(h1D[i],legLabel,"l");
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Files:");
    leg->Draw();
    
    sprintf(OutCan,"CompExp_pythiaOmega_Chan%i.gif",chan);
    c1->Print(OutCan);
    sprintf(OutCan,"CompExp_pythiaOmega_Chan%i.eps",chan);
    c1->Print(OutCan);
}
