// IMomega_pythiaOmega.C
//
// macro to plot histograms from hIMomega of Ana_pythiaOmega.C
//
//          fAna = output from Ana_pythiaOmega.C
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

void IMomega_pythiaOmega_all(string fAna)
{
    Int_t i;
	char OutCan[100];
    char strname[100];
    
    PartComb myComb;
    const Int_t nHists = myComb.Get_nPartComb();
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,1400,800);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5); 
	gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    c1->Divide(4,4);
    
	// data files contain the trees
	printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
	TH2D *h2D = (TH2D*)fm->Get("hIMomega");
    TH1D *h1D[nHists];
    
    for(i=1; i<=nHists; i++){
        c1->cd(i);
        h2D->Draw();
        sprintf(strname,"hIMomega_px_%i_%i",i,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i,i,"");
        h1D[i]->SetTitle(myComb.Get_PartComb(i-1).c_str());
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(2);
        h1D[i]->SetFillColor(2);
        h1D[i]->Draw();
    }
    
	sprintf(OutCan,"IMomega_pythiaOmega_all.gif");
	c1->Print(OutCan);
	sprintf(OutCan,"IMomega_pythiaOmega_all.eps");
	c1->Print(OutCan);
}

void IMomega_pythiaOmega_single(string fAna, Int_t chan = 1)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    
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
    printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TH2D *h2D = (TH2D*)fm->Get("hIMomega");
    TH1D *h1D;
    
    sprintf(strname,"hIMomega_px_%i_%i",chan,chan);
    h1D = (TH1D*)h2D->ProjectionX(strname,chan,chan,"");
    h1D->SetTitle(myComb.Get_PartComb(chan-1).c_str());
    h1D->GetXaxis()->CenterTitle();
    h1D->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    h1D->GetYaxis()->CenterTitle();
    h1D->GetYaxis()->SetTitle("Counts");
    h1D->GetYaxis()->SetTitleOffset(yoff);
    h1D->SetLineWidth(2);
    h1D->SetLineColor(2);
    h1D->SetFillColor(2);
    h1D->Draw();
    
    sprintf(OutCan,"IMomega_pythiaOmega_Chan%i.gif",chan);
    c1->Print(OutCan);
    sprintf(OutCan,"IMomega_pythiaOmega_Chan%i.eps",chan);
    c1->Print(OutCan);
}

