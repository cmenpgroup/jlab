// MMsq_pythiaOmega.C
//
// macro to plot histograms from hMMsq_proton of Ana_pythiaOmega.C
//
//          fAna = output from Ana_pythiaOmega.C
//
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();   // start from scratch

Float_t Lmar = 0.125;
Float_t Rmar = 0.125;
Float_t yoff = 1.5;

const Int_t NUM_HISTS = 12;

char* title[NUM_HISTS] = {
    "N(#pi^{+})>=1 && N(#pi^{-})>=1 && N(#pi^{0})>=1",
    "N(#pi^{+})==1 && N(#pi^{-})==1 && N(#pi^{0})==1",
    "N(#pi^{+})>=1 && N(#pi^{-})>=1 && N(#gamma)>=2",
    "N(#pi^{+})==1 && N(#pi^{-})==1 && N(#gamma)==2",
    "N(#pi^{+})>=1 && N(#pi^{-})>=1 && N(#pi^{0})>=1 && #omega",
    "N(#pi^{+})==1 && N(#pi^{-})==1 && N(#pi^{0})==1 && #omega",
    "N(#pi^{+})>=1 && N(#pi^{-})>=1 && N(#gamma)>=2 && #omega",
    "N(#pi^{+})==1 && N(#pi^{-})==1 && N(#gamma)==2 && #omega",
    "N(#pi^{+})>=1 && N(#pi^{-})>=1 && N(#pi^{0})>=1 && !#omega",
    "N(#pi^{+})==1 && N(#pi^{-})==1 && N(#pi^{0})==1 && !#omega",
    "N(#pi^{+})>=1 && N(#pi^{-})>=1 && N(#gamma)>=2 && !#omega",
    "N(#pi^{+})==1 && N(#pi^{-})==1 && N(#gamma)==2 && !#omega",
};

void MMsq_pythiaOmega_all(char *fAna)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,1400,800);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5); 
	gStyle->SetOptStat(0);
	c1->SetFillStyle(4000);
    c1->Divide(4,3);
    
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);  
	TFile *fm = new TFile(fAna,"READ");
	
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
	TH2D *h2D = (TH2D*)fm->Get("hMMsq_proton");
    TH1D *h1D[NUM_HISTS];
    
    for(i=1; i<=NUM_HISTS; i++){
        c1->cd(i);
        sprintf(strname,"hMMsq_px_%i_%i",i,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i,i,"");
        h1D[i]->SetTitle(title[i-1]);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetXaxis()->SetTitle("Particles per event");
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(2);
        h1D[i]->SetFillColor(2);
        h1D[i]->Draw();
    }
    
	sprintf(OutCan,"MMsq_pythiaOmega_all.gif");
	c1->Print(OutCan);
	sprintf(OutCan,"MMsq_pythiaOmega_all.eps");
	c1->Print(OutCan);
}

void MMsq_pythiaOmega_single(char *fAna, Int_t chan = 1)
{
    Int_t i;
    char Prefix[100];
    char OutCan[100];
    char strname[100];
    
    if(chan<=0 || chan>NUM_HISTS){
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
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TH2D *h2D = (TH2D*)fm->Get("hMMsq_proton");
    TH1D *h1D;
    
    sprintf(strname,"hMMsq_px_%i_%i",chan,chan);
    h1D = (TH1D*)h2D->ProjectionX(strname,chan,chan,"");
    h1D->SetTitle(title[chan-1]);
    h1D->GetXaxis()->CenterTitle();
    h1D->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    h1D->GetYaxis()->CenterTitle();
    h1D->GetYaxis()->SetTitle("Counts");
    h1D->GetYaxis()->SetTitleOffset(yoff);
    h1D->SetLineWidth(2);
    h1D->SetLineColor(2);
    h1D->SetFillColor(2);
    h1D->Draw();
    
    if(chan<10){
        sprintf(Prefix,"MMsq_pythiaOmega_Chan0%i",chan);
    }else{
        sprintf(Prefix,"MMsq_pythiaOmega_Chan%i",chan);
    }
    
    sprintf(OutCan,"%s.gif",Prefix);
    c1->Print(OutCan);
    sprintf(OutCan,"%s.eps",Prefix);
    c1->Print(OutCan);
}

