// PartComb_pythiaOmega.C
//
// macro to plot histograms from hPartPerEvt of Ana_pythiaOmega.C
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

const Int_t NUM_HISTS = 6;

char* title[NUM_HISTS] = {"#pi^{+}","#pi^{-}","#pi^{0}","photons","#omega","Combinations"};

void PartComb_pythiaOmega(char *fAna)
{
    Int_t i;
	char OutCan[100];
    char strname[100];
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,800);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5); 
	gStyle->SetOptStat(0);
	c1->SetFillStyle(4000);
    c1->Divide(3,2);
    
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);  
	TFile *fm = new TFile(fAna,"READ");
	
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
	TH2D *h2D = (TH2D*)fm->Get("hPartPerEvt");
    TH1D *h1D[NUM_HISTS];
    
    for(i=1; i<=NUM_HISTS; i++){
        c1->cd(i);
        sprintf(strname,"hPartPerEvt_px_%i_%i",i,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i,i,"");
        h1D[i]->SetTitle(title[i-1]);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetXaxis()->SetTitle("Particles per event");
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetFillColor(2);
        h1D[i]->Draw();
    }
    
	sprintf(OutCan,"PartComb_pythiaOmega.gif");
	c1->Print(OutCan);
	sprintf(OutCan,"PartComb_pythiaOmega.eps");
	c1->Print(OutCan);
}

