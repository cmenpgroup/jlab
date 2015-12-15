// PlotHists_pythiaOmega.C
//
// macro to plot eg2a histograms
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
Float_t yoff = 1.75;

// 
// PlotHists_pythiaOmega - plot histogram with labels
//                  
//                  fAna = output from eg2a DMS
//                  hname = histogram name
//
void PlotHists_pythiaOmega(char *fAna, char *hname, int iDim = 1)
{
	char OutCan[100];
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5); 
	gStyle->SetOptStat(0);
	c1->SetFillStyle(4000);
	
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);  
	TFile *fm = new TFile(fAna,"READ");
	
	c1->cd();
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
	TH1F *hist = (TH1F*)fm->Get(hname);
	hist->SetTitle();
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->CenterTitle();
	hist->GetYaxis()->SetTitleOffset(yoff);
    switch(iDim){
        case 1: hist->SetLineWidth(2); hist->Draw(); break;
        case 2: hist->Draw(); break;
//        case 2: hist->Draw("colz"); break;
        default:
            cout << "Incorrect dimension for histogram." << endl;
            exit(0);
            break;
    }

	sprintf(OutCan,"pythiaOmega_%s.gif",hname);
	c1->Print(OutCan);
	sprintf(OutCan,"pythiaOmega_%s.eps",hname);
	c1->Print(OutCan);
}

