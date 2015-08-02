// PlotSCMassSquared.C
//
// macro to analyze TOF mass-squared
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();   // start from scratch

const Int_t NPART = 5;

Int_t lcol[10] = {1,2,4,6,7,8,9,13,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};

Float_t Lmar = 0.125; // set the left margin
Float_t Rmar = 0.125; // set the right margin
Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values

Float_t xLo[5] = {-0.4,-0.4,-0.4,-0.4,-0.4};
Float_t xHi[5] = {0.4,0.4,0.4,0.4,0.4};

char *legHeader[3] = {"Cuts: ","All Cuts Except:","Anti-Cuts:"};
char *HistName[4] = {"scMassSquared","scMassSquared_elecID","scMassSquared_photID","scMassSquared_elecIDphotID"};
char *RunName[4] = {"C12","Fe56","Sn","Pb208"};
char *PartName[5] = {"Electron","#pi^{-}","#pi^{+}","Photon 1","Photon 2"};

// 
// PlotSCMassSquared_Particle - plot histogram with labels
//                  
//                  fAna = output from eg2a DMS
//                  histIndex = histogram index
//                  tgtIndex = target index
//                  chan = particle channel
//
void PlotSCMassSquared_Particle(char *fAna,  Int_t histIndex =0, Int_t tgtIndex = 0, Int_t chan = 0)
{
	char OutCan[100];
    char strname[100];
    
	// Canvas to plot histogram
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
	c1->SetBorderSize(5); 
	gStyle->SetOptStat(0);
	c1->SetFillStyle(4000);
	
	// data files contain the trees
	printf("Analyzing file %s\n",fAna);  
	TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("Kinematics");
	
	c1->cd();
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
    
	TH2D *h2D = (TH2D*)tmp->Get(HistName[histIndex]);
    TH1D *h1D;

    sprintf(strname,"%s_px_%i_%i",HistName[histIndex],chan,chan);
    h1D = (TH1D*)h2D->ProjectionX(strname,chan+1,chan+1,"");
    
	h1D->SetTitle(0);
	h1D->GetXaxis()->CenterTitle();
	h1D->GetYaxis()->CenterTitle();
    h1D->GetYaxis()->SetTitle("Counts");
	h1D->GetYaxis()->SetTitleOffset(yoff);
//    h1D->SetAxisRange(xLo[chan],xHi[chan],"X");
    h1D->SetLineWidth(2);
    h1D->Draw();

	sprintf(OutCan,"Plot_%s_%i_%s.gif",HistName[histIndex],chan,RunName[tgtIndex]);
	c1->Print(OutCan);
	sprintf(OutCan,"Plot_%s_%i_%s.eps",HistName[histIndex],chan,RunName[tgtIndex]);
	c1->Print(OutCan);
}

//
// PlotSCMassSquared_Particle - plot histogram with labels
//
//                  fAna = output from eg2a DMS
//                  histIndex = histogram index
//                  tgtIndex = target index
//                  chan = particle channel
//
void OverlaySCMassSquared_All(char *fAna,  Int_t histIndex =0, Int_t tgtIndex = 0)
{
    Int_t i;
    Int_t ymax = 0;
    Int_t ytmp;
    
    char OutCan[100];
    char strname[100];
    char legLabel[50];
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna);
    TFile *fm = new TFile(fAna,"READ");
    TDirectory *tmp = fm->GetDirectory("Kinematics");
    
    c1->cd();
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TLegend *leg = new TLegend(0.6,0.5,1.0,0.875);
    
    TH2D *h2D = (TH2D*)tmp->Get(HistName[histIndex]);
    TH1D *h1D[NPART];
    
    for(i=0; i<NPART; i++){
        sprintf(strname,"%s_px_%i_%i",HistName[histIndex],i,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
    
        h1D[i]->SetTitle(0);
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetLineWidth(2);
        h1D[i]->SetLineColor(i+1);
        ytmp = h1D[i]->GetMaximum();
        if(ytmp > ymax){
            cout<< ymax <<" "<<ytmp<<endl;
            h1D[0]->SetMaximum(ytmp*1.1);
            ymax = ytmp;
        }
        h1D[i]->Draw(fSame[i]);
        
        sprintf(legLabel,"%s",PartName[i]);
        leg->AddEntry(h1D[i],legLabel,"l");
    }
    
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Particles:");
    leg->Draw();
    
    sprintf(OutCan,"OL_%s_%s.gif",HistName[histIndex],RunName[tgtIndex]);
    c1->Print(OutCan);
    sprintf(OutCan,"OL_%s_%s.eps",HistName[histIndex],RunName[tgtIndex]);
    c1->Print(OutCan);
}



