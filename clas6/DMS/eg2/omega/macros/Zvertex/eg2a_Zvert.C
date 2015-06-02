// eg2a_Zvert.C
//
// macro to plot eg2a vertex histograms
// 
// Michael H. Wood, Canisius College
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();   // start from scratch

const Int_t NSECTORS= 6;

Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};

char hname[50];
char htitle[500];
char title[500];
char cname[50];
char ctitle[500];
char xtitle[100];
char ytitle[100];
char OutCan[100];
char OutText[100];

Float_t Lmar = 0.125;
Float_t Rmar = 0.125;
Float_t yoff = 1.5;

// 
// PlotZvert - plot histogram with labels
//                  
//                  fAna = output from eg2a DMS
//                  suffix = string to append to end of image file name
//                  target = target name
//
void PlotZvert(char *fAna="Ana.root", char *target, int sector)
{
    char strname[50];
    
    if(sector>0 && sector<=6){
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
    
        sprintf(hname,"elecZVertSector");
        TH2F *hZvert = (TH2F*)tmp->Get(hname);

        sprintf(strname,"%s_px_%i",hname,sector);
        hZvert_proj = (TH1D*)hZvert->ProjectionX(strname,sector,sector,"");
    
        sprintf(title,"%s, Sector %i",target,sector);
        hZvert_proj->SetTitle(title);
        hZvert_proj->SetXTitle("e^{-} Z Vertex (cm)");
        hZvert_proj->GetXaxis()->CenterTitle();
        hZvert_proj->SetYTitle("Counts");
        hZvert_proj->GetYaxis()->CenterTitle();
        hZvert_proj->GetYaxis()->SetTitleOffset(yoff);
        hZvert_proj->SetLineWidth(2);
        hZvert_proj->Draw();

        sprintf(OutCan,"PlotZvert_%s_S%i.gif",target,sector);
        c1->Print(OutCan);
        sprintf(OutCan,"PlotZvert_%s_S%i.eps",target,sector);
        c1->Print(OutCan);
    }else{
        cout<<"Wrong sector number "<<sector<<endl;
    }
}

//
// OverlayZvertBySector - overlay histogram of e- Z vertex by sector
//
//                  fAna = output from eg2a DMS
//                  suffix = string to append to end of image file name
//                  target = target name
//
void OverlayZvertBySector(char *fAna="Ana.root", char *target)
{
	Int_t i;
	char legLabel[50];
    char strname[50];
    
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
 
    sprintf(hname,"elecZVertSector");
    TH2F *hZvert = (TH2F*)tmp->Get(hname);
    
	TH1D *hZSec[NSECTORS];
    
	TLegend *leg = new TLegend(0.60,0.50,1.0,0.85); //declare Legend and give its location
    
    	for(i=0; i<NSECTORS ; i++){
            sprintf(strname,"%s_px_%i",hname,i+1);
            hZSec[i] = (TH1D*)hZvert->ProjectionX(strname,i+1,i+1,"");
        	hZSec[i]->SetTitle(0);
        	hZSec[i]->SetXTitle("e^{-} Z Vertex (cm)");
        	hZSec[i]->GetXaxis()->CenterTitle();
        	hZSec[i]->SetYTitle("Counts");
        	hZSec[i]->GetYaxis()->CenterTitle();
        	hZSec[i]->GetYaxis()->SetTitleOffset(yoff);
        	hZSec[i]->SetLineWidth(2);
        	hZSec[i]->SetLineColor(lcol[i]);
        	hZSec[i]->Draw(fSame[i]);
        
        	sprintf(legLabel,"Sector %i",i+1);
        	leg->AddEntry(hZSec[i],legLabel,"l");
    	}
    	leg->SetLineColor(0);
    	leg->SetFillStyle(0);
    	leg->SetHeader(target);
    	leg->Draw();
    
	sprintf(OutCan,"OverlayZvertBySector_%s.gif",target);
	c1->Print(OutCan);
	sprintf(OutCan,"OverlayZvertBySector_%s.eps",target);
	c1->Print(OutCan);
}

