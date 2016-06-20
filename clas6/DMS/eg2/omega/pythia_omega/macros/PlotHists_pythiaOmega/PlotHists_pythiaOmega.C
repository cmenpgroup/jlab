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

void Project1D_pythiaOmega_all(string fAna, string HistName, string xtitle)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    
    PartComb myComb;
    const Int_t nHists = myComb.Get_nPartComb();
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,680);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    c1->Divide(sqrt(nHists),sqrt(nHists));
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TH2D *h2D = (TH2D*)fm->Get(HistName.c_str());
    TH1D *h1D[nHists];
    
    for(i=0; i<nHists; i++){
        c1->cd(i+1);
        sprintf(strname,"hist_px_%i_%i",i,i);
        h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
        h1D[i]->SetTitle(myComb.Get_PartComb(i).c_str());
        h1D[i]->GetXaxis()->SetTitle(xtitle.c_str());
        h1D[i]->GetXaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitle("Counts");
        h1D[i]->GetYaxis()->CenterTitle();
        h1D[i]->GetYaxis()->SetTitleOffset(yoff);
        h1D[i]->SetFillColor(2);
        h1D[i]->Draw();
    }
    
    sprintf(OutCan,"pythiaOmega_%s_all.gif", HistName.c_str());
    c1->Print(OutCan);
    sprintf(OutCan,"pythiaOmega_%s_all.eps", HistName.c_str());
    c1->Print(OutCan);
}

void Plot2D_pythiaOmega_all(string fAna, string HistName)
{
    Int_t i;
    char OutCan[100];
    char strname[100];
    
    PartComb myComb;
    const Int_t nHists = myComb.Get_nPartComb();
    
    // Canvas to plot histogram
    TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,680);
    c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    c1->SetBorderSize(5);
    gStyle->SetOptStat(0);
    c1->SetFillStyle(4000);
    c1->Divide(sqrt(nHists),sqrt(nHists));
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    TH2D *h2D[nHists];
    
    for(i=0; i<nHists; i++){
        c1->cd(i+1);
        sprintf(strname,"%s_%i",HistName.c_str(),i);
        h2D[i] = (TH2D*)fm->Get(strname);
        h2D[i]->SetTitle(myComb.Get_PartComb(i).c_str());
        h2D[i]->GetXaxis()->CenterTitle();
        h2D[i]->GetYaxis()->CenterTitle();
        h2D[i]->GetYaxis()->SetTitleOffset(yoff);
        h2D[i]->SetFillColor(2);
        h2D[i]->Draw("colz");
    }
    
    sprintf(OutCan,"pythiaOmega_%s_all.gif", HistName.c_str());
    c1->Print(OutCan);
    sprintf(OutCan,"pythiaOmega_%s_all.eps", HistName.c_str());
    c1->Print(OutCan);
}

void Plot2D_pythiaOmega_sets(string fAna, string HistName)
{
    Int_t i, j;
    char ctitle[100];
    char OutCan[100];
    char strname[100];
    
    PartComb myComb;
    const Int_t nHists = myComb.Get_nPartComb();
    
    // data files contain the trees
    printf("Analyzing file %s\n",fAna.c_str());
    TFile *fm = new TFile(fAna.c_str(),"READ");
    
//    const Int_t nSets = ceil(nHists/4);
    const Int_t nSets = 4;

    cout<<"Sets "<<nSets<<endl;
    TCanvas *can[nSets];
    TH2D *h2D[nHists];
    
    Int_t hctr = 0; // histogram counter
    
    for(j=0; j<nSets; j++){

        sprintf(ctitle,"can%i",j);
        can[j] = new TCanvas(ctitle,ctitle,j*50,0,700,700);
        can[j]->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
        can[j]->SetBorderSize(5);
        gStyle->SetOptStat(0);
        can[j]->SetFillStyle(4000);
        can[j]->Divide(2,2);
        
        gPad->SetLeftMargin(Lmar);
        gPad->SetRightMargin(Rmar);
        gPad->SetFillColor(0);
        
        for(i=0; i<4; i++){
            can[j]->cd(i+1);
            sprintf(strname,"%s_%i",HistName.c_str(),hctr);
            h2D[hctr] = (TH2D*)fm->Get(strname);
            h2D[hctr]->SetTitle(myComb.Get_PartComb(hctr).c_str());
            h2D[hctr]->GetXaxis()->CenterTitle();
            h2D[hctr]->GetYaxis()->CenterTitle();
            h2D[hctr]->GetYaxis()->SetTitleOffset(1.0);
            h2D[hctr]->SetFillColor(2);
            h2D[hctr]->Draw("colz");

            hctr++;
        }
    
        sprintf(OutCan,"pythiaOmega_%s_set%i.gif", HistName.c_str(),j);
        can[j]->Print(OutCan);
        sprintf(OutCan,"pythiaOmega_%s_set%i.eps", HistName.c_str(),j);
        can[j]->Print(OutCan);
    }
}