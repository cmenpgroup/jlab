#include <vector>
#include <string>
#include "HistManager.h"
#include <iostream>

HistManager::HistManager()
{
    MAX_SECTORS = 6;

    ID_ELECTRON = 11;
    ID_PHOTON = 22;
    ID_PION_POS = 211;
    ID_PION_NEG = -211;
    ID_PROTON = 2212;

    MASS_PHOTON = 0.0;
    MASS_ELECTRON = 0.000511;
    MASS_PION_CHARGED = 0.138;
    MASS_PION_NEUTRAL = 0.135;
    MASS_PROTON = 0.938;
    MASS_DEUTERIUM = 2*MASS_PROTON;

    LIGHTSPEED = 30.0;

    BEAM_ENERGY = 5.01; // 4.5
}

void HistManager::SetMass2Pions_VS_massOmega_NCX(int index)
{ 
    mass2Pions_VS_massOmega_NCX[index] = mass2Pions_VS_massOmega_NC[index]->ProjectionX();
}
void HistManager::SetMass2Pions_VS_massOmega_EPCX(int index)
{
    mass2Pions_VS_massOmega_EPCX[index] = mass2Pions_VS_massOmega_EPC[index]->ProjectionX();
}
void HistManager::SetMass2Pions_VS_massOmega_EPOCX(int index)
{
    mass2Pions_VS_massOmega_EPOCX[index] = mass2Pions_VS_massOmega_EPOC[index]->ProjectionX();
}

//
// BookHist - routine to set up histograms
//
void HistManager::BookHist()
{

    int i, j;
    
    char hname[100];
	char htitle[100];
    
    int nIMomega = 125;
    double IMomegaLo = 0.0;
    double IMomegaHi = 2.5;

    int nPtSq_omega = 100;
    double PtSq_omegaLo = 0.0;
    double PtSq_omegaHi = 2.0;
    
    DetectedParticles myDetPart;
    ParticleList myPartList;
    EG2Target myTgt;
    OmegaMixedEvent myMixEvt;
    ElectronID myElecID;
    
    int nME_Methods = myMixEvt.Get_nLabel();
    int nElecID = myElecID.Get_nElecID();

    sprintf(hname,"StartTime");
    sprintf(htitle,"Event Start Time");
    StartTime = new TH1D(hname,htitle, 420, -40., 100.);
    
    sprintf(hname,"q2");
    sprintf(htitle,"Q^{2}");
    q2 = new TH1D(hname,htitle, 100, 0., 4.);

    sprintf(hname,"q2_VS_theta");
    sprintf(htitle,"Q^{2} vs. 4E_{e'}sin^{2}(0.5*#theta_{e'})");
    q2_VS_theta = new TH2D(hname,htitle, 200, 0., 1.0, 200, 0., 4.);
    
    sprintf(hname,"nu_EnergyTransfer");
    sprintf(htitle,"\nu");
    nu_EnergyTransfer = new TH1D(hname,htitle, 100, 0., 5.);

    sprintf(hname,"elecZVert");
    sprintf(htitle,"Z Vertex of Electron");
	elecZVert = new TH1D(hname,htitle, 300, -35, -20);

    sprintf(hname,"elecZVertSector");
    sprintf(htitle,"Z Vertex of Electron vs Sector");
    elecZVertSector = new TH2D(hname, htitle, 300, -40, -10, 6, 0.5, 6.5);

    sprintf(hname,"elecZVertSector_Corr");
    sprintf(htitle,"Z Vertex of Electron vs Sector, Corrected");
    elecZVertSector_Corr = new TH2D(hname, htitle, 300, -40, -10, 6, 0.5, 6.5);
    
    sprintf(hname,"elecZVert_VS_Phi");
    sprintf(htitle,"Z Vertex  vs. #phi, Electrons");
    elecZVert_VS_Phi = new TH2D(hname,htitle, 360, -180., 180., 300, -35., -20.);

    sprintf(hname,"elecZVert_VS_Phi_Corr");
    sprintf(htitle,"Z Vertex (corrected) vs. #phi, Electrons");
    elecZVert_VS_Phi_Corr = new TH2D(hname,htitle, 360, -180., 180., 300, -35., -20.);

    sprintf(hname,"Xvert");
    sprintf(htitle,"X vertex");
    Xvert = new TH2D(hname,htitle, 300, -10, 10,5,-0.5,4.5);

    sprintf(hname,"Yvert");
    sprintf(htitle,"Y vertex");
    Yvert = new TH2D(hname,htitle, 300, -10, 10,5,-0.5,4.5);
    
    sprintf(hname,"ZVertDiff");
    sprintf(htitle,"Difference Between Z Vertices of electron and other particle");
    ZVertDiff = new TH2D(hname,htitle, 300, -10, 10,4,0.5,4.5);
    
    sprintf(hname,"Beta_VS_Momentum");
    sprintf(htitle,"Beta vs Momentum");
	Beta_VS_Momentum = new TH2D(hname,htitle, 500, 0, 5, 115, 0, 1.15);

    sprintf(hname,"Beta_VS_Momentum_Recalc");
    sprintf(htitle,"Beta vs Momentum, Recalculated");
    Beta_VS_Momentum_Recalc = new TH2D(hname,htitle, 500, 0, 5, 115, 0, 1.15);
    
    sprintf(hname,"Beta_Recalc");
    sprintf(htitle,"Beta, Recalculated");
    Beta_Recalc = new TH2D(hname,htitle, 320, 0.4, 2.0, 5, -0.5, 4.5);

    sprintf(hname,"scMassSquared");
    sprintf(htitle,"TOF Mass^{2}");
    scMassSquared = new TH2D(hname,htitle, 200, -0.3, 0.3, 5, -0.5, 4.5);

    sprintf(hname,"scMassSquared_elecID");
    sprintf(htitle,"TOF Mass^{2}, e- ID cuts");
    scMassSquared_elecID = new TH2D(hname,htitle, 200, -0.3, 0.3, 5, -0.5, 4.5);
    
    sprintf(hname,"scMassSquared_photID");
    sprintf(htitle,"TOF Mass^{2}, photon ID cuts");
    scMassSquared_photID = new TH2D(hname,htitle, 200, -0.3, 0.3, 5, -0.5, 4.5);

    sprintf(hname,"scMassSquared_elecIDphotID");
    sprintf(htitle,"TOF Mass^{2}, e- and photon ID cuts");
    scMassSquared_elecIDphotID = new TH2D(hname,htitle, 200, -0.3, 0.3, 5, -0.5, 4.5);
    
    sprintf(hname,"TotalMomentum");
    sprintf(htitle,"Total Momentum");
	TotalMomentum = new TH2D(hname,htitle, 600, 0, 6, 7, -0.5, 6.5);
    
    sprintf(hname,"OpAng_2Photons");
    sprintf(htitle,"Opening Angle Between Photons");
	OpAng_2Photons = new TH1D(hname,htitle, 180, 0, 180);

    sprintf(hname,"OpAng_elecPhoton1");
    sprintf(htitle,"Opening Angle Between e^{-} and #gamma_{1}");
	OpAng_elecPhoton1 = new TH1D(hname,htitle, 180, 0, 180);

    sprintf(hname,"OpAng_elecPhoton2");
    sprintf(htitle,"Opening Angle Between e^{-} and #gamma_{2}");
	OpAng_elecPhoton2 = new TH1D(hname,htitle, 180, 0, 180);
    
    // particle ID histogram
    sprintf(hname,"CCnphe");
    sprintf(htitle,"CC Number of Photo-electrons");
    CCnphe = new TH2D(hname,htitle, 100, 0, 100, 5, -0.5, 4.5);

    sprintf(hname,"ECu");
    sprintf(htitle,"EC U-view");
    ECu = new TH2D(hname,htitle, 450, 0, 450, 5, -0.5, 4.5);
    
    sprintf(hname,"ECv");
    sprintf(htitle,"EC V-view");
    ECv = new TH2D(hname,htitle, 450, 0, 450, 5, -0.5, 4.5);
    
    sprintf(hname,"ECw");
    sprintf(htitle,"EC W-view");
    ECw = new TH2D(hname,htitle, 450, 0, 450, 5, -0.5, 4.5);
    
    sprintf(hname,"dtime_ECSC");
    sprintf(htitle,"#Delta t(EC-SC)");
    dtime_ECSC = new TH2D(hname,htitle, 100, -5.0, 5.0, 5, -0.5, 4.5);
    
    for(i=0; i<myPartList.Get_nPartLabel(); i++){
        sprintf(hname,"Theta_VS_Phi_%s",myPartList.Get_PartLabel(i).c_str());
        sprintf(htitle,"Theta vs Phi for %s",myPartList.Get_PartLabel(i).c_str());
        Theta_VS_Phi[i] = new TH2D(hname,htitle, 180, 0, 180, 360, -180, 180);
    }
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        sprintf(hname,"Xvert_VS_Yvert_%s",myDetPart.Get_DetPartLabel(i).c_str());
        sprintf(htitle,"X Vertex vs Y Vertex, %s",myDetPart.Get_DetPartLabel(i).c_str());
        Xvert_VS_Yvert[i] = new TH2D(hname,htitle, 100, -0.05, 0.05, 100, -0.05, 0.05);

    	sprintf(hname,"ECtot_VS_P_%s",myDetPart.Get_DetPartLabel(i).c_str());
    	sprintf(htitle,"ECtot vs P, %s",myDetPart.Get_DetPartLabel(i).c_str());
    	ECtot_VS_P[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 1.0);

    	sprintf(hname,"ECin_VS_ECout_%s",myDetPart.Get_DetPartLabel(i).c_str());
    	sprintf(htitle,"ECin vs ECout, %s",myDetPart.Get_DetPartLabel(i).c_str());
    	ECin_VS_ECout[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.25);

        sprintf(hname,"ECtotP_VS_P_%s",myDetPart.Get_DetPartLabel(i).c_str());
        sprintf(htitle,"ECtot/P vs P, %s",myDetPart.Get_DetPartLabel(i).c_str());
        ECtotP_VS_P[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);
        
        sprintf(hname,"dBeta_VS_Momentum_%s",myDetPart.Get_DetPartLabel(i).c_str());
        sprintf(htitle,"#Delta #beta vs Momentum, %s",myDetPart.Get_DetPartLabel(i).c_str());
        dBeta_VS_Momentum[i] = new TH2D(hname,htitle, 500, 0, 5, 200, -0.25, 0.25);

        sprintf(hname,"dBeta_VS_Momentum_EPC_%s",myDetPart.Get_DetPartLabel(i).c_str());
        sprintf(htitle,"#Delta #beta vs Momentum with e- and photon ID, %s",myDetPart.Get_DetPartLabel(i).c_str());
        dBeta_VS_Momentum_EPC[i] = new TH2D(hname,htitle, 500, 0, 5, 200, -0.25, 0.25);
    }

    // electron ID histogram
    sprintf(hname,"Mom_elecID");
    sprintf(htitle,"Momentum");
    Mom_elecID = new TH2D(hname,htitle, 500, 0, 5.0, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"CCnphe_elecID");
    sprintf(htitle,"CC Number of Photo-electrons");
    CCnphe_elecID = new TH2D(hname,htitle, 100, 0, 100, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"ECu_elecID");
    sprintf(htitle,"EC U-view");
    ECu_elecID = new TH2D(hname,htitle, 450, 0, 450, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"ECv_elecID");
    sprintf(htitle,"EC V-view");
    ECv_elecID = new TH2D(hname,htitle, 450, 0, 450, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"ECw_elecID");
    sprintf(htitle,"EC W-view");
    ECw_elecID = new TH2D(hname,htitle, 450, 0, 450, nElecID, -0.5, nElecID - 0.5);
    
    sprintf(hname,"dtime_ECSC_elecID");
    sprintf(htitle,"#Delta t(EC-SC)");
    dtime_ECSC_elecID = new TH2D(hname,htitle, 100, -5.0, 5.0, nElecID, -0.5, nElecID - 0.5);
    
    for(i=0; i<myElecID.Get_nElecID(); i++){
        sprintf(hname,"ECtot_VS_P_elecID_0%i",i);
        sprintf(htitle,"ECtot vs P, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECtot_VS_P_elecID[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 1.0);
        
        sprintf(hname,"ECin_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"ECin vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECin_VS_ECout_elecID[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.25);
        
        sprintf(hname,"ECtotP_VS_P_elecID_0%i",i);
        sprintf(htitle,"ECtot/P vs P, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECtotP_VS_P_elecID[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);

        sprintf(hname,"Mom_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"Mom. vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        Mom_VS_ECout_elecID[i] = new TH2D(hname,htitle, 500, 0, 5.0, 100, 0, 0.25);

        sprintf(hname,"ECu_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"EC U-view vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECu_VS_ECout_elecID[i] = new TH2D(hname,htitle, 450, 0, 450, 100, 0, 0.25);

        sprintf(hname,"ECv_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"EC V-view vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECv_VS_ECout_elecID[i] = new TH2D(hname,htitle, 450, 0, 450, 100, 0, 0.25);
        
        sprintf(hname,"ECw_VS_ECout_elecID_0%i",i);
        sprintf(htitle,"EC W-view vs ECout, %s",myElecID.Get_elecIDLabel(i).c_str());
        ECw_VS_ECout_elecID[i] = new TH2D(hname,htitle, 450, 0, 450, 100, 0, 0.25);
    }
    
    sprintf(hname,"ECin_VS_ECout_elecID_All");
    sprintf(htitle,"ECin vs ECout, all e- cuts");
    ECin_VS_ECout_elecID_All = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.25);
    
    sprintf(hname,"Theta_VS_Phi_ECoutCut");
    sprintf(htitle,"Theta vs Phi for e-, EC_{out} < 0.01 GeV");
    Theta_VS_Phi_ECoutCut = new TH2D(hname,htitle, 80, 0, 80, 360, -180, 180);
    
    sprintf(hname,"q2_ECoutCut");
    sprintf(htitle,"Q^{2}, EC_{out} < 0.01 GeV");
    q2_ECoutCut = new TH1D(hname,htitle, 100, 0., 4.);
    
    sprintf(hname,"elecZVert_ECoutCut");
    sprintf(htitle,"Z Vertex of Electron, EC_{out} < 0.01 GeV");
    elecZVert_ECoutCut = new TH1D(hname,htitle, 300, -35, -20);
    
    sprintf(hname,"Beta_VS_Momentum_ECoutCut");
    sprintf(htitle,"Beta vs Momentum, EC_{out} < 0.01 GeV");
    Beta_VS_Momentum_ECoutCut = new TH2D(hname,htitle, 500, 0, 5, 115, 0, 1.15);
    
    sprintf(hname,"ECtot_VS_P_ECoutCut");
    sprintf(htitle,"ECtot vs P, EC_{out} < 0.01 GeV");
    ECtot_VS_P_ECoutCut = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 1.0);
    
    sprintf(hname,"ECtotP_VS_P_ECoutCut");
    sprintf(htitle,"ECtot/P vs P, EC_{out} < 0.01 GeV");
    ECtotP_VS_P_ECoutCut = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);

    sprintf(hname,"ECtotMinusECin_ECoutCut");
    sprintf(htitle,"ECtot - ECin, EC_{out} < 0.01 GeV");
    ECtotMinusECin_ECoutCut = new TH1D(hname,htitle, 100,-1.0,1.0);
    
    sprintf(hname,"Theta_VS_Phi_AntiECoutCut");
    sprintf(htitle,"Theta vs Phi for e-, EC_{out} >= 0.01 GeV");
    Theta_VS_Phi_AntiECoutCut = new TH2D(hname,htitle, 80, 0, 80, 360, -180, 180);
    
    sprintf(hname,"q2_AntiECoutCut");
    sprintf(htitle,"Q^{2}, EC_{out} >= 0.01 GeV");
    q2_AntiECoutCut = new TH1D(hname,htitle, 100, 0., 4.);
    
    sprintf(hname,"elecZVert_AntiECoutCut");
    sprintf(htitle,"Z Vertex of Electron, EC_{out} >= 0.01 GeV");
    elecZVert_AntiECoutCut = new TH1D(hname,htitle, 300, -35, -20);
    
    sprintf(hname,"Beta_VS_Momentum_AntiECoutCut");
    sprintf(htitle,"Beta vs Momentum, EC_{out} >= 0.01 GeV");
    Beta_VS_Momentum_AntiECoutCut = new TH2D(hname,htitle, 500, 0, 5, 115, 0, 1.15);
    
    sprintf(hname,"ECtot_VS_P_AntiECoutCut");
    sprintf(htitle,"ECtot vs P, EC_{out} >= 0.01 GeV");
    ECtot_VS_P_AntiECoutCut = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 1.0);
    
    sprintf(hname,"ECtotP_VS_P_AntiECoutCut");
    sprintf(htitle,"ECtot/P vs P, EC_{out} >= 0.01 GeV");
    ECtotP_VS_P_AntiECoutCut = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);
    
    sprintf(hname,"ECtotMinusECin_AntiECoutCut");
    sprintf(htitle,"ECtot - ECin, EC_{out} >= 0.01 GeV");
    ECtotMinusECin_AntiECoutCut = new TH1D(hname,htitle, 100,-1.0,1.0);
    
    sprintf(hname,"ECin_VS_ECout_ECfid");
    sprintf(htitle,"ECin vs ECout, EC fid. cut");
    ECin_VS_ECout_ECfid = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.25);
    
    for(i=0; i<myTgt.Get_nIndex(); i++){
        sprintf(hname,"Xvert_VS_Yvert_AllCuts_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"X Vertex vs Y Vertex, All Cuts, %s",myTgt.Get_Label(i).c_str());
        Xvert_VS_Yvert_AllCuts[i] = new TH2D(hname,htitle, 100, -0.05, 0.05, 100, -0.05, 0.05);

        sprintf(hname,"Xvert_VS_Yvert_Omega_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"X Vertex vs Y Vertex, #omega, %s",myTgt.Get_Label(i).c_str());
        Xvert_VS_Yvert_Omega[i] = new TH2D(hname,htitle, 100, -0.05, 0.05, 100, -0.05, 0.05);

		sprintf(hname,"hW_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"W of Reaction, %s",myTgt.Get_Label(i).c_str());
		hW[i] = new TH1D(hname, htitle, 250, 0, 5);
        
        sprintf(hname,"hMx_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"M_{x} of Reaction, %s",myTgt.Get_Label(i).c_str());
        hMx[i] = new TH1D(hname, htitle, 250, 0, 5);
        
		sprintf(hname,"z_fracE_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Fractional Energy, %s",myTgt.Get_Label(i).c_str());
		z_fracE[i] = new TH1D(hname, htitle, 150, 0, 1.5);
        
		sprintf(hname,"LongMom_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Longitudinal Momentum of Reconstructed Particle, %s",myTgt.Get_Label(i).c_str());
		LongMom[i] = new TH1D(hname, htitle, 500, 0, 5);
        
		sprintf(hname,"TransMom_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Transverse Momentum of Reconstructed Particle, %s",myTgt.Get_Label(i).c_str());
		TransMom[i] = new TH1D(hname, htitle, 500, 0, 5);
        
		sprintf(hname,"IM2Photons_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of Pi0 - Single Cut, %s",myTgt.Get_Label(i).c_str());
		IM2Photons[i] = new TH2D(hname, htitle, 100, 0., 1., 9, -0.5, 8.5);
        
        sprintf(hname,"OpAng_VS_IM2Photons_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. Reconstructed Mass of Pi0, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_IM2Photons[i] = new TH2D(hname, htitle, 100, 0., 1., 100, 0, 100.);
        
		sprintf(hname,"OpAng_VS_E_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. #pi^{0} Total Energy, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_E[i] = new TH2D(hname, htitle, 350, 0., 3.5, 100, 0, 100.);

        sprintf(hname,"OpAng_VS_E_MassPi0Cut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs. #pi^{0} Total Energy with IM(#pi^{0}) Cut, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_E_MassPi0Cut[i] = new TH2D(hname, htitle, 350, 0., 3.5, 100, 0, 100.);

        sprintf(hname,"IM2Pions_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #pi^{+}#pi^{-} vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
        IM2Pions_VS_IMOmega[i] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

        sprintf(hname,"IM2Pions_VS_IMOmega_AllCuts_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #pi^{+}#pi^{-} vs Reconstructed Mass of #omega, all cuts, %s",myTgt.Get_Label(i).c_str());
        IM2Pions_VS_IMOmega_AllCuts[i] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"IM2Photons_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #pi^{0} vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		IM2Photons_VS_IMOmega[i] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);
        
        sprintf(hname,"W_VS_IMOmega_AllCuts_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"W vs Reconstructed Mass of #omega, All Cuts except W, %s",myTgt.Get_Label(i).c_str());
        W_VS_IMOmega_AllCuts[i] = new TH2D(hname, htitle, 250, 0., 5.0, nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"Q2_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Q2 vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		Q2_VS_IMOmega[i] = new TH2D(hname, htitle, 100, 0., 4.0, nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"Pt_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Omega Trans. Mom. vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		Pt_VS_IMOmega[i] = new TH2D(hname, htitle, 500, 0., 5., nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"Pl_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Omega Long. Mom. vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		Pl_VS_IMOmega[i] = new TH2D(hname, htitle, 500, 0., 5., nIMomega, IMomegaLo, IMomegaHi);
        
		sprintf(hname,"OpAng_VS_IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Opening Angle vs Reconstructed Mass of #omega, %s",myTgt.Get_Label(i).c_str());
		OpAng_VS_IMOmega[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, 100, 0., 100.);
        
		sprintf(hname,"MissMom_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Missing Momentum, %s",myTgt.Get_Label(i).c_str());
		MissMom[i] = new TH1D(hname, htitle, 600, 0, 6);
        
		sprintf(hname,"MMsq_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Missing Mass Squared, %s",myTgt.Get_Label(i).c_str());
		MMsq[i] = new TH1D(hname, htitle, 700, 0, 7);

        sprintf(hname,"IMOmega_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - Single Cut, %s",myTgt.Get_Label(i).c_str());
		IMOmega[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, 9, -0.5, 8.5);
        
        sprintf(hname,"IMOmega_woCut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"Reconstructed Mass of #omega - All Cuts, %s",myTgt.Get_Label(i).c_str());
		IMOmega_woCut[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, 9, -0.5, 8.5);

        sprintf(hname,"IMOmega_antiCut_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - Anti-Cuts, %s",myTgt.Get_Label(i).c_str());
        IMOmega_antiCut[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, 9, -0.5, 8.5);
        
        sprintf(hname,"PtSq_Omega_AllCuts_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"#omega Meson - All ID Cuts, %s",myTgt.Get_Label(i).c_str());
		PtSq_Omega_AllCuts[i] = new TH1D(hname, htitle, nPtSq_omega, PtSq_omegaLo, PtSq_omegaHi);

        sprintf(hname,"PtSq_Omega_AllCuts_IMOmegaCut_%s",myTgt.Get_Label(i).c_str());
		sprintf(htitle,"#omega Meson - All ID Cuts & IM(#omega) Cut, %s",myTgt.Get_Label(i).c_str());
		PtSq_Omega_AllCuts_IMOmegaCut[i] = new TH1D(hname, htitle, nPtSq_omega, PtSq_omegaLo, PtSq_omegaHi);

        sprintf(hname,"PtSq_Omega_AllCuts_IMOmegaSBCut_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"#omega Meson - All ID Cuts & IM(#omega) Cut, sideband, %s",myTgt.Get_Label(i).c_str());
        PtSq_Omega_AllCuts_IMOmegaSBCut[i] = new TH1D(hname, htitle, nPtSq_omega, PtSq_omegaLo, PtSq_omegaHi);
            
        sprintf(hname,"IM2Photons_ME%i_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of Pi0, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IM2Photons_ME[i] = new TH2D(hname, htitle, 100, 0., 1., nME_Methods, -0.5, nME_Methods-0.5);
            
        sprintf(hname,"IM2Photons_OpAng_ElecPhoton_Cut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of Pi0 with e^{-} #gamma Opening Angle Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IM2Photons_OpAng_ElecPhoton_Cut_ME[i] = new TH2D(hname, htitle, 100, 0., 1., nME_Methods, -0.5, nME_Methods-0.5);
            
        sprintf(hname,"IMOmega_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_OpAng_ElecPhoton_Cut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - e^{-} #gamma Opening Angle Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_OpAng_ElecPhoton_Cut_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_MassPi0Cut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - Pi0 Mass Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_MassPi0Cut_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_ZVertCut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - Z Vertex Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_ZVertCut_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_QsqCut_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - Q^{2} Cut, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_QsqCut_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
        
        sprintf(hname,"IMOmega_AllCuts_ME_%s",myTgt.Get_Label(i).c_str());
        sprintf(htitle,"Reconstructed Mass of #omega - All Cuts, Mixed Evt, %s",myTgt.Get_Label(i).c_str());
        IMOmega_AllCuts_ME[i] = new TH2D(hname, htitle, nIMomega, IMomegaLo, IMomegaHi, nME_Methods, -0.5, nME_Methods-0.5);
    }
    
    for(i=0; i < 5; i++) {
        sprintf(hname, "ECinP_VS_ECout/P_%.1f_P_%.1f",i*0.5+0.5, i*0.5+1.0);
        sprintf(htitle,"ECin/P vs ECout/P, %.1f < P <= %.1f",i*0.5+0.5, i*0.5+1.0);
        ECinP_VS_ECoutP_Range[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.5);
    }

    for(i=0; i<MAX_SECTORS; i++){
        sprintf(hname,"ECtotP_VS_P_Sector%i",i+1);
        sprintf(htitle,"ECtot/P vs P, Sector %i",i+1);
        ECtotP_VS_P_Sector[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);

        sprintf(hname,"ECinP_VS_ECoutP_Sector%i",i+1);
        sprintf(htitle,"ECin/P vs ECout/P, Sector %i",i+1);
        ECinP_VS_ECoutP[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.5);

        sprintf(hname,"ECinP_VS_ECoutP_cut_Sector%i",i+1);
        sprintf(htitle,"ECin/P vs ECout/P cut, Sector %i",i+1);
        ECinP_VS_ECoutP_cut[i] = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.5);

        sprintf(hname,"ECtotP_VS_P_ECPCut%i",i+1);
        sprintf(htitle,"ECtot/P vs P, Sector %i",i+1);
        ECtotP_VS_P_ECPCut[i] = new TH2D(hname,htitle, 500, 0, 5, 100, 0, 0.5);
        
        sprintf(hname,"EC_XvsY_local_Sector%i",i+1);
        sprintf(htitle,"EC local X vs local Y, Sector %i",i+1);
        EC_XvsY_local_Sector[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_ECoutCut%i",i+1);
        sprintf(htitle,"EC local X vs local Y, EC_{out} cut, Sector %i",i+1);
        EC_XvsY_local_ECoutCut[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);

        sprintf(hname,"EC_XvsY_local_FidCut%i",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut, Sector %i",i+1);
        EC_XvsY_local_FidCut[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);

        sprintf(hname,"EC_XvsY_local_AntiFidCut%i",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut (anti), Sector %i",i+1);
        EC_XvsY_local_AntiFidCut[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
    }
    
    sprintf(hname,"RelativityOpAngPhotonsA");
    sprintf(htitle,"0/180 Degree Decay Photons");
	RelativityOpAngPhotonsA = new TH2D(hname,htitle, 500, 0, 5, 200, 0, 200);

    sprintf(hname,"RelativityOpAngPhotonsB");
    sprintf(htitle,"90/-90 Degree Decay Photons");
	RelativityOpAngPhotonsB = new TH2D(hname,htitle, 500, 0, 5, 180, 0, 180);
    
    sprintf(hname,"GammaPi0");
    sprintf(htitle,"Gamma of Pi0");
	GammaPi0 = new TH1D(hname,htitle, 240, 1, 25);

    sprintf(hname,"BetaPi0");
    sprintf(htitle,"Beta of Pi0");
	BetaPi0 = new TH1D(hname,htitle, 100, 0, 1);

    sprintf(hname,"MomentumPhoton1");
    sprintf(htitle,"Total Momentum of Photon 1");
	MomentumPhoton1 = new TH1D(hname,htitle, 500, 0, 5);

    sprintf(hname,"MomentumPhoton2");
    sprintf(htitle,"Total Momentum of Photon 2");
	MomentumPhoton2 = new TH1D(hname,htitle, 500, 0, 5);

    sprintf(hname,"MomentumPhoton1_cut");
    sprintf(htitle,"Total Momentum of Photon 1 Cut");
	MomentumPhoton1_cut = new TH1D(hname,htitle, 500, 0, 5);

    sprintf(hname,"MomentumPhoton2_cut");
    sprintf(htitle,"Total Momentum of Photon 2 Cut");
	MomentumPhoton2_cut = new TH1D(hname,htitle, 500, 0, 5);

    sprintf(hname,"BetaPhoton1");
    sprintf(htitle,"Beta of Photon 1");
	BetaPhoton1 = new TH1D(hname,htitle, 100, 0.8, 2.1);

    sprintf(hname,"BetaPhoton2");
    sprintf(htitle,"Beta of Photon 2");
	BetaPhoton2 = new TH1D(hname,htitle, 100, 0.8, 2.1);

    sprintf(hname,"BetaPhoton1_cut");
    sprintf(htitle,"Beta of Photon 1 Cut");
	BetaPhoton1_cut = new TH1D(hname,htitle, 100, 0.8, 2.1);

    sprintf(hname,"BetaPhoton2_cut");
    sprintf(htitle,"Beta of Photon 2 Cut");
	BetaPhoton2_cut = new TH1D(hname,htitle, 100, 0.8, 2.1);

    sprintf(hname,"ECuPhoton1");
    sprintf(htitle,"ECu of Photon 1");
	ECuPhoton1 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECuPhoton2");
    sprintf(htitle,"ECu of Photon 2");
	ECuPhoton2 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECvPhoton1");
    sprintf(htitle,"ECv of Photon 1");
	ECvPhoton1 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECvPhoton2");
    sprintf(htitle,"ECv of Photon 2");
	ECvPhoton2 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECwPhoton1");
    sprintf(htitle,"ECw of Photon 1");
	ECwPhoton1 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECwPhoton2");
    sprintf(htitle,"ECw of Photon 2");
	ECwPhoton2 = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECuPhoton1_cut");
    sprintf(htitle,"ECu of Photon 1 Cut");
	ECuPhoton1_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECuPhoton2_cut");
    sprintf(htitle,"ECu of Photon 2 Cut");
	ECuPhoton2_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECvPhoton1_cut");
    sprintf(htitle,"ECv of Photon 1 Cut");
	ECvPhoton1_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECvPhoton2_cut");
    sprintf(htitle,"ECv of Photon 2 Cut");
	ECvPhoton2_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECwPhoton1_cut");
    sprintf(htitle,"ECw of Photon 1 Cut");
	ECwPhoton1_cut = new TH1D(hname,htitle, 450, 0, 450);

    sprintf(hname,"ECwPhoton2_cut");
    sprintf(htitle,"ECw of Photon 2 Cut");
	ECwPhoton2_cut = new TH1D(hname,htitle, 450, 0, 450);

    for(i=0; i<MAX_SECTORS; i++){
        sprintf(hname,"EC_XvsY_local_Sector_%i_Photon1",i+1);
        sprintf(htitle,"EC local X vs local Y, Sector %i, Photon 1",i+1);
        EC_XvsY_local_Sector_Photon1[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_FidCut%i_Photon1",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut, Sector %i, Photon 1",i+1);
        EC_XvsY_local_FidCut_Photon1[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_AntiFidCut%i_Photon1",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut (anti), Sector %i, Photon1",i+1);
        EC_XvsY_local_AntiFidCut_Photon1[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);

        sprintf(hname,"EC_XvsY_local_Sector_%i_Photon2",i+1);
        sprintf(htitle,"EC local X vs local Y, Sector %i, Photon 2",i+1);
        EC_XvsY_local_Sector_Photon2[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_FidCut%i_Photon2",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut, Sector %i, Photon 2",i+1);
        EC_XvsY_local_FidCut_Photon2[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
        
        sprintf(hname,"EC_XvsY_local_AntiFidCut%i_Photon2",i+1);
        sprintf(htitle,"EC local X vs local Y, EC fid. cut (anti), Sector %i, Photon2",i+1);
        EC_XvsY_local_AntiFidCut_Photon2[i] = new TH2D(hname,htitle, 200, -200, 200, 200, -200,200);
    }
    
    sprintf(hname,"ECtime_ECl_Start_Photon1");
    sprintf(htitle,"ECtime - StartTime - EClength/c of Photon 1");
    ECtime_ECl_Start_Photon1 = new TH1D(hname,htitle, 400, -80., 20.);
    
    sprintf(hname,"ECtime_ECl_Start_Photon2");
    sprintf(htitle,"ECtime - StartTime - EClength/c of Photon 2");
    ECtime_ECl_Start_Photon2 = new TH1D(hname,htitle, 400, -80., 20.);
    
    sprintf(hname,"ECtime_ECl_Photon1");
    sprintf(htitle,"ECtime - EClength/c of Photon 1");
    ECtime_ECl_Photon1 = new TH1D(hname,htitle, 200, -3, 2);

    sprintf(hname,"ECtime_ECl_Photon2");
    sprintf(htitle,"ECtime - EClength/c of Photon 2");
    ECtime_ECl_Photon2 = new TH1D(hname,htitle, 200, -3, 2);

    sprintf(hname,"ECtime_ECl_Photon1_cut");
    sprintf(htitle,"ECtime - EClength/c of Photon 1 Cut");
    ECtime_ECl_Photon1_cut = new TH1D(hname,htitle, 200, -3, 2);

    sprintf(hname,"ECtime_ECl_Photon2_cut");
    sprintf(htitle,"ECtime - EClength/c of Photon 2 Cut");
    ECtime_ECl_Photon2_cut = new TH1D(hname,htitle, 200, -3, 2);

    sprintf(hname,"ECtimePhoton1");
    sprintf(htitle,"ECtime of Photon 1");
    ECtimePhoton1 = new TH1D(hname,htitle, 300, 10, 25);

    sprintf(hname,"ECtimePhoton2");
    sprintf(htitle,"ECtime of Photon 2");
    ECtimePhoton2 = new TH1D(hname,htitle, 300, 10, 25);

    sprintf(hname,"ECpathPhoton1");
    sprintf(htitle,"ECpath of Photon 1");
    ECpathPhoton1 = new TH1D(hname,htitle, 200, 450, 650);

    sprintf(hname,"ECpathPhoton2");
    sprintf(htitle,"ECpath of Photon 2");
    ECpathPhoton2 = new TH1D(hname,htitle, 200, 450, 650);

    sprintf(hname,"ECpathtimePhoton1");
    sprintf(htitle,"ECpath/c of Photon 1");
    ECpathtimePhoton1 = new TH1D(hname,htitle, 300, 17, 20);
    
    sprintf(hname,"ECpathtimePhoton2");
    sprintf(htitle,"ECpath/c of Photon 2");
    ECpathtimePhoton2 = new TH1D(hname,htitle, 300, 17, 20);
    
    sprintf(hname,"ECtotP_vs_P_Photon1");
    sprintf(htitle,"ECtot / P vs P of Photon 1");
    ECtotP_vs_P_Photon1 = new TH2D(hname,htitle, 500, 0, 5, 500, 0, 0.5);

    sprintf(hname,"ECtotP_vs_P_Photon2");
    sprintf(htitle,"ECtot / P vs P of Photon 2");
    ECtotP_vs_P_Photon2 = new TH2D(hname,htitle, 500, 0, 5, 500, 0, 0.5);

    sprintf(hname,"ECin_vs_ECout_Photon1");
    sprintf(htitle,"EC_{in} vs EC_{out} of Photon 1");
    ECin_vs_ECout_Photon1 = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.35);

    sprintf(hname,"ECin_vs_ECout_Photon2");
    sprintf(htitle,"EC_{in} vs EC_{out} of Photon 2");
    ECin_vs_ECout_Photon2 = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.35);

    sprintf(hname,"ECtotP_vs_P_InOutZeroCut_Photon1");
    sprintf(htitle,"ECtot / P vs P of Photon 1");
    ECtotP_vs_P_InOutZeroCut_Photon1 = new TH2D(hname,htitle, 500, 0, 5, 500, 0, 0.5);
    
    sprintf(hname,"ECtotP_vs_P_InOutZeroCut_Photon2");
    sprintf(htitle,"ECtot / P vs P of Photon 2");
    ECtotP_vs_P_InOutZeroCut_Photon2 = new TH2D(hname,htitle, 500, 0, 5, 500, 0, 0.5);
    
    sprintf(hname,"ECin_vs_ECout_InOutZeroCut_Photon1");
    sprintf(htitle,"EC_{in} vs EC_{out} of Photon 1");
    ECin_vs_ECout_InOutZeroCut_Photon1 = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.35);
    
    sprintf(hname,"ECin_vs_ECout_InOutZeroCut_Photon2");
    sprintf(htitle,"EC_{in} vs EC_{out} of Photon 2");
    ECin_vs_ECout_InOutZeroCut_Photon2 = new TH2D(hname,htitle, 100, 0, 0.5, 100, 0, 0.35);
    
    sprintf(hname,"Recon_Pi0_Mass_PhotID_Cuts");
    sprintf(htitle,"Reconstructed #pi^{0} Mass - Photon ID Cuts");
    Pi0Mass_PhotIDcuts = new TH1D(hname,htitle, 100, 0, 1);

    sprintf(hname,"Recon_Omega_Mass_All_Cuts");
    sprintf(htitle,"Reconstructed #omega Mass - All Cuts");
    OmegaMass_AllCuts = new TH1D(hname,htitle,100,0,3);

    sprintf(hname,"mass2Pions_VS_massOmega_NC1");
    sprintf(htitle,"IM #pi^{+} #pi^{-} vs IM #pi^{+} #pi^{-} #gamma #gamma - No Cuts");
    mass2Pions_VS_massOmega_NC[0] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_NC2");
    sprintf(htitle,"IM #pi^{+} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - No Cuts");
    mass2Pions_VS_massOmega_NC[1] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_NC3");
    sprintf(htitle,"IM #pi^{-} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - No Cuts");
    mass2Pions_VS_massOmega_NC[2] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);
//TODO
    sprintf(hname,"mass2Pions_VS_massOmega_NCX1");
    sprintf(htitle,"IM #pi^{+} #pi^{-} vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection No Cuts");
    mass2Pions_VS_massOmega_NCX[0] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_NCX2");
    sprintf(htitle,"IM #pi^{+} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection No Cuts");
    mass2Pions_VS_massOmega_NCX[1] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_NCX3");
    sprintf(htitle,"IM #pi^{-} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection No Cuts");
    mass2Pions_VS_massOmega_NCX[2] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);
//endtodo

    sprintf(hname,"mass2Pions_VS_massOmega_EPC1");
    sprintf(htitle,"IM #pi^{+} #pi^{-} vs IM #pi^{+} #pi^{-} #gamma #gamma - Electron and Photon Cuts");
    mass2Pions_VS_massOmega_EPC[0] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_EPC2");
    sprintf(htitle,"IM #pi^{+} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - Electron and Photon Cuts");
    mass2Pions_VS_massOmega_EPC[1] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_EPC3");
    sprintf(htitle,"IM #pi^{-} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - Electron and Photon Cuts");
    mass2Pions_VS_massOmega_EPC[2] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

//TODO
    sprintf(hname,"mass2Pions_VS_massOmega_EPCX1");
    sprintf(htitle,"IM #pi^{+} #pi^{-} vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection Electron and Photon Cuts");
    mass2Pions_VS_massOmega_EPCX[0] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_EPCX2");
    sprintf(htitle,"IM #pi^{+} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection Electron and Photon Cuts");
    mass2Pions_VS_massOmega_EPCX[1] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_EPCX3");
    sprintf(htitle,"IM #pi^{-} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection Electron and Photon Cuts");
    mass2Pions_VS_massOmega_EPCX[2] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);
//endtodo

    sprintf(hname,"mass2Pions_VS_massOmega_EPOC1");
    sprintf(htitle,"IM #pi^{+} #pi^{-} vs IM #pi^{+} #pi^{-} #gamma #gamma - Electron, Photon, and Omega Cuts");
    mass2Pions_VS_massOmega_EPOC[0] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_EPOC2");
    sprintf(htitle,"IM #pi^{+} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - Electron, Photon, and Omega Cuts");
    mass2Pions_VS_massOmega_EPOC[1] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_EPOC3");
    sprintf(htitle,"IM #pi^{-} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - Electron, Photon, and Omega Cuts");
    mass2Pions_VS_massOmega_EPOC[2] = new TH2D(hname, htitle, 100, 0, 1., nIMomega, IMomegaLo, IMomegaHi);

//TODO
    sprintf(hname,"mass2Pions_VS_massOmega_EPOCX1");
    sprintf(htitle,"IM #pi^{+} #pi^{-} vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection Electron, Photon, and Omega Cuts");
    mass2Pions_VS_massOmega_EPOCX[0] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_EPOCX2");
    sprintf(htitle,"IM #pi^{+} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection Electron, Photon, and Omega Cuts");
    mass2Pions_VS_massOmega_EPOCX[1] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);

    sprintf(hname,"mass2Pions_VS_massOmega_EPOCX3");
    sprintf(htitle,"IM #pi^{-} #gamma #gamma vs IM #pi^{+} #pi^{-} #gamma #gamma - X Projection Electron, Photon, and Omega Cuts");
    mass2Pions_VS_massOmega_EPOCX[2] = new TH1D(hname, htitle, 100, 0, 1.);//, nIMomega, IMomegaLo, IMomegaHi);
//endtodo
}

//
// WriteHist - routine to write histograms to the output file
//
void HistManager::WriteHist(string RootFile){
    
    int i, j;
    
    DetectedParticles myDetPart;
    ParticleList myPartList;
    EG2Target myTgt;
    ElectronID myElecID;
    
	TFile *out = new TFile(RootFile.c_str(), "recreate");
	out->cd();
  
    // create a directory for check on kinematics
    TDirectory *cdKine = out->mkdir("Kinematics");
    cdKine->cd();

    StartTime->GetXaxis()->SetTitle("Event Start Time (ns)");
    StartTime->GetYaxis()->SetTitle("Counts");
    StartTime->Write();
    
    q2->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    q2->GetYaxis()->SetTitle("Counts");
	q2->Write();

    q2_VS_theta->GetYaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    q2_VS_theta->GetXaxis()->SetTitle("4E_{e'}sin^{2}(0.5*#theta_{e'})");
	q2_VS_theta->Write();
    
    nu_EnergyTransfer->GetXaxis()->SetTitle("\nu (GeV)");
    nu_EnergyTransfer->GetYaxis()->SetTitle("Counts");
	nu_EnergyTransfer->Write();
    
    elecZVert->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert->GetYaxis()->SetTitle("Counts");
	elecZVert->Write();

    elecZVertSector->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVertSector->GetYaxis()->SetTitle("Sector");
    elecZVertSector->Write();

    elecZVertSector_Corr->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVertSector_Corr->GetYaxis()->SetTitle("Sector");
    elecZVertSector_Corr->Write();
    
    elecZVert_VS_Phi->GetXaxis()->SetTitle("#phi (deg.)");
    elecZVert_VS_Phi->GetYaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_VS_Phi->Write();

    elecZVert_VS_Phi_Corr->GetXaxis()->SetTitle("#phi (deg.)");
    elecZVert_VS_Phi_Corr->GetYaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_VS_Phi_Corr->Write();
    
    Xvert->GetXaxis()->SetTitle("X vertex (cm)");
    Xvert->GetYaxis()->SetTitle("Particle");
    Xvert->Write();

    Yvert->GetXaxis()->SetTitle("Y vertex (cm)");
    Yvert->GetYaxis()->SetTitle("Particle");
    Yvert->Write();
    
    ZVertDiff->GetXaxis()->SetTitle(" #Delta Z (cm)");
    ZVertDiff->GetYaxis()->SetTitle("Particle");
	ZVertDiff->Write();
    
    Beta_VS_Momentum->GetXaxis()->SetTitle("Momentum (GeV/c)");
    Beta_VS_Momentum->GetYaxis()->SetTitle("#beta");
	Beta_VS_Momentum->Write();

    Beta_VS_Momentum_Recalc->GetXaxis()->SetTitle("Momentum (GeV/c)");
    Beta_VS_Momentum_Recalc->GetYaxis()->SetTitle("#beta");
    Beta_VS_Momentum_Recalc->Write();

    Beta_Recalc->GetXaxis()->SetTitle("#beta");
    Beta_Recalc->GetYaxis()->SetTitle("Particle");
    Beta_Recalc->Write();

    scMassSquared->GetXaxis()->SetTitle("TOF M^{2} (GeV^{2}/c^{2})");
    scMassSquared->GetYaxis()->SetTitle("Particle");
    scMassSquared->Write();

    scMassSquared_elecID->GetXaxis()->SetTitle("TOF M^{2} (GeV^{2}/c^{2})");
    scMassSquared_elecID->GetYaxis()->SetTitle("Particle");
    scMassSquared_elecID->Write();
    
    scMassSquared_photID->GetXaxis()->SetTitle("TOF M^{2} (GeV^{2}/c^{2})");
    scMassSquared_photID->GetYaxis()->SetTitle("Particle");
    scMassSquared_photID->Write();
    
    scMassSquared_elecIDphotID->GetXaxis()->SetTitle("TOF M^{2} (GeV^{2}/c^{2})");
    scMassSquared_elecIDphotID->GetYaxis()->SetTitle("Particle");
    scMassSquared_elecIDphotID->Write();
    
    TotalMomentum->GetXaxis()->SetTitle("Momentum (GeV/c)");
    TotalMomentum->GetYaxis()->SetTitle("Particle");
	TotalMomentum->Write();
    
    OpAng_2Photons->GetXaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
    OpAng_2Photons->GetYaxis()->SetTitle("Counts");
	OpAng_2Photons->Write();

    OpAng_elecPhoton1->GetXaxis()->SetTitle("Opening Angle between e^{-} and #gamma_{1} (deg.)");
    OpAng_elecPhoton1->GetYaxis()->SetTitle("Counts");
	OpAng_elecPhoton1->Write();

    OpAng_elecPhoton2->GetXaxis()->SetTitle("Opening Angle between e^{-} and #gamma_{2} (deg.)");
    OpAng_elecPhoton2->GetYaxis()->SetTitle("Counts");
    OpAng_elecPhoton2->Write();
 
    for(i=0; i<myPartList.Get_nPartLabel(); i++){
        Theta_VS_Phi[i]->GetXaxis()->SetTitle("#theta (deg.)");
        Theta_VS_Phi[i]->GetYaxis()->SetTitle("#phi (deg.)");
        Theta_VS_Phi[i]->Write();
    }
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        Xvert_VS_Yvert[i]->GetXaxis()->SetTitle("X vertex (cm)");
        Xvert_VS_Yvert[i]->GetYaxis()->SetTitle("Y vertex (cm)");
        Xvert_VS_Yvert[i]->Write();
    }

    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        dBeta_VS_Momentum[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        dBeta_VS_Momentum[i]->GetYaxis()->SetTitle("#Delta #beta");
        dBeta_VS_Momentum[i]->Write();
    }

    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        dBeta_VS_Momentum_EPC[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        dBeta_VS_Momentum_EPC[i]->GetYaxis()->SetTitle("#Delta #beta");
        dBeta_VS_Momentum_EPC[i]->Write();
    }
    
    // create a directory for check on detector info
    TDirectory *cdDetectors = out->mkdir("Detectors");
    cdDetectors->cd();
    
    CCnphe->GetXaxis()->SetTitle("Number of Photo-electrons");
    CCnphe->GetYaxis()->SetTitle("Particle");
    CCnphe->Write();

    ECu->GetXaxis()->SetTitle("EC U (cm)");
    ECu->GetYaxis()->SetTitle("Particle");
    ECu->Write();
    
    ECv->GetXaxis()->SetTitle("EC V (cm)");
    ECv->GetYaxis()->SetTitle("Particle");
    ECv->Write();
    
    ECw->GetXaxis()->SetTitle("EC W (cm)");
    ECw->GetYaxis()->SetTitle("Particle");
    ECw->Write();

    dtime_ECSC->GetXaxis()->SetTitle("#Delta t(EC-SC) (ns)");
    dtime_ECSC->GetYaxis()->SetTitle("Particle");
    dtime_ECSC->Write();
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        ECtot_VS_P[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
    	ECtot_VS_P[i]->GetYaxis()->SetTitle("EC total energy");
        ECtot_VS_P[i]->Write();
    }
    
    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        ECtotP_VS_P[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P[i]->Write();
    }

    for(i=0; i<myDetPart.Get_nDetPartLabel(); i++){
        ECin_VS_ECout[i]->GetXaxis()->SetTitle("EC inner energy");
    	ECin_VS_ECout[i]->GetYaxis()->SetTitle("EC outer energy");
        ECin_VS_ECout[i]->Write();
    }

    // create a directory for electron ID
    TDirectory *cdTgt[myTgt.Get_nIndex()];
    
    for(i=0; i<myTgt.Get_nIndex(); i++){
        cdTgt[i] = out->mkdir(myTgt.Get_Label(i).c_str());
        cdTgt[i]->cd();

        Xvert_VS_Yvert_AllCuts[i]->GetXaxis()->SetTitle("X vertex (cm)");
        Xvert_VS_Yvert_AllCuts[i]->GetYaxis()->SetTitle("Y vertex (cm)");
        Xvert_VS_Yvert_AllCuts[i]->Write();

        Xvert_VS_Yvert_Omega[i]->GetXaxis()->SetTitle("X vertex (cm)");
        Xvert_VS_Yvert_Omega[i]->GetYaxis()->SetTitle("Y vertex (cm)");
        Xvert_VS_Yvert_Omega[i]->Write();
        
        hMx[i]->GetXaxis()->SetTitle("M_{x} (GeV)");
        hMx[i]->GetYaxis()->SetTitle("Counts");
        hMx[i]->Write();
        
        hW[i]->GetXaxis()->SetTitle("W (GeV)");
        hW[i]->GetYaxis()->SetTitle("Counts");
        hW[i]->Write();

        z_fracE[i]->GetXaxis()->SetTitle("z");
        z_fracE[i]->GetYaxis()->SetTitle("Counts");
        z_fracE[i]->Write();
        
        LongMom[i]->GetXaxis()->SetTitle("#omega Longitudinal Momentum (GeV/c)");
        LongMom[i]->GetYaxis()->SetTitle("Counts");
		LongMom[i]->Write();

        TransMom[i]->GetXaxis()->SetTitle("#omega Transverse Momentum (GeV/c)");
        TransMom[i]->GetYaxis()->SetTitle("Counts");
		TransMom[i]->Write();

        IM2Photons[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons[i]->GetYaxis()->SetTitle("Cut Index");
        IM2Photons[i]->Write();

        OpAng_VS_IM2Photons[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        OpAng_VS_IM2Photons[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
		OpAng_VS_IM2Photons[i]->Write();
        
        OpAng_VS_E[i]->GetXaxis()->SetTitle("#pi^{0} Total Energy (GeV)");
        OpAng_VS_E[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
        OpAng_VS_E[i]->Write();
        
        OpAng_VS_E_MassPi0Cut[i]->GetXaxis()->SetTitle("#pi^{0} Total Energy (GeV)");
        OpAng_VS_E_MassPi0Cut[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
        OpAng_VS_E_MassPi0Cut[i]->Write();

        IM2Pions_VS_IMOmega[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} Inv. Mass (GeV/c^{2})");
        IM2Pions_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        IM2Pions_VS_IMOmega[i]->Write();

        IM2Pions_VS_IMOmega_AllCuts[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} Inv. Mass (GeV/c^{2})");
        IM2Pions_VS_IMOmega_AllCuts[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        IM2Pions_VS_IMOmega_AllCuts[i]->Write();
        
        IM2Photons_VS_IMOmega[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        IM2Photons_VS_IMOmega[i]->Write();

        W_VS_IMOmega_AllCuts[i]->GetXaxis()->SetTitle("W (GeV)");
        W_VS_IMOmega_AllCuts[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        W_VS_IMOmega_AllCuts[i]->Write();
        
        Q2_VS_IMOmega[i]->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
        Q2_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        Q2_VS_IMOmega[i]->Write();
        
        Pt_VS_IMOmega[i]->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");
        Pt_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
		Pt_VS_IMOmega[i]->Write();

        Pl_VS_IMOmega[i]->GetXaxis()->SetTitle("Longitudinal Momentum (GeV/c)");
        Pl_VS_IMOmega[i]->GetYaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
		Pl_VS_IMOmega[i]->Write();
		
        OpAng_VS_IMOmega[i]->GetXaxis()->SetTitle("#omega Inv. Mass (GeV/c^{2})");
        OpAng_VS_IMOmega[i]->GetYaxis()->SetTitle("Opening Angle between #gamma_{1} and #gamma_{2} (deg.)");
        OpAng_VS_IMOmega[i]->Write();
        
        MissMom[i]->GetXaxis()->SetTitle("Missing Momentum (GeV/c)");
        MissMom[i]->GetYaxis()->SetTitle("Counts");
		MissMom[i]->Write();
        
        MMsq[i]->GetXaxis()->SetTitle("Missing Mass Squared (GeV/c)^{2}");
        MMsq[i]->GetYaxis()->SetTitle("Counts");
		MMsq[i]->Write();

        IMOmega[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega[i]->GetYaxis()->SetTitle("Cut index");
		IMOmega[i]->Write();
        
        IMOmega_woCut[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_woCut[i]->GetYaxis()->SetTitle("Cut Index");
        IMOmega_woCut[i]->Write();

        IMOmega_antiCut[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_antiCut[i]->GetYaxis()->SetTitle("Cut Index");
        IMOmega_antiCut[i]->Write();
        
        PtSq_Omega_AllCuts[i]->GetXaxis()->SetTitle("P^{2}_{T}  (GeV/c)^{2}");
        PtSq_Omega_AllCuts[i]->GetYaxis()->SetTitle("Counts");
        PtSq_Omega_AllCuts[i]->Write();

        PtSq_Omega_AllCuts_IMOmegaCut[i]->GetXaxis()->SetTitle("P^{2}_{T}  (GeV/c)^{2}");
        PtSq_Omega_AllCuts_IMOmegaCut[i]->GetYaxis()->SetTitle("Counts");
        PtSq_Omega_AllCuts_IMOmegaCut[i]->Write();

        PtSq_Omega_AllCuts_IMOmegaSBCut[i]->GetXaxis()->SetTitle("P^{2}_{T}  (GeV/c)^{2}");
        PtSq_Omega_AllCuts_IMOmegaSBCut[i]->GetYaxis()->SetTitle("Counts");
        PtSq_Omega_AllCuts_IMOmegaSBCut[i]->Write();
        
        IM2Photons_ME[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IM2Photons_ME[i]->Write();
            
        IM2Photons_OpAng_ElecPhoton_Cut_ME[i]->GetXaxis()->SetTitle("#gamma #gamma Inv. Mass (GeV/c^{2})");
        IM2Photons_OpAng_ElecPhoton_Cut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IM2Photons_OpAng_ElecPhoton_Cut_ME[i]->Write();
            
        IMOmega_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_ME[i]->Write();
            
        IMOmega_OpAng_ElecPhoton_Cut_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_OpAng_ElecPhoton_Cut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_OpAng_ElecPhoton_Cut_ME[i]->Write();
            
        IMOmega_MassPi0Cut_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_MassPi0Cut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_MassPi0Cut_ME[i]->Write();
            
        IMOmega_ZVertCut_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_ZVertCut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_ZVertCut_ME[i]->Write();
            
        IMOmega_QsqCut_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_QsqCut_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_QsqCut_ME[i]->Write();
            
        IMOmega_AllCuts_ME[i]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} #gamma #gamma Inv. Mass (GeV/c^{2})");
        IMOmega_AllCuts_ME[i]->GetYaxis()->SetTitle("Mixing Method");
        IMOmega_AllCuts_ME[i]->Write();
    }
    
    out->cd();

    // create a directory for check on relativisitic kinematics
    TDirectory *cdRel = out->mkdir("Relativity");
    cdRel->cd();
    
    RelativityOpAngPhotonsA->Write();
	RelativityOpAngPhotonsB->Write();
	BetaPi0->Write();
	GammaPi0->Write();

    // create a directory for electron ID
    TDirectory *cdElecID = out->mkdir("ElectronID");
    cdElecID->cd();
    
    Mom_elecID->GetXaxis()->SetTitle("Momentum (GeV)");
    Mom_elecID->GetYaxis()->SetTitle("e- ID Cut");
    Mom_elecID->Write();
    
    CCnphe_elecID->GetXaxis()->SetTitle("Number of Photo-electrons");
    CCnphe_elecID->GetYaxis()->SetTitle("e- ID Cut");
    CCnphe_elecID->Write();
    
    ECu_elecID->GetXaxis()->SetTitle("EC U (cm)");
    ECu_elecID->GetYaxis()->SetTitle("e- ID Cut");
    ECu_elecID->Write();
    
    ECv_elecID->GetXaxis()->SetTitle("EC V (cm)");
    ECv_elecID->GetYaxis()->SetTitle("e- ID Cut");
    ECv_elecID->Write();
    
    ECw_elecID->GetXaxis()->SetTitle("EC W (cm)");
    ECw_elecID->GetYaxis()->SetTitle("e- ID Cut");
    ECw_elecID->Write();
    
    dtime_ECSC_elecID->GetXaxis()->SetTitle("#Delta t(EC-SC) (ns)");
    dtime_ECSC_elecID->GetYaxis()->SetTitle("e- ID Cut");
    dtime_ECSC_elecID->Write();
    
    for(i=0; i<myElecID.Get_nElecID(); i++){
        ECtot_VS_P_elecID[i]->GetXaxis()->SetTitle("Momentum (GeV)");
        ECtot_VS_P_elecID[i]->GetYaxis()->SetTitle("EC total energy");
        ECtot_VS_P_elecID[i]->Write();
        
        ECtotP_VS_P_elecID[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P_elecID[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P_elecID[i]->Write();
        
        ECin_VS_ECout_elecID[i]->GetXaxis()->SetTitle("EC inner energy");
        ECin_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        ECin_VS_ECout_elecID[i]->Write();
        
        Mom_VS_ECout_elecID[i]->GetXaxis()->SetTitle("Momentum (GeV)");
        Mom_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        Mom_VS_ECout_elecID[i]->Write();

        ECu_VS_ECout_elecID[i]->GetXaxis()->SetTitle("EC U (cm)");
        ECu_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        ECu_VS_ECout_elecID[i]->Write();

        ECv_VS_ECout_elecID[i]->GetXaxis()->SetTitle("EC V (cm)");
        ECv_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        ECv_VS_ECout_elecID[i]->Write();
        
        ECw_VS_ECout_elecID[i]->GetXaxis()->SetTitle("EC W (cm)");
        ECw_VS_ECout_elecID[i]->GetYaxis()->SetTitle("EC outer energy");
        ECw_VS_ECout_elecID[i]->Write();
    }
    
    q2_ECoutCut->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    q2_ECoutCut->GetYaxis()->SetTitle("Counts");
    q2_ECoutCut->Write();
    
    elecZVert_ECoutCut->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_ECoutCut->GetYaxis()->SetTitle("Counts");
    elecZVert_ECoutCut->Write();
    
    Beta_VS_Momentum_ECoutCut->GetXaxis()->SetTitle("Momentum (GeV/c)");
    Beta_VS_Momentum_ECoutCut->GetYaxis()->SetTitle("#beta");
    Beta_VS_Momentum_ECoutCut->Write();
    
    Theta_VS_Phi_ECoutCut->GetXaxis()->SetTitle("#theta (deg.)");
    Theta_VS_Phi_ECoutCut->GetYaxis()->SetTitle("#phi (deg.)");
    Theta_VS_Phi_ECoutCut->Write();

    ECtot_VS_P_ECoutCut->GetXaxis()->SetTitle("Momentum (GeV)");
    ECtot_VS_P_ECoutCut->GetYaxis()->SetTitle("EC total energy");
    ECtot_VS_P_ECoutCut->Write();
    
    ECtotP_VS_P_ECoutCut->GetXaxis()->SetTitle("Momentum (GeV/c)");
    ECtotP_VS_P_ECoutCut->GetYaxis()->SetTitle("EC_{total}/Mom.");
    ECtotP_VS_P_ECoutCut->Write();

    ECtotMinusECin_ECoutCut->GetXaxis()->SetTitle("EC_{total} - EC_{in} (GeV)");
    ECtotMinusECin_ECoutCut->GetYaxis()->SetTitle("Counts");
    ECtotMinusECin_ECoutCut->Write();
    
    q2_AntiECoutCut->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    q2_AntiECoutCut->GetYaxis()->SetTitle("Counts");
    q2_AntiECoutCut->Write();
    
    elecZVert_AntiECoutCut->GetXaxis()->SetTitle("e^{-} Z vertex (cm)");
    elecZVert_AntiECoutCut->GetYaxis()->SetTitle("Counts");
    elecZVert_AntiECoutCut->Write();
    
    Beta_VS_Momentum_AntiECoutCut->GetXaxis()->SetTitle("Momentum (GeV/c)");
    Beta_VS_Momentum_AntiECoutCut->GetYaxis()->SetTitle("#beta");
    Beta_VS_Momentum_AntiECoutCut->Write();
    
    Theta_VS_Phi_AntiECoutCut->GetXaxis()->SetTitle("#theta (deg.)");
    Theta_VS_Phi_AntiECoutCut->GetYaxis()->SetTitle("#phi (deg.)");
    Theta_VS_Phi_AntiECoutCut->Write();
    
    ECtot_VS_P_AntiECoutCut->GetXaxis()->SetTitle("Momentum (GeV)");
    ECtot_VS_P_AntiECoutCut->GetYaxis()->SetTitle("EC total energy");
    ECtot_VS_P_AntiECoutCut->Write();
    
    ECtotP_VS_P_AntiECoutCut->GetXaxis()->SetTitle("Momentum (GeV/c)");
    ECtotP_VS_P_AntiECoutCut->GetYaxis()->SetTitle("EC_{total}/Mom.");
    ECtotP_VS_P_AntiECoutCut->Write();
    
    ECtotMinusECin_AntiECoutCut->GetXaxis()->SetTitle("EC_{total} - EC_{in} (GeV)");
    ECtotMinusECin_AntiECoutCut->GetYaxis()->SetTitle("Counts");
    ECtotMinusECin_AntiECoutCut->Write();
    
    ECin_VS_ECout_ECfid->GetXaxis()->SetTitle("EC inner energy");
    ECin_VS_ECout_ECfid->GetYaxis()->SetTitle("EC outer energy");
    ECin_VS_ECout_ECfid->Write();

    ECin_VS_ECout_elecID_All->GetXaxis()->SetTitle("EC inner energy");
    ECin_VS_ECout_elecID_All->GetYaxis()->SetTitle("EC outer energy");
    ECin_VS_ECout_elecID_All->Write();

    for(i=0; i<MAX_SECTORS; i++){
        ECtotP_VS_P_Sector[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P_Sector[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P_Sector[i]->Write();
    }

    for(i=0; i<MAX_SECTORS; i++){
        ECinP_VS_ECoutP[i]->GetXaxis()->SetTitle("EC_{in}/Mom.");
        ECinP_VS_ECoutP[i]->GetYaxis()->SetTitle("EC_{out}/Mom.");
        ECinP_VS_ECoutP[i]->Write();
        ECinP_VS_ECoutP_cut[i]->GetXaxis()->SetTitle("EC_{in}/Mom.");
        ECinP_VS_ECoutP_cut[i]->GetYaxis()->SetTitle("EC_{out}/Mom.");
        ECinP_VS_ECoutP_cut[i]->Write();
    }

    for(i=0; i<5; i++) {
        ECinP_VS_ECoutP_Range[i]->GetXaxis()->SetTitle("EC_{in}/Mom.");
        ECinP_VS_ECoutP_Range[i]->GetYaxis()->SetTitle("EC_{out}/Mom.");
        ECinP_VS_ECoutP_Range[i]->Write();
    }

    for(i=0; i<MAX_SECTORS; i++){
        ECtotP_VS_P_ECPCut[i]->GetXaxis()->SetTitle("Momentum (GeV/c)");
        ECtotP_VS_P_ECPCut[i]->GetYaxis()->SetTitle("EC_{total}/Mom.");
        ECtotP_VS_P_ECPCut[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_Sector[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_Sector[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_Sector[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_ECoutCut[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_ECoutCut[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_ECoutCut[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_FidCut[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_FidCut[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_FidCut[i]->Write();
    }

    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_AntiFidCut[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_AntiFidCut[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_AntiFidCut[i]->Write();
    }

    // create a directory for photon id
    TDirectory *cdPhotID = out->mkdir("PhotonID");
    cdPhotID->cd();
    
    MomentumPhoton1->GetXaxis()->SetTitle("P (GeV/c)");
    MomentumPhoton1->GetYaxis()->SetTitle("Counts");
	MomentumPhoton1->Write();

    MomentumPhoton2->GetXaxis()->SetTitle("P (GeV/c)");
    MomentumPhoton2->GetYaxis()->SetTitle("Counts");
	MomentumPhoton2->Write();

    MomentumPhoton1_cut->GetXaxis()->SetTitle("P (GeV/c)");
    MomentumPhoton1_cut->GetYaxis()->SetTitle("Counts");
	MomentumPhoton1_cut->Write();

    MomentumPhoton2_cut->GetXaxis()->SetTitle("P (GeV/c)");
    MomentumPhoton2_cut->GetYaxis()->SetTitle("Counts");
	MomentumPhoton2_cut->Write();

    BetaPhoton1->GetXaxis()->SetTitle("#beta");
    BetaPhoton1->GetYaxis()->SetTitle("Counts");
	BetaPhoton1->Write();

    BetaPhoton2->GetXaxis()->SetTitle("#beta");
    BetaPhoton2->GetYaxis()->SetTitle("Counts");
	BetaPhoton2->Write();

    BetaPhoton1_cut->GetXaxis()->SetTitle("#beta");
    BetaPhoton1_cut->GetYaxis()->SetTitle("Counts");
	BetaPhoton1_cut->Write();

    BetaPhoton2_cut->GetXaxis()->SetTitle("#beta");
    BetaPhoton2_cut->GetYaxis()->SetTitle("Counts");
	BetaPhoton2_cut->Write();

    ECuPhoton1->GetXaxis()->SetTitle("EC U (cm)");
    ECuPhoton1->GetYaxis()->SetTitle("Counts");
	ECuPhoton1->Write();

    ECuPhoton2->GetXaxis()->SetTitle("EC U (cm)");
    ECuPhoton2->GetYaxis()->SetTitle("Counts");
	ECuPhoton2->Write();

    ECvPhoton1->GetXaxis()->SetTitle("EC V (cm)");
    ECvPhoton1->GetYaxis()->SetTitle("Counts");
	ECvPhoton1->Write();

    ECvPhoton2->GetXaxis()->SetTitle("EC V (cm)");
    ECvPhoton2->GetYaxis()->SetTitle("Counts");
	ECvPhoton2->Write();

    ECwPhoton1->GetXaxis()->SetTitle("EC W (cm)");
    ECwPhoton1->GetYaxis()->SetTitle("Counts");
	ECwPhoton1->Write();

    ECwPhoton2->GetXaxis()->SetTitle("EC W (cm)");
    ECwPhoton2->GetYaxis()->SetTitle("Counts");
	ECwPhoton2->Write();

    ECuPhoton1_cut->GetXaxis()->SetTitle("EC U (cm)");
    ECuPhoton1_cut->GetYaxis()->SetTitle("Counts");
	ECuPhoton1_cut->Write();

    ECuPhoton2_cut->GetXaxis()->SetTitle("EC U (cm)");
    ECuPhoton2_cut->GetYaxis()->SetTitle("Counts");
	ECuPhoton2_cut->Write();

    ECvPhoton1_cut->GetXaxis()->SetTitle("EC V (cm)");
    ECvPhoton1_cut->GetYaxis()->SetTitle("Counts");
	ECvPhoton1_cut->Write();

    ECvPhoton2_cut->GetXaxis()->SetTitle("EC V (cm)");
    ECvPhoton2_cut->GetYaxis()->SetTitle("Counts");
	ECvPhoton2_cut->Write();

    ECwPhoton1_cut->GetXaxis()->SetTitle("EC W (cm)");
    ECwPhoton1_cut->GetYaxis()->SetTitle("Counts");
	ECwPhoton1_cut->Write();

    ECwPhoton2_cut->GetXaxis()->SetTitle("EC W (cm)");
    ECwPhoton2_cut->GetYaxis()->SetTitle("Counts");
	ECwPhoton2_cut->Write();

    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_Sector_Photon1[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_Sector_Photon1[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_Sector_Photon1[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_FidCut_Photon1[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_FidCut_Photon1[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_FidCut_Photon1[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_AntiFidCut_Photon1[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_AntiFidCut_Photon1[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_AntiFidCut_Photon1[i]->Write();
    }

    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_Sector_Photon2[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_Sector_Photon2[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_Sector_Photon2[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_FidCut_Photon2[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_FidCut_Photon2[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_FidCut_Photon2[i]->Write();
    }
    
    for(i=0; i<MAX_SECTORS; i++){
        EC_XvsY_local_AntiFidCut_Photon2[i]->GetXaxis()->SetTitle("EC X_{local} (cm)");
        EC_XvsY_local_AntiFidCut_Photon2[i]->GetYaxis()->SetTitle("EC Y_{local} (cm)");
        EC_XvsY_local_AntiFidCut_Photon2[i]->Write();
    }
    
    ECtime_ECl_Start_Photon1->GetXaxis()->SetTitle("t_{EC} - t_{start} - l_{EC}/c (cm)");
    ECtime_ECl_Start_Photon1->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Start_Photon1->Write();
    
    ECtime_ECl_Start_Photon2->GetXaxis()->SetTitle("t_{EC} - t_{start} - l_{EC}/c (cm)");
    ECtime_ECl_Start_Photon2->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Start_Photon2->Write();
    
    ECtime_ECl_Photon1->GetXaxis()->SetTitle("t_{EC} - l_{EC}/c  (ns)");
    ECtime_ECl_Photon1->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Photon1->Write();

    ECtime_ECl_Photon2->GetXaxis()->SetTitle("t_{EC} - l_{EC}/c  (ns)");
    ECtime_ECl_Photon2->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Photon2->Write();

    ECtime_ECl_Photon1_cut->GetXaxis()->SetTitle("t_{EC} - l_{EC}/c (ns)");
    ECtime_ECl_Photon1_cut->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Photon1_cut->Write();

    ECtime_ECl_Photon2_cut->GetXaxis()->SetTitle("t_{EC} - l_{EC}/c  (ns)");
    ECtime_ECl_Photon2_cut->GetYaxis()->SetTitle("Counts");
    ECtime_ECl_Photon2_cut->Write();

    ECtimePhoton1->GetXaxis()->SetTitle("EC_{time} (ns)");
    ECtimePhoton1->GetYaxis()->SetTitle("Counts");
    ECtimePhoton1->Write();

    ECtimePhoton2->GetXaxis()->SetTitle("EC_{time} (ns)");
    ECtimePhoton2->GetYaxis()->SetTitle("Counts");
    ECtimePhoton2->Write();

    ECpathPhoton1->GetXaxis()->SetTitle("EC_{path} (cm)");
    ECpathPhoton1->GetYaxis()->SetTitle("Counts");
    ECpathPhoton1->Write();

    ECpathPhoton2->GetXaxis()->SetTitle("EC_{path} (cm)");
    ECpathPhoton2->GetYaxis()->SetTitle("Counts");
    ECpathPhoton2->Write();

    ECpathtimePhoton1->GetXaxis()->SetTitle("EC_{path}/c (ns)");
    ECpathtimePhoton1->GetYaxis()->SetTitle("Counts");
    ECpathtimePhoton1->Write();
    
    ECpathtimePhoton2->GetXaxis()->SetTitle("EC_{path}/c (ns)");
    ECpathtimePhoton2->GetYaxis()->SetTitle("Counts");
    ECpathtimePhoton2->Write();
    
    ECtotP_vs_P_Photon1->GetXaxis()->SetTitle("P (GeV/c)");
    ECtotP_vs_P_Photon1->GetYaxis()->SetTitle("EC_{tot} / P");
    ECtotP_vs_P_Photon1->Write();

    ECtotP_vs_P_Photon2->GetXaxis()->SetTitle("P (GeV/c)");
    ECtotP_vs_P_Photon2->GetYaxis()->SetTitle("EC_{tot} / P");
    ECtotP_vs_P_Photon2->Write();

    ECin_vs_ECout_Photon1->GetXaxis()->SetTitle("EC_{in}");
    ECin_vs_ECout_Photon1->GetYaxis()->SetTitle("EC_{out}");
    ECin_vs_ECout_Photon1->Write();

    ECin_vs_ECout_Photon2->GetXaxis()->SetTitle("EC_{in}");
    ECin_vs_ECout_Photon2->GetYaxis()->SetTitle("EC_{out}");
    ECin_vs_ECout_Photon2->Write();

    ECtotP_vs_P_InOutZeroCut_Photon1->GetXaxis()->SetTitle("P (GeV/c)");
    ECtotP_vs_P_InOutZeroCut_Photon1->GetYaxis()->SetTitle("EC_{tot} / P");
    ECtotP_vs_P_InOutZeroCut_Photon1->Write();
    
    ECtotP_vs_P_InOutZeroCut_Photon2->GetXaxis()->SetTitle("P (GeV/c)");
    ECtotP_vs_P_InOutZeroCut_Photon2->GetYaxis()->SetTitle("EC_{tot} / P");
    ECtotP_vs_P_InOutZeroCut_Photon2->Write();
    
    ECin_vs_ECout_InOutZeroCut_Photon1->GetXaxis()->SetTitle("EC_{in}");
    ECin_vs_ECout_InOutZeroCut_Photon1->GetYaxis()->SetTitle("EC_{out}");
    ECin_vs_ECout_InOutZeroCut_Photon1->Write();
    
    ECin_vs_ECout_InOutZeroCut_Photon2->GetXaxis()->SetTitle("EC_{in}");
    ECin_vs_ECout_InOutZeroCut_Photon2->GetYaxis()->SetTitle("EC_{out}");
    ECin_vs_ECout_InOutZeroCut_Photon2->Write();
    
    // create a directory for reconstructed particle cuts
    TDirectory *cdReconC = out->mkdir("ReconCuts");
    cdReconC->cd();

    Pi0Mass_PhotIDcuts->GetXaxis()->SetTitle("Reconstructed #pi^{0} mass (GeV/c^{2})");
    Pi0Mass_PhotIDcuts->GetYaxis()->SetTitle("Counts");
    Pi0Mass_PhotIDcuts->Write();

    OmegaMass_AllCuts->GetXaxis()->SetTitle("Reconstructed #omega mass (GeV/c^{2})");
    OmegaMass_AllCuts->GetYaxis()->SetTitle("Counts");
    OmegaMass_AllCuts->Write();

    mass2Pions_VS_massOmega_NC[0]->GetXaxis()->SetTitle("IM #pi^{+} #pi^{-} (GeV/c^{2})");
    mass2Pions_VS_massOmega_NC[0]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_NC[0]->Write();

    mass2Pions_VS_massOmega_NC[1]->GetXaxis()->SetTitle("IM #pi^{+} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_NC[1]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_NC[1]->Write();

    mass2Pions_VS_massOmega_NC[2]->GetXaxis()->SetTitle("IM #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_NC[2]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_NC[2]->Write();

    mass2Pions_VS_massOmega_NCX[0]->GetXaxis()->SetTitle("IM #pi^{+} #pi^{-} (GeV/c^{2})");
    mass2Pions_VS_massOmega_NCX[0]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_NCX[0]->Write();

    mass2Pions_VS_massOmega_NCX[1]->GetXaxis()->SetTitle("IM #pi^{+} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_NCX[1]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_NCX[1]->Write();

    mass2Pions_VS_massOmega_NCX[2]->GetXaxis()->SetTitle("IM #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_NCX[2]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_NCX[2]->Write();

    mass2Pions_VS_massOmega_EPC[0]->GetXaxis()->SetTitle("IM #pi^{+} #pi^{-} (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPC[0]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPC[0]->Write();

    mass2Pions_VS_massOmega_EPC[1]->GetXaxis()->SetTitle("IM #pi^{+} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPC[1]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPC[1]->Write();

    mass2Pions_VS_massOmega_EPC[2]->GetXaxis()->SetTitle("IM #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPC[2]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPC[2]->Write();

    mass2Pions_VS_massOmega_EPCX[0]->GetXaxis()->SetTitle("IM #pi^{+} #pi^{-} (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPCX[0]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_EPCX[0]->Write();

    mass2Pions_VS_massOmega_EPCX[1]->GetXaxis()->SetTitle("IM #pi^{+} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPCX[1]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_EPCX[1]->Write();

    mass2Pions_VS_massOmega_EPCX[2]->GetXaxis()->SetTitle("IM #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPCX[2]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_EPCX[2]->Write();

    mass2Pions_VS_massOmega_EPOC[0]->GetXaxis()->SetTitle("IM #pi^{+} #pi^{-} (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOC[0]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOC[0]->Write();

    mass2Pions_VS_massOmega_EPOC[1]->GetXaxis()->SetTitle("IM #pi^{+} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOC[1]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOC[1]->Write();

    mass2Pions_VS_massOmega_EPOC[2]->GetXaxis()->SetTitle("IM #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOC[2]->GetYaxis()->SetTitle("IM #pi^{+} #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOC[2]->Write();

    mass2Pions_VS_massOmega_EPOCX[0]->GetXaxis()->SetTitle("IM #pi^{+} #pi^{-} (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOCX[0]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_EPOCX[0]->Write();

    mass2Pions_VS_massOmega_EPOCX[1]->GetXaxis()->SetTitle("IM #pi^{+} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOCX[1]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_EPOCX[1]->Write();

    mass2Pions_VS_massOmega_EPOCX[2]->GetXaxis()->SetTitle("IM #pi^{-} #gamma #gamma (GeV/c^{2})");
    mass2Pions_VS_massOmega_EPOCX[2]->GetYaxis()->SetTitle("Counts");
    mass2Pions_VS_massOmega_EPOCX[2]->Write();
    
    out->Close();
}
