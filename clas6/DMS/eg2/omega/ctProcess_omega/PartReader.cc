#include <cmath>
#include <vector>
#include <string>
#include "PartReader.h"
#include <iostream>
PartReader::PartReader(TTree *tree, string branchName)
{
    b_part = tree->GetBranch(branchName.c_str());
    l_sector = b_part->GetLeaf("Sector");
    l_charge = b_part->GetLeaf("Charge");
    l_pid = b_part->GetLeaf("Pid");
    l_beta = b_part->GetLeaf("Beta");
    l_px = b_part->GetLeaf("Px");
    l_py = b_part->GetLeaf("Py");
    l_pz = b_part->GetLeaf("Pz");
    l_mom = b_part->GetLeaf("Mom");
    l_mass2 = b_part->GetLeaf("Mass2");
    l_x = b_part->GetLeaf("X");
    l_y = b_part->GetLeaf("Y");
    l_z = b_part->GetLeaf("Z");
    l_ecx = b_part->GetLeaf("ECx");
    l_ecy = b_part->GetLeaf("ECy");
    l_ecz = b_part->GetLeaf("ECz");
    l_ecu = b_part->GetLeaf("ECu");
    l_ecv = b_part->GetLeaf("ECv");
    l_ecw = b_part->GetLeaf("ECw");
    l_ectot = b_part->GetLeaf("ECtot");
    l_ecin = b_part->GetLeaf("ECin");
    l_ecout = b_part->GetLeaf("ECout");
    l_ectime = b_part->GetLeaf("ECtime");
    l_ecpath = b_part->GetLeaf("ECpath");
    l_echit_m2 = b_part->GetLeaf("EChit_M2");
    l_echit_m3 = b_part->GetLeaf("EChit_M3");
    l_echit_m4 = b_part->GetLeaf("EChit_M4");
    l_chi2ec = b_part->GetLeaf("Chi2EC");
    l_sctime = b_part->GetLeaf("SCtime");
    l_scpath = b_part->GetLeaf("SCpath");
    l_ccnphe = b_part->GetLeaf("CCnphe");
    l_t = b_part->GetLeaf("T");
    l_xf = b_part->GetLeaf("Xf");
    l_mx2 = b_part->GetLeaf("Mx2");
    l_pt = b_part->GetLeaf("Pt");
    l_zh = b_part->GetLeaf("Zh");
    l_thetapq = b_part->GetLeaf("ThetaPQ");
    l_phipq = b_part->GetLeaf("PhiPQ");
    l_timecorr4 = b_part->GetLeaf("TimeCorr4");
}

void PartReader::ReadEntry(int num)
{
    b_part->GetEntry(num);
}

TLorentzVector PartReader::GetLorentzVector(double mass)
{
    TLorentzVector TLVect;
    TLVect.SetXYZM(this->Get_Px(),this->Get_Py(),this->Get_Pz(),mass);
    
    return TLVect;
}

TVector3 PartReader::GetVertex()
{
    TVector3 Vect;
    Vect.SetXYZ(this->Get_X(),this->Get_Y(),this->Get_Z());
    
    return Vect;
}

// Return the Time-of-Flight mass squared
double PartReader::Get_TOF_MassSquared()
{
    double ret = -99.0;
    double fBeta = this->Get_Beta();
    double fBetaSq = fBeta*fBeta;
    
    if(fBetaSq) ret = this->Get_Mom()*this->Get_Mom()*(1.0-fBetaSq)/fBetaSq;
    
    return ret;
}


// Return the expected particle beta given the mass and momentum
//
//          fMass = particle mass
//
double PartReader::Get_BetaFromMass(double fMass){
    
    double ret = -99.0;
    double fMom = this->Get_Mom();
    
    if(fMom) ret = 1.0/sqrt((fMass*fMass)/(fMom*fMom) + 1.0);

    return ret;
}

// Return the difference in beta between measured beta and that given the mass and momentum
//
//          fMass = particle mass
//
double PartReader::Get_BetaDifference(double fMass){
    
    double ret = -99.0;
    double fBetaFromMass = this->Get_BetaFromMass(fMass);
    
    if(fBetaFromMass!=-99.0) ret = this->Get_Beta() - fBetaFromMass;
    
    return ret;
}
