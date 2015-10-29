#include <cmath>
#include <vector>
#include <string>
#include "ElectronID.h"
#include <iostream>
ElectronID::ElectronID()
{
    elecIDLabel.push_back("No Cuts");
    elecIDLabel.push_back("Momentum");
    elecIDLabel.push_back("EC U-view");
    elecIDLabel.push_back("EC V-view");
    elecIDLabel.push_back("EC W-view");
    elecIDLabel.push_back("ECtot/P VS P");
    elecIDLabel.push_back("ECin");
    elecIDLabel.push_back("CC Nphe");
    elecIDLabel.push_back("dt(EC-SC)");
    
    RangeElecMom.push_back(0.64); // Lower limit on e- momentum (in GeV)
    RangeElecMom.push_back(1000.0); // Upper limit on e- momentum (in GeV)
    
    RangeECu.push_back(40); // Lower limit on EC U-view (in cm)
    RangeECu.push_back(400); // Upper limit on EC U-view (in cm)

    RangeECv.push_back(0); // Lower limit on EC V-view (in cm)
    RangeECv.push_back(360); // Upper limit on EC V-view (in cm)

    RangeECw.push_back(0); // Lower limit on EC W-view (in cm)
    RangeECw.push_back(390); // Upper limit on EC W-view (in cm)

    RangeECin.push_back(0.06); // Lower limit on EC inner energy (in GeV)
    RangeECin.push_back(10.0); // Upper limit on EC inner energy (in GeV)
    
    RangeCCnphe.push_back(28); // Lower limit on CC num. photo-electrons
    RangeCCnphe.push_back(9999); // Upper limit on CC num. photo-electrons

    double dtCentroid = 0.0;
    double dtWidth = 0.6;
    double dtNsigmas = 3.0;
    double dtLo = dtCentroid - dtNsigmas*dtWidth;
    double dtHi = dtCentroid + dtNsigmas*dtWidth;
    Range_dtECSC.push_back(dtLo); // Lower limit on time difference between EC and SC (in ns)
    Range_dtECSC.push_back(dtHi); // Upper limit on time difference between EC and SC (in ns)

    EC_SamplingFrac_C[0][0] = 0.226726; EC_SamplingFrac_C[0][1] = 0.0379557; EC_SamplingFrac_C[0][2] = -0.00855326; EC_SamplingFrac_C[0][3] = 7.27022e-09; EC_SamplingFrac_C[0][4] = 0.0370079;
    EC_SamplingFrac_C[1][0] = 0.222333; EC_SamplingFrac_C[1][1] = 0.0581705; EC_SamplingFrac_C[1][2] = -0.0131283; EC_SamplingFrac_C[1][3] = 3.12094e-12; EC_SamplingFrac_C[1][4] = 0.0413565;
    EC_SamplingFrac_C[2][0] = 0.245212; EC_SamplingFrac_C[2][1] = 0.0213835; EC_SamplingFrac_C[2][2] = -0.00277372; EC_SamplingFrac_C[2][3] = 8.27916e-08; EC_SamplingFrac_C[2][4] = 0.0426498;
    EC_SamplingFrac_C[3][0] = 0.238399; EC_SamplingFrac_C[3][1] = 0.0301926; EC_SamplingFrac_C[3][2] = -0.00720393; EC_SamplingFrac_C[3][3] = -3.81029e-09; EC_SamplingFrac_C[3][4] = 0.0309331;
    EC_SamplingFrac_C[4][0] = 0.241834; EC_SamplingFrac_C[4][1] = 0.0442975; EC_SamplingFrac_C[4][2] = -0.0105584; EC_SamplingFrac_C[4][3] = 9.74651e-09; EC_SamplingFrac_C[4][4] = 0.0303602;
    EC_SamplingFrac_C[5][0] = 0.245868; EC_SamplingFrac_C[5][1] = 0.0545128; EC_SamplingFrac_C[5][2] = -0.0149168; EC_SamplingFrac_C[5][3] = 1.43097e-08; EC_SamplingFrac_C[5][4] = 0.0483305;

/*    EC_SamplingFrac_C[0][0] = 2.52E-1; EC_SamplingFrac_C[0][1] = 1.22E-2; EC_SamplingFrac_C[0][2] = -7.94E-3; EC_SamplingFrac_C[0][3] = 9.55E-3; EC_SamplingFrac_C[0][4] = 3.41E-2;
    EC_SamplingFrac_C[1][0] = 2.78E-1; EC_SamplingFrac_C[1][1] = 1.87E-2; EC_SamplingFrac_C[1][2] = -2.38E-3; EC_SamplingFrac_C[1][3] = 1.399E-2; EC_SamplingFrac_C[1][4] = 3.75E-2;
    EC_SamplingFrac_C[2][0] = 2.62E-1; EC_SamplingFrac_C[2][1] = 2.31E-2; EC_SamplingFrac_C[2][2] = -3.54E-3; EC_SamplingFrac_C[2][3] = 9.32E-3; EC_SamplingFrac_C[2][4] = 2.90E-2;
    EC_SamplingFrac_C[3][0] = 2.51E-1; EC_SamplingFrac_C[3][1] = 2.01E-2; EC_SamplingFrac_C[3][2] = -3.32E-3; EC_SamplingFrac_C[3][3] = 8.21E-3; EC_SamplingFrac_C[3][4] = 2.99E-2;
    EC_SamplingFrac_C[4][0] = 2.63E-1; EC_SamplingFrac_C[4][1] = 9.55E-2; EC_SamplingFrac_C[4][2] = -1.02E-3; EC_SamplingFrac_C[4][3] = 2.25E-2; EC_SamplingFrac_C[4][4] = 3.06E-2;
    EC_SamplingFrac_C[5][0] = 2.55E-1; EC_SamplingFrac_C[5][1] = 2.32E-2; EC_SamplingFrac_C[5][2] = -3.05E-3; EC_SamplingFrac_C[5][3] = 1.17E-2; EC_SamplingFrac_C[5][4] = 3.64E-2;
*/
    EC_SamplingFrac_Fe[0][0] = 2.22E-1; EC_SamplingFrac_Fe[0][1] = 2.23E-2; EC_SamplingFrac_Fe[0][2] = -2.41E-3; EC_SamplingFrac_Fe[0][3] = 9.23E-3; EC_SamplingFrac_Fe[0][4] = 2.98E-2;
    EC_SamplingFrac_Fe[1][0] = 2.34E-1; EC_SamplingFrac_Fe[1][1] = 1.95E-2; EC_SamplingFrac_Fe[1][2] = -2.08E-3; EC_SamplingFrac_Fe[1][3] = 8.66E-3; EC_SamplingFrac_Fe[1][4] = 3.09E-2;
    EC_SamplingFrac_Fe[2][0] = 2.52E-1; EC_SamplingFrac_Fe[2][1] = 2.42E-2; EC_SamplingFrac_Fe[2][2] = -3.39E-3; EC_SamplingFrac_Fe[2][3] = 1.08E-2; EC_SamplingFrac_Fe[2][4] = 2.64E-2;
    EC_SamplingFrac_Fe[3][0] = 2.51E-1; EC_SamplingFrac_Fe[3][1] = 2.08E-2; EC_SamplingFrac_Fe[3][2] = -3.27E-3; EC_SamplingFrac_Fe[3][3] = 7.22E-3; EC_SamplingFrac_Fe[3][4] = 2.98E-2;
    EC_SamplingFrac_Fe[4][0] = 2.72E-1; EC_SamplingFrac_Fe[4][1] = 1.18E-2; EC_SamplingFrac_Fe[4][2] = -1.87E-3; EC_SamplingFrac_Fe[4][3] = 1.84E-2; EC_SamplingFrac_Fe[4][4] = 3.48E-2;
    EC_SamplingFrac_Fe[5][0] = 2.52E-1; EC_SamplingFrac_Fe[5][1] = 2.28E-2; EC_SamplingFrac_Fe[5][2] = -3.11E-3; EC_SamplingFrac_Fe[5][3] = 4.11E-3; EC_SamplingFrac_Fe[5][4] = 3.55E-2;

    EC_SamplingFrac_Pb[0][0] = 2.53E-1; EC_SamplingFrac_Pb[0][1] = 1.38E-2; EC_SamplingFrac_Pb[0][2] = -1.40E-3; EC_SamplingFrac_Pb[0][3] = 7.67E-3; EC_SamplingFrac_Pb[0][4] = 3.54E-2;
    EC_SamplingFrac_Pb[1][0] = 2.49E-1; EC_SamplingFrac_Pb[1][1] = 1.47E-2; EC_SamplingFrac_Pb[1][2] = -1.49E-3; EC_SamplingFrac_Pb[1][3] = 7.53E-3; EC_SamplingFrac_Pb[1][4] = 3.38E-2;
    EC_SamplingFrac_Pb[2][0] = 2.54E-1; EC_SamplingFrac_Pb[2][1] = 2.26E-2; EC_SamplingFrac_Pb[2][2] = -3.05E-3; EC_SamplingFrac_Pb[2][3] = 8.13E-3; EC_SamplingFrac_Pb[2][4] = 2.77E-2;
    EC_SamplingFrac_Pb[3][0] = 2.55E-1; EC_SamplingFrac_Pb[3][1] = 1.90E-2; EC_SamplingFrac_Pb[3][2] = -3.05E-3; EC_SamplingFrac_Pb[3][3] = 7.20E-3; EC_SamplingFrac_Pb[3][4] = 3.04E-2;
    EC_SamplingFrac_Pb[4][0] = 2.76E-1; EC_SamplingFrac_Pb[4][1] = 1.11E-2; EC_SamplingFrac_Pb[4][2] = -1.76E-3; EC_SamplingFrac_Pb[4][3] = 1.81E-2; EC_SamplingFrac_Pb[4][4] = 3.53E-2;
    EC_SamplingFrac_Pb[5][0] = 2.62E-1; EC_SamplingFrac_Pb[5][1] = 1.92E-2; EC_SamplingFrac_Pb[5][2] = -2.62E-3; EC_SamplingFrac_Pb[5][3] = 1.99E-3; EC_SamplingFrac_Pb[5][4] = 3.76E-2;
    
    // parameters for ECin/P vs ECout/P cut
    // Sector 1
    Param_ECinP_ECoutP[0][0] = 0.192702;
    Param_ECinP_ECoutP[0][1] = 0.0420221;
    Param_ECinP_ECoutP[0][2] = -0.586568;
    Param_ECinP_ECoutP[0][3] = -0.104084;

    // Sector 2
    Param_ECinP_ECoutP[1][0] = 0.2118;
    Param_ECinP_ECoutP[1][1] = 0.0478416;
    Param_ECinP_ECoutP[1][2] = -0.631934;
    Param_ECinP_ECoutP[1][3] = -0.120601;

    // Sector 3
    Param_ECinP_ECoutP[2][0] = 0.179616;
    Param_ECinP_ECoutP[2][1] = 0.052625;
    Param_ECinP_ECoutP[2][2] = -0.515158;
    Param_ECinP_ECoutP[2][3] = -0.149614;

    // Sector 4
    Param_ECinP_ECoutP[3][0] = 0.198436;
    Param_ECinP_ECoutP[3][1] = 0.0387745;
    Param_ECinP_ECoutP[3][2] = -0.621269;
    Param_ECinP_ECoutP[3][3] = -0.0961029;
    
    // Sector 5
    Param_ECinP_ECoutP[4][0] = 0.231849;
    Param_ECinP_ECoutP[4][1] = 0.0349574;
    Param_ECinP_ECoutP[4][2] = -0.734145;
    Param_ECinP_ECoutP[4][3] = -0.0709758;
    
    // Sector 6
    Param_ECinP_ECoutP[5][0] = 0.177784;
    Param_ECinP_ECoutP[5][1] = 0.0620507;
    Param_ECinP_ECoutP[5][2] = -0.451267;
    Param_ECinP_ECoutP[5][3] = -0.17344;
    
    this->InitCuts();
}

void ElectronID::InitCuts()
{
    cuts_ElecID = false;
    cuts_ElecID_Mom =false;
    cuts_ElecID_ECoverP = false;
    cuts_ElecID_ECin = false;
    cuts_ElecID_ECfid = false;
    cuts_ElecID_dtECSC = false;
    cuts_ElecID_CCnphe = false;
    cuts_ElecID_ECinP_ECoutP = false;
}

// check the cut on electron momentum
bool ElectronID::Check_ElecMom(double mom)
{
    bool ret = (mom >= this->Get_ElecMom_lo() && mom < this->Get_ElecMom_hi()) ? true : false;
    
    return ret;
}

// set the momentum cut
void ElectronID::SetCut_ElecMom(double mom)
{
    cuts_ElecID_Mom = this->Check_ElecMom(mom);
}

// check the cut on electron EC U-view
bool ElectronID::Check_ElecECu(double ecu)
{
    bool ret = (ecu >= this->Get_ElecECu_lo() && ecu < this->Get_ElecECu_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron EC V-view
bool ElectronID::Check_ElecECv(double ecv)
{
    bool ret = (ecv >= this->Get_ElecECv_lo() && ecv < this->Get_ElecECv_hi()) ? true : false;
    
    return ret;
}

// check the cut on electron EC W-view
bool ElectronID::Check_ElecECw(double ecw)
{
    bool ret = (ecw >= this->Get_ElecECw_lo() && ecw < this->Get_ElecECw_hi()) ? true : false;
    
    return ret;
}

// set the EC fid. cut
void ElectronID::SetCut_ElecECfid(double ecu, double ecv, double ecw)
{
    cuts_ElecID_ECfid = (this->Check_ElecECu(ecu) && this->Check_ElecECv(ecv) && this->Check_ElecECw(ecw));
}

// check the cut on electron EC inner energy
bool ElectronID::Check_ElecECin(double ecin)
{
    bool ret = (ecin >= this->Get_ElecECin_lo() && ecin < this->Get_ElecECin_hi()) ? true : false;
    
    return ret;
}

// set the ECin cut
void ElectronID::SetCut_ElecECin(double ecin)
{
    cuts_ElecID_ECin = this->Check_ElecECin(ecin);
}

// check the cut on electron time difference between EC and SC
bool ElectronID::Check_Elec_dtECSC(double dt)
{
    bool ret = (dt >= this->Get_Elec_dtECSC_lo() && dt < this->Get_Elec_dtECSC_hi()) ? true : false;
    
    return ret;
}

// set the electron timing cut
void ElectronID::SetCut_Elec_dtECSC(double dt)
{
    cuts_ElecID_dtECSC = this->Check_Elec_dtECSC(dt);
}

// check the cut on electron CC num. of photo-electrons
bool ElectronID::Check_ElecCCnphe(double nphe)
{
    bool ret = (nphe >= this->Get_ElecCCnphe_lo() && nphe < this->Get_ElecCCnphe_hi()) ? true : false;
    
    return ret;
}

// set the electron CC cut
void ElectronID::SetCut_ElecCCnphe(double nphe)
{
    cuts_ElecID_CCnphe = this->Check_ElecCCnphe(nphe);
}

double ElectronID::Get_EC_SamplingFraction(int coeff, int sector, int targMass)
{
    double ret = 0.0;
    
    if(sector>=1 && sector<=6){ //check that the sector is between 1 and 6
        if(coeff>=0 && coeff<5){
            switch (targMass){
                case 12: ret = this->EC_SamplingFrac_C[sector-1][coeff]; break;
                case 56: ret = this->EC_SamplingFrac_Fe[sector-1][coeff]; break;
                case 208: ret = this->EC_SamplingFrac_Pb[sector-1][coeff]; break;
                default:
                    cout<<"ElectronID::Get_EC_SamplingFraction: Target Mass "<< targMass <<" is unknown."<<endl;
                    ret = 0.0;
                    break;
            }
        }
        else{
            cout<<"ElectronID::Get_EC_SamplingFraction: Coefficient "<<coeff<<" is out of range."<<endl;
        }
    }
    else{
        cout<<"ElectronID::Get_EC_SamplingFraction: Sector "<<sector<<" is out of range."<<endl;
    }
    return ret;
}

// check the cut on electron EC inner energy
bool ElectronID::Check_ElecECoverP(double mom, double ectot, int sector, int targMass)
{
    bool ret = false; // initialize to false
    
    double a = this->Get_EC_SamplingFraction(0,sector,targMass);
    double b = this->Get_EC_SamplingFraction(1,sector,targMass);
    double c = this->Get_EC_SamplingFraction(2,sector,targMass);
    double d = this->Get_EC_SamplingFraction(3,sector,targMass);
    double f = this->Get_EC_SamplingFraction(4,sector,targMass);

    double centroid = a + b*mom + c*mom*mom;
    double sigma = sqrt(d*d + f*f/sqrt(mom));
    double Nsigma = 2.5;
    
    double diff = fabs(ectot/mom - centroid);
    
    ret = (diff < Nsigma*sigma) ? true : false;

    return ret;
}

// set the electron EC/P vs P cut
void ElectronID::SetCut_ElecECoverP(double mom, double ectot, int sector, int targMass)
{
    cuts_ElecID_ECoverP = this->Check_ElecECoverP(mom,ectot,sector,targMass);
}

double ElectronID::Get_ECinP_VS_ECoutP_Parameters(int coeff, int sector)
{
    double ret = 0.0;
    
    if(sector>=1 && sector<=6){ //check that the sector is between 1 and 6
        if(coeff>=0 && coeff<4){
            ret = this->Param_ECinP_ECoutP[sector-1][coeff];
        }
        else{
            cout<<"ElectronID::Get_ECinP_VS_ECout_Parameters: Coefficient "<<coeff<<" is out of range."<<endl;
        }
    }
    else{
        cout<<"ElectronID::Get_ECinP_VS_ECout_Parameters: Sector "<<sector<<" is out of range."<<endl;
    }
    return ret;
}

// check the cut on electron EC inner energy
bool ElectronID::Check_ElecECinP_VS_ECoutP(double mom, double ecin, double ecout, int sector)
{
    bool ret = false; // initialize to false
    
    double Nsigma = 2.5;
    
    double aMean = this->Get_ECinP_VS_ECoutP_Parameters(0,sector);
    double aSigma = this->Get_ECinP_VS_ECoutP_Parameters(1,sector);
    double bMean = this->Get_ECinP_VS_ECoutP_Parameters(2,sector);
    double bSigma = this->Get_ECinP_VS_ECoutP_Parameters(3,sector);
    
    double aAbove = aMean + Nsigma * aSigma;
    double bAbove = bMean + Nsigma * bSigma;
    double aBelow = aMean - Nsigma * aSigma;
    double bBelow = bMean - Nsigma * bSigma;
    
    double upper = ecout/mom - (bBelow*ecin/mom + aBelow);
    double lower = ecout/mom - (bAbove*ecin/mom + aAbove);
    
    ret = ((upper > 0) && (lower < 0)) ? true : false;
    
    return ret;
}

// set the electron ECin/P vs ECout/P cut
void ElectronID::SetCut_ElecECinP_VS_ECoutP(double mom, double ecin, double ecout, int sector)
{
    cuts_ElecID_ECinP_ECoutP = this->Check_ElecECinP_VS_ECoutP(mom,ecin,ecout,sector);
}

// set the total electron cut
void ElectronID::SetCut_Elec_All()
{
    cuts_ElecID = (this->GetCut_ElecECoverP() && this->GetCut_ElecCCnphe() && this->GetCut_Elec_dtECSC() && this->GetCut_ElecECin() && this->GetCut_ElecECfid() && this->GetCut_ElecMom());

}

// print the cut information
void ElectronID::Print_ElectronID()
{
    int ii;
    cout<<"Electron ID Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nElecID();ii++){
        cout << this->Get_elecIDLabel(ii) << "\t";
        if (this->Get_elecIDLabel(ii).compare("Momentum")==0) {
            cout << "[" << this->Get_ElecMom_lo() << "," << this->Get_ElecMom_hi() << "] (GeV)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("EC U-view")==0) {
            cout << "[" << this->Get_ElecECu_lo() << "," << this->Get_ElecECu_hi() << "] (cm)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("EC V-view")==0) {
            cout << "[" << this->Get_ElecECv_lo() << "," << this->Get_ElecECv_hi() << "] (cm)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("EC W-view")==0) {
            cout << "[" << this->Get_ElecECw_lo() << "," << this->Get_ElecECw_hi() << "] (cm)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("ECin")==0) {
            cout << "[" << this->Get_ElecECin_lo() << "," << this->Get_ElecECin_hi() << "] (GeV)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("CC Nphe")==0) {
            cout << "[" << this->Get_ElecCCnphe_lo() << "," << this->Get_ElecCCnphe_hi() << "]" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("dt(EC-SC)")==0) {
            cout << "[" << this->Get_Elec_dtECSC_lo() << "," << this->Get_Elec_dtECSC_hi() << "] (ns)" << endl;
        }else if (this->Get_elecIDLabel(ii).compare("ECtot/P VS P")==0) {
            cout << "Depends on the momentum dependent sampling fraction" << endl;
        }else{
            cout << endl;
        }
    }
    cout << endl;
}
