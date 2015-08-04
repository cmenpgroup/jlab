#include <cmath>
#include <vector>
#include <string>
#include "ChargedPionID.h"
#include <iostream>
ChargedPionID::ChargedPionID()
{
    ChargedPionIDLabel.push_back("No Cuts");
    ChargedPionIDLabel.push_back("SCMassSq");
    ChargedPionIDLabel.push_back("DiffBeta");

    RangeChargedPionDiffBeta.push_back(-0.035); // Lower limit on beta_TOF - beta_ideal
    RangeChargedPionDiffBeta.push_back(0.025); // Upper limit on beta_TOF - beta_ideal
    
    double Centroid, Width, Nsigmas, Lo, Hi;
    Centroid = 0.017956;
    Width = 0.005;
    Nsigmas = 3.0;
    Lo = Centroid - Nsigmas*Width;
    Hi = Centroid + Nsigmas*Width;
    RangeChargedPionSCMassSq.push_back(Lo);
    RangeChargedPionSCMassSq.push_back(Hi);

}

// check the cut on TOF mass squared
bool ChargedPionID::Check_ChargedPionSCMassSq(double MassSq)
{
    bool ret = (MassSq >= this->Get_ChargedPionSCMassSq_lo() && MassSq < this->Get_ChargedPionSCMassSq_hi()) ? true : false;
    
    return ret;
}

// check the cut on beta difference
bool ChargedPionID::Check_ChargedPionDiffBeta(double dBeta)
{
    bool ret = (dBeta >= this->Get_ChargedPionDiffBeta_lo() && dBeta < this->Get_ChargedPionDiffBeta_hi()) ? true : false;
    
    return ret;
}

// print the cut information
void ChargedPionID::Print_ChargedPionID()
{
    int ii;
    cout<<"ChargedPion ID Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nChargedPionID();ii++){
        cout << this->Get_ChargedPionIDLabel(ii) << "\t";
        if (this->Get_ChargedPionIDLabel(ii).compare("SCMassSq")==0) {
            cout << "[" << this->Get_ChargedPionSCMassSq_lo() << "," << this->Get_ChargedPionSCMassSq_hi() << "] (GeV)" << endl;
        }else if (this->Get_ChargedPionIDLabel(ii).compare("DiffBeta")==0) {
            cout << "[" << this->Get_ChargedPionDiffBeta_lo() << "," << this->Get_ChargedPionDiffBeta_hi() << "]" << endl;
        }else{
            cout << endl;
        }
    }
    cout << endl;
}
