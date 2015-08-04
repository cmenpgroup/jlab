#ifndef CHARGEDPIONID_H
#define CHARGEDPIONID_H
#include <vector>
#include <string>

using namespace std;

class ChargedPionID
{
    vector<string> ChargedPionIDLabel;
    vector<double> RangeChargedPionSCMassSq;
    vector<double> RangeChargedPionDiffBeta;
    
public:
    ChargedPionID();
    int Get_nChargedPionID() {return ChargedPionIDLabel.size();};
    string Get_ChargedPionIDLabel(int num) {return ChargedPionIDLabel[num];};
    double Get_ChargedPionSCMassSq_lo() {return RangeChargedPionSCMassSq[0];};
    double Get_ChargedPionSCMassSq_hi() {return RangeChargedPionSCMassSq[1];};
    double Get_ChargedPionDiffBeta_lo() {return RangeChargedPionDiffBeta[0];};
    double Get_ChargedPionDiffBeta_hi() {return RangeChargedPionDiffBeta[1];};
    
    bool Check_ChargedPionSCMassSq(double MassSq);
    bool Check_ChargedPionDiffBeta(double dBeta);
    
    void Print_ChargedPionID();
};
#endif
