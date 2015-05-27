#ifndef EG2CUTS_H
#define EG2CUTS_H
#include <vector>
#include <string>

using namespace std;

class EG2Cuts
{
    vector<string> CutsLabel;
    vector<double> RangeZdiff_ElecPim;
    vector<double> RangeZdiff_ElecPip;
    vector<double> RangeMassPi0;
    vector<double> RangeQSquared;
    vector<double> RangeOpAng_ElecPhoton;
    vector<double> RangeBetaPhoton;
    vector<double> RangeWcut;
    vector<double> RangeElecR;
    vector<double> RangeMassOmega;
    vector<double> RangeMassOmega_sb;
public:
    EG2Cuts();
    int Get_nCuts() {return CutsLabel.size();};
	string Get_CutsLabel(int num) {return CutsLabel[num];};
    double Get_Zdiff_ElecPim_lo() {return RangeZdiff_ElecPim[0];};
    double Get_Zdiff_ElecPim_hi() {return RangeZdiff_ElecPim[1];};
    double Get_Zdiff_ElecPip_lo() {return RangeZdiff_ElecPip[0];};
    double Get_Zdiff_ElecPip_hi() {return RangeZdiff_ElecPip[1];};
    double Get_MassPi0_lo() {return RangeMassPi0[0];};
    double Get_MassPi0_hi() {return RangeMassPi0[1];};
    double Get_QSquared_lo() {return RangeQSquared[0];};
    double Get_QSquared_hi() {return RangeQSquared[1];};
    double Get_BetaPhoton_lo() {return RangeBetaPhoton[0];};
    double Get_BetaPhoton_hi() {return RangeBetaPhoton[1];};
    double Get_Wcut_lo() {return RangeWcut[0];};
    double Get_Wcut_hi() {return RangeWcut[1];};
    double Get_OpAng_ElecPhoton_lo() {return RangeOpAng_ElecPhoton[0];};
    double Get_OpAng_ElecPhoton_hi() {return RangeOpAng_ElecPhoton[1];};
    double Get_ElectronR_lo() {return RangeElecR[0];};
    double Get_ElectronR_hi() {return RangeElecR[1];};
    double Get_MassOmega_lo() {return RangeMassOmega[0];};
    double Get_MassOmega_hi() {return RangeMassOmega[1];};
    double Get_MassOmega_sb_lo() {return RangeMassOmega_sb[0];};
    double Get_MassOmega_sb_hi() {return RangeMassOmega_sb[1];};
    bool Check_Zdiff_ElecPim(double zdiff);
    bool Check_Zdiff_ElecPip(double zdiff);
    bool Check_MassPi0(double mass);
    bool Check_QSquared(double Qsq);
    bool Check_OpAng_ElecPhoton(double OpAng);
    bool Check_Wcut(double W);
    bool Check_ElectronR(double vr);
    bool Check_BetaPhoton(double beta);
    bool Check_MassOmega(double mass);
    bool Check_MassOmega_sb(double mass);
    void Print_Cuts();
};
#endif
