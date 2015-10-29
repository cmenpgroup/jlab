#ifndef ELECTRONID_H
#define ELECTRONID_H
#include <vector>
#include <string>

using namespace std;

class ElectronID
{
    vector<string> elecIDLabel;
    vector<double> RangeElecMom;
    vector<double> RangeECu;
    vector<double> RangeECv;
    vector<double> RangeECw;
    vector<double> RangeECin;
    vector<double> Range_dtECSC;
    vector<double> RangeCCnphe;
    

    // parameters to calculate the EC sampling fraction of total energy vs P
    double EC_SamplingFrac_C[6][5];
    double EC_SamplingFrac_Fe[6][5];
    double EC_SamplingFrac_Pb[6][5];
    
    // parameters to calculate the ECin/P vs ECout/P cut
    double Param_ECinP_ECoutP[6][4];
    
    bool cuts_ElecID;
    bool cuts_ElecID_Mom;
    bool cuts_ElecID_ECoverP;
    bool cuts_ElecID_ECin;
    bool cuts_ElecID_ECfid;
    bool cuts_ElecID_dtECSC;
    bool cuts_ElecID_CCnphe;
    bool cuts_ElecID_ECinP_ECoutP;
    
public:
    ElectronID();
    int Get_nElecID() {return elecIDLabel.size();};
    string Get_elecIDLabel(int num) {return elecIDLabel[num];};
    double Get_ElecMom_lo() {return RangeElecMom[0];};
    double Get_ElecMom_hi() {return RangeElecMom[1];};
    double Get_ElecECu_lo() {return RangeECu[0];};
    double Get_ElecECu_hi() {return RangeECu[1];};
    double Get_ElecECv_lo() {return RangeECv[0];};
    double Get_ElecECv_hi() {return RangeECv[1];};
    double Get_ElecECw_lo() {return RangeECw[0];};
    double Get_ElecECw_hi() {return RangeECw[1];};
    double Get_ElecECin_lo() {return RangeECin[0];};
    double Get_ElecECin_hi() {return RangeECin[1];};
    double Get_ElecCCnphe_lo() {return RangeCCnphe[0];};
    double Get_ElecCCnphe_hi() {return RangeCCnphe[1];};
    double Get_Elec_dtECSC_lo() {return Range_dtECSC[0];};
    double Get_Elec_dtECSC_hi() {return Range_dtECSC[1];};
    double Get_EC_SamplingFraction(int coeff, int sector, int targMass);
    double Get_ECinP_VS_ECoutP_Parameters(int coeff, int sector);
    
    bool Check_ElecMom(double mom);
    bool Check_ElecECu(double ecu);
    bool Check_ElecECv(double ecv);
    bool Check_ElecECw(double ecw);
    bool Check_ElecECin(double ecin);
    bool Check_Elec_dtECSC(double dt);
    bool Check_ElecCCnphe(double nphe);
    bool Check_ElecECoverP(double mom, double ectot, int sector, int targMass);
    bool Check_ElecECinP_VS_ECoutP(double mom, double ecin, double ecout, int sector);
    
    void InitCuts();
    void SetCut_ElecECinP_VS_ECoutP(double mom, double ecin, double ecout, int sector);
    bool GetCut_ElecECinP_VS_ECoutP() {return cuts_ElecID_ECinP_ECoutP;};
    void SetCut_ElecECoverP(double mom, double ectot, int sector, int targMass);
    bool GetCut_ElecECoverP() {return cuts_ElecID_ECoverP;};
    void SetCut_ElecCCnphe(double nphe);
    bool GetCut_ElecCCnphe() {return cuts_ElecID_CCnphe;};
    void SetCut_Elec_dtECSC(double dt);
    bool GetCut_Elec_dtECSC() {return cuts_ElecID_dtECSC;};
    void SetCut_ElecECin(double ecin);
    bool GetCut_ElecECin() {return cuts_ElecID_ECin;};
    void SetCut_ElecECfid(double ecu, double ecv, double ecw);
    bool GetCut_ElecECfid() {return cuts_ElecID_ECfid;};
    void SetCut_ElecMom(double mom);
    bool GetCut_ElecMom() {return cuts_ElecID_Mom;};
    void SetCut_Elec_All();
    bool GetCut_Elec_All() {return cuts_ElecID;};
    
    void Print_ElectronID();
};
#endif
