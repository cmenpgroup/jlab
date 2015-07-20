#ifndef PHOTONID_H
#define PHOTONID_H
#include <vector>
#include <string>

using namespace std;

class PhotonID
{
    vector<string> PhotonIDLabel;
    vector<double> RangePhotonMom;
    vector<double> RangePhotonBeta;
    vector<double> RangeECu;
    vector<double> RangeECv;
    vector<double> RangeECw;
    vector<double> RangeTiming1;
    vector<double> RangeTiming2;
    vector<double> RangePhotonECinTimesECout;
    
public:
    PhotonID();
    int Get_nPhotonID() {return PhotonIDLabel.size();};
    string Get_PhotonIDLabel(int num) {return PhotonIDLabel[num];};
    double Get_PhotonMom_lo() {return RangePhotonMom[0];};
    double Get_PhotonMom_hi() {return RangePhotonMom[1];};
    double Get_PhotonBeta_lo() {return RangePhotonBeta[0];};
    double Get_PhotonBeta_hi() {return RangePhotonBeta[1];};
    double Get_PhotonECu_lo() {return RangeECu[0];};
    double Get_PhotonECu_hi() {return RangeECu[1];};
    double Get_PhotonECv_lo() {return RangeECv[0];};
    double Get_PhotonECv_hi() {return RangeECv[1];};
    double Get_PhotonECw_lo() {return RangeECw[0];};
    double Get_PhotonECw_hi() {return RangeECw[1];};
    double Get_PhotonTiming1_lo() {return RangeTiming1[0];};
    double Get_PhotonTiming1_hi() {return RangeTiming1[1];};
    double Get_PhotonTiming2_lo() {return RangeTiming2[0];};
    double Get_PhotonTiming2_hi() {return RangeTiming2[1];};
    double Get_PhotonECinTimesECout_lo() {return RangePhotonECinTimesECout[0];};
    double Get_PhotonECinTimesECout_hi() {return RangePhotonECinTimesECout[1];};
    
    bool Check_PhotonMom(double mom);
    bool Check_PhotonBeta(double beta);
    bool Check_PhotonECu(double ecu);
    bool Check_PhotonECv(double ecv);
    bool Check_PhotonECw(double ecw);
    bool Check_PhotonTiming(double dt, int num);
    bool Check_PhotonECinTimesECout(double ECin, double ECout);
    
    void Print_PhotonID();
};
#endif
