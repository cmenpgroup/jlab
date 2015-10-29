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
    vector<double> RangeSCMsq1;
    vector<double> RangeSCMsq2;
    vector<double> RangeECu;
    vector<double> RangeECv;
    vector<double> RangeECw;
    vector<double> RangeTiming1;
    vector<double> RangeTiming2;
    vector<double> RangePhotonECinTimesECout;
    
    bool cuts_photID1_mom;
    bool cuts_photID2_mom;
    bool cuts_photID_mom;
    bool cuts_photID1_beta;
    bool cuts_photID2_beta;
    bool cuts_photID_beta;
    bool cuts_photID1_scMsq;
    bool cuts_photID2_scMsq;
    bool cuts_photID_scMsq;
    bool cuts_photID1_fidu;
    bool cuts_photID2_fidu;
    bool cuts_photID1_fidv;
    bool cuts_photID2_fidv;
    bool cuts_photID1_fidw;
    bool cuts_photID2_fidw;
    bool cuts_photID1_fid;
    bool cuts_photID2_fid;
    bool cuts_photID_fid;
    bool cuts_photID1_time;
    bool cuts_photID2_time;
    bool cuts_photID_time;
    bool cuts_photID1_ECinTimesECout;
    bool cuts_photID2_ECinTimesECout;
    bool cuts_photID_ECinTimesECout;
    bool cuts_photID;
    
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
    double Get_PhotonSCMsq_lo(int num);
    double Get_PhotonSCMsq_hi(int num);
    
    bool Check_PhotonMom(double mom);
    bool Check_PhotonBeta(double beta);
    bool Check_PhotonECu(double ecu);
    bool Check_PhotonECv(double ecv);
    bool Check_PhotonECw(double ecw);
    bool Check_PhotonTiming(double dt, int num);
    bool Check_PhotonSCMsq(double scMsq, int num);
    bool Check_PhotonECinTimesECout(double ECin, double ECout);
    
    void InitCuts();
    void SetCut_Photon_All();
    bool GetCut_Photon_All() {return cuts_photID;};
    void SetCut_PhotonECinTimesECout_All();
    bool GetCut_PhotonECinTimesECout_All() {return cuts_photID_ECinTimesECout;};
    bool GetCut_PhotonECinTimesECout(int num);
    void SetCut_PhotonECinTimesECout(double ECin, double ECout, int num);
    void SetCut_PhotonTiming_All();
    bool GetCut_PhotonTiming_All() {return cuts_photID_time;};
    bool GetCut_PhotonTiming(int num);
    void SetCut_PhotonTiming(double dt, int num);
    void SetCut_PhotonSCMsq_All();
    bool GetCut_PhotonSCMsq_All() {return cuts_photID_scMsq;};
    bool GetCut_PhotonSCMsq(int num);
    void SetCut_PhotonSCMsq(double scMsq, int num);
    void SetCut_PhotonECfid_All();
    bool GetCut_PhotonECfid_All() {return cuts_photID_fid;};
    bool GetCut_PhotonECfid(int num);
    void SetCut_PhotonECfid(int num);
    bool GetCut_PhotonECw(int num);
    void SetCut_PhotonECw(double ecw, int num);
    bool GetCut_PhotonECv(int num);
    void SetCut_PhotonECv(double ecu, int num);
    bool GetCut_PhotonECu(int num);
    void SetCut_PhotonECu(double ecu, int num);
    void SetCut_PhotonBeta_All();
    bool GetCut_PhotonBeta_All() {return cuts_photID_beta;};
    bool GetCut_PhotonBeta(int num);
    void SetCut_PhotonBeta(double beta, int num);
    void SetCut_PhotonMom_All();
    bool GetCut_PhotonMom_All() {return cuts_photID_mom;};
    bool GetCut_PhotonMom(int num);
    void SetCut_PhotonMom(double mom, int num);
    
    void Print_PhotonID();
};
#endif
