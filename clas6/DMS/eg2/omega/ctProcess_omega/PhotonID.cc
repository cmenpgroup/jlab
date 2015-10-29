#include <cmath>
#include <vector>
#include <string>
#include "PhotonID.h"
#include <iostream>
PhotonID::PhotonID()
{
    this->InitCuts(); // initialize the photon cuts
    
    PhotonIDLabel.push_back("No Cuts");
    PhotonIDLabel.push_back("Momentum");
    PhotonIDLabel.push_back("Beta");
    PhotonIDLabel.push_back("EC U-view");
    PhotonIDLabel.push_back("EC V-view");
    PhotonIDLabel.push_back("EC W-view");
    PhotonIDLabel.push_back("Timing");
    PhotonIDLabel.push_back("TOF Msq");
    PhotonIDLabel.push_back("ECinTimesECout");
    
//    RangePhotonMom.push_back(0.3); // Lower limit on photon momentum (in GeV)
    RangePhotonMom.push_back(0.15); // Lower limit on photon momentum (in GeV)
    RangePhotonMom.push_back(1000.0); // Upper limit on photon momentum (in GeV)

    RangePhotonBeta.push_back(0.95); // Lower limit on photon beta
    RangePhotonBeta.push_back(1.95); // Upper limit on photon beta

    RangeSCMsq1.push_back(-0.041537); // Lower limit on photon 1 TOF M^2
    RangeSCMsq1.push_back(0.031435); // Upper limit on photon 1 TOF M^2

    RangeSCMsq2.push_back(-0.044665); // Lower limit on photon 2 TOF M^2
    RangeSCMsq2.push_back(0.035156); // Upper limit on photon 2 TOF M^2
    
    RangeECu.push_back(40); // Lower limit on EC U-view (in cm)
    RangeECu.push_back(410); // Upper limit on EC U-view (in cm)

    RangeECv.push_back(0); // Lower limit on EC V-view (in cm)
    RangeECv.push_back(370); // Upper limit on EC V-view (in cm)

    RangeECw.push_back(0); // Lower limit on EC W-view (in cm)
    RangeECw.push_back(410); // Upper limit on EC W-view (in cm)

    RangePhotonECinTimesECout.push_back(0.0); // Lower limit on photon ECin times ECout (in GeV^2)
    RangePhotonECinTimesECout.push_back(1000.0); // Upper limit on photon ECin times ECout (in GeV^2)
    
    double dtCentroid, dtWidth, dtNsigmas, dtLo, dtHi;
    for (int i=1; i<=2; i++){
        switch (i){
            case 1: // difference between ECtime and ECpath/c (in ns) for photon 1
                dtCentroid = 54.0;
                dtWidth = 5.3;
                dtNsigmas = 3.0;
                dtLo = dtCentroid - dtNsigmas*dtWidth;
                dtHi = dtCentroid + dtNsigmas*dtWidth;
                RangeTiming1.push_back(dtLo);
                RangeTiming1.push_back(dtHi);
                break;
            case 2: // difference between ECtime and ECpath/c (in ns) for photon 2
                dtCentroid = 54.0;
                dtWidth = 5.3;
                dtNsigmas = 3.0;
                dtLo = dtCentroid - dtNsigmas*dtWidth;
                dtHi = dtCentroid + dtNsigmas*dtWidth;
                RangeTiming2.push_back(dtLo);
                RangeTiming2.push_back(dtHi);
                break;
            default:
                cout<<"PhotonID::PhotonTiming Initialization-> Wrong Photon number."<<endl;
                break;
        }
    }
}

// initialize the photon cuts
void PhotonID::InitCuts()
{
    cuts_photID1_mom = false;
    cuts_photID2_mom = false;
    cuts_photID_mom = false;
    cuts_photID1_beta = false;
    cuts_photID2_beta = false;
    cuts_photID_beta = false;
    cuts_photID1_scMsq = false;
    cuts_photID2_scMsq = false;
    cuts_photID_scMsq = false;
    cuts_photID1_fidu = false;
    cuts_photID2_fidu = false;
    cuts_photID1_fidv = false;
    cuts_photID2_fidv = false;
    cuts_photID1_fidw = false;
    cuts_photID2_fidw = false;
    cuts_photID1_fid = false;
    cuts_photID2_fid = false;
    cuts_photID_fid = false;
    cuts_photID1_time = false;
    cuts_photID2_time = false;
    cuts_photID_time = false;
    cuts_photID1_ECinTimesECout = false;
    cuts_photID2_ECinTimesECout = false;
    cuts_photID_ECinTimesECout = false;
    cuts_photID = false;
}

// check the cut on Photon momentum
bool PhotonID::Check_PhotonMom(double mom)
{
    bool ret = (mom >= this->Get_PhotonMom_lo() && mom < this->Get_PhotonMom_hi()) ? true : false;
    
    return ret;
}

// set the value of the momentum cut for individual photons
void PhotonID::SetCut_PhotonMom(double mom, int num)
{
    switch (num) {
        case 1: cuts_photID1_mom = this->Check_PhotonMom(mom); break;
        case 2: cuts_photID2_mom = this->Check_PhotonMom(mom); break;
        default:
            cout<<"PhotonID::SetCut_PhotonMom, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the momentum cut for individual photons
bool PhotonID::GetCut_PhotonMom(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_mom; break;
        case 2: ret = cuts_photID2_mom; break;
        default:
            cout<<"PhotonID::GetCut_PhotonMom, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the momentum cut for both photons
void PhotonID::SetCut_PhotonMom_All()
{
    cuts_photID_mom = (this->GetCut_PhotonMom(1) && this->GetCut_PhotonMom(2));
}

// check the cut on Photon beta
bool PhotonID::Check_PhotonBeta(double beta)
{
    bool ret = (beta >= this->Get_PhotonBeta_lo() && beta < this->Get_PhotonBeta_hi()) ? true : false;
    
    return ret;
}

// set the value of the beta cut for individual photons
void PhotonID::SetCut_PhotonBeta(double beta, int num)
{
    switch (num) {
        case 1: cuts_photID1_beta = this->Check_PhotonBeta(beta); break;
        case 2: cuts_photID2_beta = this->Check_PhotonBeta(beta); break;
        default:
            cout<<"PhotonID::SetCut_PhotonBeta, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the beta cut for individual photons
bool PhotonID::GetCut_PhotonBeta(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_beta; break;
        case 2: ret = cuts_photID2_beta; break;
        default:
            cout<<"PhotonID::GetCut_PhotonBeta, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the beta cut for both photons
void PhotonID::SetCut_PhotonBeta_All()
{
    cuts_photID_beta = (this->GetCut_PhotonBeta(1) && this->GetCut_PhotonBeta(2));
}

// check the cut on Photon EC U-view
bool PhotonID::Check_PhotonECu(double ecu)
{
    bool ret = (ecu >= this->Get_PhotonECu_lo() && ecu < this->Get_PhotonECu_hi()) ? true : false;
    
    return ret;
}

// set the value of the ECu cut for individual photons
void PhotonID::SetCut_PhotonECu(double ecu, int num)
{
    switch (num) {
        case 1: cuts_photID1_fidu = this->Check_PhotonECu(ecu); break;
        case 2: cuts_photID2_fidu = this->Check_PhotonECu(ecu); break;
        default:
            cout<<"PhotonID::SetCut_PhotonECu, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the ECu cut for individual photons
bool PhotonID::GetCut_PhotonECu(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_fidu; break;
        case 2: ret = cuts_photID2_fidu; break;
        default:
            cout<<"PhotonID::GetCut_PhotonECu, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// check the cut on Photon EC V-view
bool PhotonID::Check_PhotonECv(double ecv)
{
    bool ret = (ecv >= this->Get_PhotonECv_lo() && ecv < this->Get_PhotonECv_hi()) ? true : false;
    
    return ret;
}

// set the value of the ECv cut for individual photons
void PhotonID::SetCut_PhotonECv(double ecv, int num)
{
    switch (num) {
        case 1: cuts_photID1_fidv = this->Check_PhotonECv(ecv); break;
        case 2: cuts_photID2_fidv = this->Check_PhotonECv(ecv); break;
        default:
            cout<<"PhotonID::SetCut_PhotonECv, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the ECv cut for individual photons
bool PhotonID::GetCut_PhotonECv(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_fidv; break;
        case 2: ret = cuts_photID2_fidv; break;
        default:
            cout<<"PhotonID::GetCut_PhotonECv, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// check the cut on Photon EC W-view
bool PhotonID::Check_PhotonECw(double ecw)
{
    bool ret = (ecw >= this->Get_PhotonECw_lo() && ecw < this->Get_PhotonECw_hi()) ? true : false;
    
    return ret;
}

// set the value of the ECw cut for individual photons
void PhotonID::SetCut_PhotonECw(double ecw, int num)
{
    switch (num) {
        case 1: cuts_photID1_fidw = this->Check_PhotonECw(ecw); break;
        case 2: cuts_photID2_fidw = this->Check_PhotonECw(ecw); break;
        default:
            cout<<"PhotonID::SetCut_PhotonECw, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the ECw cut for individual photons
bool PhotonID::GetCut_PhotonECw(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_fidw; break;
        case 2: ret = cuts_photID2_fidw; break;
        default:
            cout<<"PhotonID::GetCut_PhotonECw, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the EC fid. cut (U, V, W) for individual photons
void PhotonID::SetCut_PhotonECfid(int num)
{
    bool fidU = this->GetCut_PhotonECu(num);
    bool fidV = this->GetCut_PhotonECv(num);
    bool fidW = this->GetCut_PhotonECw(num);
    
    switch (num) {
        case 1: cuts_photID1_fid = (fidU && fidV && fidW); break;
        case 2: cuts_photID2_fid = (fidU && fidV && fidW); break;
        default:
            cout<<"PhotonID::SetCut_PhotonECfid, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the ECfid cut for individual photons
bool PhotonID::GetCut_PhotonECfid(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_fid; break;
        case 2: ret = cuts_photID2_fid; break;
        default:
            cout<<"PhotonID::GetCut_PhotonECw, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the EC fid. cut (U, V, W) for individual photons
void PhotonID::SetCut_PhotonECfid_All()
{
    cuts_photID_fid = (this->GetCut_PhotonECfid(1) && this->GetCut_PhotonECfid(2));
}

// check the cut on Photon time difference between ECtime and ECpath/c
bool PhotonID::Check_PhotonTiming(double dt, int num)
{
    bool ret;
    switch (num){
        case 1:
            ret = (dt >= this->Get_PhotonTiming1_lo() && dt < this->Get_PhotonTiming1_hi()) ? true : false;
            break;
        case 2:
            ret = (dt >= this->Get_PhotonTiming2_lo() && dt < this->Get_PhotonTiming2_hi()) ? true : false;
            break;
        default:
            cout<<"PhotonID::Check_PhotonTiming -> Wrong Photon number."<<endl;
            ret = false;
            break;
    }
    return ret;
}

// set the value of the timing cut for individual photons
void PhotonID::SetCut_PhotonTiming(double dt, int num)
{
    switch (num) {
        case 1: cuts_photID1_time = this->Check_PhotonTiming(dt,num); break;
        case 2: cuts_photID2_time = this->Check_PhotonTiming(dt,num); break;
        default:
            cout<<"PhotonID::SetCut_PhotonECw, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the timing cut for individual photons
bool PhotonID::GetCut_PhotonTiming(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_time; break;
        case 2: ret = cuts_photID2_time; break;
        default:
            cout<<"PhotonID::GetCut_PhotonTiming, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the timing cut (U, V, W) for individual photons
void PhotonID::SetCut_PhotonTiming_All()
{
    cuts_photID_time = (this->GetCut_PhotonTiming(1) && this->GetCut_PhotonTiming(2));
}

// return the Photon TOF M^2 lower limit
double PhotonID::Get_PhotonSCMsq_lo(int num)
{
    double ret;
    switch (num){
        case 1: ret = RangeSCMsq1[0]; break;
        case 2: ret = RangeSCMsq2[0]; break;
        default:
            cout<<"PhotonID::Get_PhotonSCMsq_lo -> Wrong Photon number."<<endl;
            ret = false;
            break;
    }
    return ret;
}

// return the Photon TOF M^2 upper limit
double PhotonID::Get_PhotonSCMsq_hi(int num)
{
    double ret;
    switch (num){
        case 1: ret = RangeSCMsq1[1]; break;
        case 2: ret = RangeSCMsq2[1]; break;
        default:
            cout<<"PhotonID::Get_PhotonSCMsq_hi -> Wrong Photon number."<<endl;
            ret = false;
            break;
    }
    return ret;
}

// check the cut on Photon TOF M^2
bool PhotonID::Check_PhotonSCMsq(double scMsq, int num)
{
    bool ret = (scMsq >= this->Get_PhotonSCMsq_lo(num) && scMsq < this->Get_PhotonSCMsq_hi(num)) ? true : false;
    return ret;
}

// set the value of the TOF M^2 cut for individual photons
void PhotonID::SetCut_PhotonSCMsq(double scMsq, int num)
{
    switch (num) {
        case 1: cuts_photID1_scMsq = this->Check_PhotonSCMsq(scMsq,num); break;
        case 2: cuts_photID2_scMsq = this->Check_PhotonSCMsq(scMsq,num); break;
        default:
            cout<<"PhotonID::SetCut_PhotonSCMsq, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the TOF M^2 cut for individual photons
bool PhotonID::GetCut_PhotonSCMsq(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_scMsq; break;
        case 2: ret = cuts_photID2_scMsq; break;
        default:
            cout<<"PhotonID::GetCut_PhotonSCMsq, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the TOF M^2 for individual photons
void PhotonID::SetCut_PhotonSCMsq_All()
{
    cuts_photID_scMsq = (this->GetCut_PhotonSCMsq(1) && this->GetCut_PhotonSCMsq(2));
}

// check the cut on product of ECin and ECout
bool PhotonID::Check_PhotonECinTimesECout(double ECin, double ECout)
{
    bool ret = ((ECin*ECout > this->Get_PhotonECinTimesECout_lo()) && (ECin*ECout < this->Get_PhotonECinTimesECout_hi())) ? true : false;
    
    return ret;
}

// set the value of the ECin*ECout cut for individual photons
void PhotonID::SetCut_PhotonECinTimesECout(double ECin, double ECout, int num)
{
    switch (num) {
        case 1: cuts_photID1_ECinTimesECout = this->Check_PhotonECinTimesECout(ECin,ECout); break;
        case 2: cuts_photID2_ECinTimesECout = this->Check_PhotonECinTimesECout(ECin,ECout); break;
        default:
            cout<<"PhotonID::SetCut_PhotonECinTimesECout, Wrong photon number "<<num<<endl;
            break;
    }
}

// return the value of the ECin*ECout cut for individual photons
bool PhotonID::GetCut_PhotonECinTimesECout(int num)
{
    bool ret;
    switch (num) {
        case 1: ret = cuts_photID1_ECinTimesECout; break;
        case 2: ret = cuts_photID2_ECinTimesECout; break;
        default:
            cout<<"PhotonID::GetCut_PhotonECinTimesECout, Wrong photon number "<<num<<endl;
            break;
    }
    return ret;
}

// set the value of the ECin*ECout cut for both photons
void PhotonID::SetCut_PhotonECinTimesECout_All()
{
    cuts_photID_ECinTimesECout = (this->GetCut_PhotonECinTimesECout(1) && this->GetCut_PhotonECinTimesECout(2));
}

// set the value of the ECin*ECout cut for both photons
void PhotonID::SetCut_Photon_All()
{
//    cuts_photID = (this->GetCut_PhotonECinTimesECout_All() && this->GetCut_PhotonTiming_All() && this->GetCut_PhotonECfid_All() && this->GetCut_PhotonMom_All() && this->GetCut_PhotonBeta_All());
    cuts_photID = (this->GetCut_PhotonTiming_All() && this->GetCut_PhotonECfid_All() && this->GetCut_PhotonMom_All() && this->GetCut_PhotonBeta_All());

}

// print the cut information
void PhotonID::Print_PhotonID()
{
    int ii;
    cout<<"Photon ID Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nPhotonID();ii++){
        cout << this->Get_PhotonIDLabel(ii) << "\t";
        if (this->Get_PhotonIDLabel(ii).compare("Momentum")==0) {
            cout << "[" << this->Get_PhotonMom_lo() << "," << this->Get_PhotonMom_hi() << "] (GeV)" << endl;
        }else if (this->Get_PhotonIDLabel(ii).compare("Beta")==0) {
            cout << "[" << this->Get_PhotonBeta_lo() << "," << this->Get_PhotonBeta_hi() << "]" << endl;
        }else if (this->Get_PhotonIDLabel(ii).compare("EC U-view")==0) {
            cout << "[" << this->Get_PhotonECu_lo() << "," << this->Get_PhotonECu_hi() << "] (cm)" << endl;
        }else if (this->Get_PhotonIDLabel(ii).compare("EC V-view")==0) {
            cout << "[" << this->Get_PhotonECv_lo() << "," << this->Get_PhotonECv_hi() << "] (cm)" << endl;
        }else if (this->Get_PhotonIDLabel(ii).compare("EC W-view")==0) {
            cout << "[" << this->Get_PhotonECw_lo() << "," << this->Get_PhotonECw_hi() << "] (cm)" << endl;
        }else if (this->Get_PhotonIDLabel(ii).compare("Timing")==0) {
            cout<<endl;
            cout<<"Photon 1: ["<<this->Get_PhotonTiming1_lo()<<","<<this->Get_PhotonTiming1_hi()<<"] (ns)"<<endl;
            cout<<"Photon 2: ["<<this->Get_PhotonTiming2_lo()<<","<<this->Get_PhotonTiming2_hi()<<"] (ns)"<<endl;
        }else if (this->Get_PhotonIDLabel(ii).compare("TOF Msq")==0) {
            cout<<endl;
            cout<<"Photon 1: ["<<this->Get_PhotonSCMsq_lo(1)<<","<<this->Get_PhotonSCMsq_hi(1)<<"] (GeV/c^2)^2"<<endl;
            cout<<"Photon 2: ["<<this->Get_PhotonSCMsq_lo(2)<<","<<this->Get_PhotonSCMsq_hi(2)<<"] (GeV/c^2)^2"<<endl;
        }else if (this->Get_PhotonIDLabel(ii).compare("ECinTimesECout")==0) {
            cout << "[" << this->Get_PhotonECinTimesECout_lo() << "," << this->Get_PhotonECinTimesECout_hi() << "] (GeV^2)" << endl;
        }else{
            cout << endl;
        }
    }
    cout << endl;
}
