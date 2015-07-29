#include <cmath>
#include <vector>
#include <string>
#include "PhotonID.h"
#include <iostream>
PhotonID::PhotonID()
{
    PhotonIDLabel.push_back("No Cuts");
    PhotonIDLabel.push_back("Momentum");
    PhotonIDLabel.push_back("Beta");
    PhotonIDLabel.push_back("EC U-view");
    PhotonIDLabel.push_back("EC V-view");
    PhotonIDLabel.push_back("EC W-view");
    PhotonIDLabel.push_back("Timing");
    PhotonIDLabel.push_back("ECinTimesECout");
    
    RangePhotonMom.push_back(0.3); // Lower limit on photon momentum (in GeV)
    RangePhotonMom.push_back(1000.0); // Upper limit on photon momentum (in GeV)

    RangePhotonBeta.push_back(0.95); // Lower limit on photon beta
    RangePhotonBeta.push_back(1.95); // Upper limit on photon beta
    
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
                dtCentroid = -0.0882054;
                dtWidth = 0.640051;
                dtNsigmas = 3.0;
                dtLo = dtCentroid - dtNsigmas*dtWidth;
                dtHi = dtCentroid + dtNsigmas*dtWidth;
                RangeTiming1.push_back(dtLo);
                RangeTiming1.push_back(dtHi);
                break;
            case 2: // difference between ECtime and ECpath/c (in ns) for photon 2
                dtCentroid = -0.166546;
                dtWidth = 0.710022;
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

// check the cut on Photon momentum
bool PhotonID::Check_PhotonMom(double mom)
{
    bool ret = (mom >= this->Get_PhotonMom_lo() && mom < this->Get_PhotonMom_hi()) ? true : false;
    
    return ret;
}

// check the cut on Photon beta
bool PhotonID::Check_PhotonBeta(double beta)
{
    bool ret = (beta >= this->Get_PhotonBeta_lo() && beta < this->Get_PhotonBeta_hi()) ? true : false;
    
    return ret;
}

// check the cut on Photon EC U-view
bool PhotonID::Check_PhotonECu(double ecu)
{
    bool ret = (ecu >= this->Get_PhotonECu_lo() && ecu < this->Get_PhotonECu_hi()) ? true : false;
    
    return ret;
}

// check the cut on Photon EC V-view
bool PhotonID::Check_PhotonECv(double ecv)
{
    bool ret = (ecv >= this->Get_PhotonECv_lo() && ecv < this->Get_PhotonECv_hi()) ? true : false;
    
    return ret;
}

// check the cut on Photon EC W-view
bool PhotonID::Check_PhotonECw(double ecw)
{
    bool ret = (ecw >= this->Get_PhotonECw_lo() && ecw < this->Get_PhotonECw_hi()) ? true : false;
    
    return ret;
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

// check the cut on product of ECin and ECout
bool PhotonID::Check_PhotonECinTimesECout(double ECin, double ECout)
{
    bool ret = ((ECin*ECout > this->Get_PhotonECinTimesECout_lo()) && (ECin*ECout < this->Get_PhotonECinTimesECout_hi())) ? true : false;
    
    return ret;
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
        }else if (this->Get_PhotonIDLabel(ii).compare("ECinTimesECout")==0) {
            cout << "[" << this->Get_PhotonECinTimesECout_lo() << "," << this->Get_PhotonECinTimesECout_hi() << "] (GeV^2)" << endl;
        }else{
            cout << endl;
        }
    }
    cout << endl;
}
