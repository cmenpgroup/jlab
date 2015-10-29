#ifndef EG2CUTS_H
#define EG2CUTS_H
#include <vector>
#include <string>

using namespace std;

class EG2Cuts
{
    vector<string> CutsLabel;
    vector<double> RangeZDiff_ElecPim;
    vector<double> RangeZDiff_ElecPip;
    vector<double> RangeMassPi0;
    vector<double> RangeMassPipPim;
    vector<double> RangeQSquared;
    vector<double> RangeOpAng_ElecPhoton;
    vector<double> RangeBetaPhoton;
    vector<double> RangeWcut;
    vector<double> RangeElecR;
    vector<double> RangeMassOmega;
    vector<double> RangeMassOmega_sb;
    
    int topo_nelec;
    int topo_npim;
    int topo_npip;
    int topo_ngam;
    
    bool cuts_omega_MPi0;
    bool cuts_omega_MPipPim;
    bool cuts_omega_ZDiff_ElecPim;
    bool cuts_omega_ZDiff_ElecPip;
    bool cuts_omega_ZDiff;
    bool cuts_omega_Q2;
    bool cuts_omega_W;
    bool cuts_omega_OpAng_ElecPhot1;
    bool cuts_omega_OpAng_ElecPhot2;
    bool cuts_omega_OpAng_ElecPhot;
    bool cuts_omega_NumDetPart;
    bool cuts_omega_MPipPimPi0;
    bool cuts_omega_MPipPimPi0_sb;
    bool cuts_omega_All;
    
    bool cuts_omega_woMPi0;
    bool cuts_omega_woMPipPim;
    bool cuts_omega_woZDiff;
    bool cuts_omega_woQ2;
    bool cuts_omega_woOpAng_ElecPhot;
    bool cuts_omega_woNumDetPart;
    bool cuts_omega_woW;
    
public:
    EG2Cuts();
    int Get_nCuts() {return CutsLabel.size();};
	string Get_CutsLabel(int num) {return CutsLabel[num];};
    double Get_ZDiff_ElecPim_lo() {return RangeZDiff_ElecPim[0];};
    double Get_ZDiff_ElecPim_hi() {return RangeZDiff_ElecPim[1];};
    double Get_ZDiff_ElecPip_lo() {return RangeZDiff_ElecPip[0];};
    double Get_ZDiff_ElecPip_hi() {return RangeZDiff_ElecPip[1];};
    double Get_MassPi0_lo() {return RangeMassPi0[0];};
    double Get_MassPi0_hi() {return RangeMassPi0[1];};
    double Get_MassPipPim_lo() {return RangeMassPipPim[0];};
    double Get_MassPipPim_hi() {return RangeMassPipPim[1];};
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
    int Get_Topo_nElec() {return topo_nelec;};
    int Get_Topo_nPim() {return topo_npim;};
    int Get_Topo_nPip() {return topo_npip;};
    int Get_Topo_nGam() {return topo_ngam;};
    void Set_Topo_nElec(int nElec) {topo_nelec = nElec;};
    void Set_Topo_nPim(int nPim) {topo_npim = nPim;};
    void Set_Topo_nPip(int nPip) {topo_npip = nPip;};
    void Set_Topo_nGam(int nGam) {topo_ngam = nGam;};
    bool Check_ZDiff_ElecPim(double zdiff);
    bool Check_ZDiff_ElecPip(double zdiff);
    bool Check_MassPi0(double mass);
    bool Check_MassPipPim(double mass);
    bool Check_QSquared(double Qsq);
    bool Check_OpAng_ElecPhoton(double OpAng);
    bool Check_Wcut(double W);
    bool Check_NumDetPart(int nElec, int nPim, int nPip, int nGam);
    bool Check_ElectronR(double vr);
    bool Check_BetaPhoton(double beta);
    bool Check_MassOmega(double mass);
    bool Check_MassOmega_sb(double mass);
    void Print_Cuts();
    
    void InitCuts();
    void SetCut_ZDiff_ElecPion(double zdiff, int num);
    bool GetCut_ZDiff_ElecPion(int num);
    void SetCut_ZDiff_ElecPion_All();
    bool GetCut_ZDiff_ElecPion_All() {return cuts_omega_ZDiff;};
    void SetCut_MassPi0(double mass);
    bool GetCut_MassPi0() {return cuts_omega_MPi0;};
    void SetCut_MassPipPim(double mass);
    bool GetCut_MassPipPim() {return cuts_omega_MPipPim;};
    void SetCut_QSquared(double Qsq);
    bool GetCut_QSquared() {return cuts_omega_Q2;};
    void SetCut_OpAng_ElecPhoton(double OpAng, int num);
    bool GetCut_OpAng_ElecPhoton(int num);
    void SetCut_OpAng_ElecPhoton_All();
    bool GetCut_OpAng_ElecPhoton_All() {return cuts_omega_OpAng_ElecPhot;};
    void SetCut_Wcut(double W);
    bool GetCut_Wcut() {return cuts_omega_W;};
    void SetCut_NumDetPart(int nElec, int nPim, int nPip, int nGam);
    bool GetCut_NumDetPart() {return cuts_omega_NumDetPart;};
    void SetCut_MassOmega(double mass);
    bool GetCut_MassOmega() {return cuts_omega_MPipPimPi0;};
    void SetCut_MassOmega_sb(double mass);
    bool GetCut_MassOmega_sb() {return cuts_omega_MPipPimPi0_sb;};
    void SetCut_OmegaID();
    bool GetCut_OmegaID() {return cuts_omega_All;};
    void SetCut_OmegaID_woMassPi0();
    bool GetCut_OmegaID_woMassPi0() {return cuts_omega_woMPi0;};
    void SetCut_OmegaID_woQSquared();
    bool GetCut_OmegaID_woQSquared() {return cuts_omega_woQ2;};
    void SetCut_OmegaID_woWcut();
    bool GetCut_OmegaID_woWcut() {return cuts_omega_woW;};
    void SetCut_OmegaID_woZDiff();
    bool GetCut_OmegaID_woZDiff() {return cuts_omega_woZDiff;};
    void SetCut_OmegaID_woOpAng_ElecPhoton();
    bool GetCut_OmegaID_woOpAng_ElecPhoton() {return cuts_omega_woOpAng_ElecPhot;};
    void SetCut_OmegaID_woMassPipPim();
    bool GetCut_OmegaID_woMassPipPim() {return cuts_omega_woMPipPim;};
    void SetCut_OmegaID_woNumDetPart();
    bool GetCut_OmegaID_woNumDetPart() {return cuts_omega_woNumDetPart;};
};
#endif
