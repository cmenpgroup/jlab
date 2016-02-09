#ifndef filter_hh
#define filter_hh 1

#include "cpp_headers.hh"

class Filter
{
    int KScut;
    double Q2cut_lo;
    double Q2cut_hi;
    vector<int> partType; // list of final particles
    vector<int> partQty;  // list of the number of each final particle
    vector<int> partCtr;  // counter for each final particle per event
    
    public:
    Filter();
    ~Filter();

    void init();
    int Get_nPartType() {return partType.size();};
    int Get_nPartQty() {return partQty.size();};
    int Get_nPartCtr() {return partCtr.size();};
    int GetPartType(int num) {return partType[num];};
    int GetPartQty(int num) {return partQty[num];};
    int GetPartCtr(int num) {return partCtr[num];};
    double GetQ2cut_lo() {return Q2cut_lo;};
    double GetQ2cut_hi() {return Q2cut_hi;};
    void ZeroPartCtr();
    bool CheckPartSize();
    void SetKScut(int ks);
    void SetQ2cut_lo(double Q2);
    void SetQ2cut_hi(double Q2);
    int GetKScut() {return KScut;};
    bool Cut();
};

#endif
