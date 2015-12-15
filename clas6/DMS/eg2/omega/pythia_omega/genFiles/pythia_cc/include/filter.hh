#ifndef filter_hh
#define filter_hh 1

#include "cpp_headers.hh"

class Filter
{
    int KScut;
    vector<int> partType; // list of final particles
    vector<int> partQty;  // list of the number of each final particle
    vector<int> partCtr;  // counter for each final particle per event
    
    public:
    Filter();
    ~Filter();

    void init();
    int Get_nPartType() {return partType.size();};
    int Get_nPartQty() {return partQty.size();};
    int GetPartType(int num) {return partType[num];};
    int GetPartQty(int num) {return partQty[num];};
    bool CheckPartSize();
    void SetKScut(int ks);
    int GetKScut() {return KScut;};
    bool Cut();
};

#endif
