#ifndef filter_hh
#define filter_hh 1

#include "cpp_headers.hh"

class Filter
{
    int KScut;
    vector<int> partType;
    vector<int> partQty;
    
    public:
    Filter();
    ~Filter();

    void init();
    int Get_nPartType() {return partType.size();};
    int Get_nPartQty() {return partQty.size();};
    int GetPartType(int num) {return partType[num];};
    int GetPartQty(int num) {return partQty[num];};
    void SetKScut(int ks);
    int GetKScut() {return KScut;};
    bool Cut();
};

#endif
