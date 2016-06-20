#ifndef PARTCOMB_PYTHIACPP_OMEGA_H
#define PARTCOMB_PYTHIACPP_OMEGA_H
#include <vector>
#include <string>

using namespace std;

class PartComb_pythiaCPP_Omega
{
    int nCtr;
    int nCombPhoton;
    vector<string> LabelPartComb;
    
public:
    PartComb_pythiaCPP_Omega();
    int Get_nCtr() {return nCtr;};
    int Get_nCombPhoton() {return nCombPhoton;};
    int Get_nPartComb() {return LabelPartComb.size();};
    string Get_PartComb(int num) {return LabelPartComb[num];};
};

#endif