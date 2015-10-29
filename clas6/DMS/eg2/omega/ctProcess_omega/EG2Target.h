#ifndef EG2TARGET_H
#define EG2TARGET_H
#include <vector>
#include <string>

using namespace std;

class EG2Target
{
    vector<string> Label;
    vector<int> Index;
    vector<double> RangeLD2;
    vector<double> RangeNuc;
public:
    EG2Target();
    int Get_nLabel() {return Label.size();};
    int Get_nIndex() {return Index.size();};
	string Get_Label(int num) {return Label[num];};
    double Get_LD2_lo() {return RangeLD2[0];};
    double Get_LD2_hi() {return RangeLD2[1];};
    double Get_Nuc_lo() {return RangeNuc[0];};
    double Get_Nuc_hi() {return RangeNuc[1];};
    int Get_Index(double z);
	void Print_Info();
};
#endif
