#include <vector>
#include <string>
#include "EG2Target.h"
#include <iostream>
EG2Target::EG2Target()
{
    Label.push_back("NoTarget");
    Label.push_back("LD2");
    Label.push_back("Nuc");
    
    Index.push_back(0);
    Index.push_back(1);
    Index.push_back(2);
    
    RangeLD2.push_back(-32.0);
    RangeLD2.push_back(-28.0);

    RangeNuc.push_back(-26.0);
    RangeNuc.push_back(-23.0);
}

// Return the CLAS eg2 target index from vertex Z.
//
// Return 0 = outside target limits, 1 = liquid deuterium, 2 = nuclear
// z must be given in cm
//
int EG2Target::Get_Index(double z)
{
    
    int ret = 0; // init the return variable
    
    if (z >= this->Get_LD2_lo() && z < this->Get_LD2_hi()) {
        ret = Index[1]; // deuterium target
    } else if (z >= this->Get_Nuc_lo() && z < this->Get_Nuc_hi()) {
        ret = Index[2]; // nuclear target
    } else {
        ret = Index[0]; // no target
    }
    
    return ret;
}

void EG2Target::Print_Info()
{
	int ii;
    cout<<"EG2 Target Info"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nLabel();ii++){
        cout << ii+1 << "\t" << this->Get_Label(ii) << endl;
    }
    
    cout << "LD2 target: " << this->Get_LD2_lo() << " , " << this->Get_LD2_hi() << " (cm)" << endl;
    cout << "Nuclear target: " << this->Get_Nuc_lo() << " , " << this->Get_Nuc_hi() << " (cm)" << endl;
}
