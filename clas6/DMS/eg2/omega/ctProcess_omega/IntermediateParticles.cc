#include <vector>
#include <string>
#include "IntermediateParticles.h"
#include <iostream>
IntermediateParticles::IntermediateParticles()
{
	IntPartLabel.push_back("Pi0");
}

void IntermediateParticles::Print_IntPartLabel()
{
	int ii;
    cout<<"Intermediate Particles"<<endl;
    cout<<"======================"<<endl;
	for(ii=0;ii<this->Get_nIntPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_IntPartLabel(ii) << endl;
	}
}

int IntermediateParticles::Get_nIntPartLabel()
{
    return IntPartLabel.size();
}

string IntermediateParticles::Get_IntPartLabel(int num)
{
    return IntPartLabel[num];
}
