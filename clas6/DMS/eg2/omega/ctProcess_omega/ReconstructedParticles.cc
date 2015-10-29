#include <vector>
#include <string>
#include "ReconstructedParticles.h"
#include <iostream>
ReconstructedParticles::ReconstructedParticles()
{
	RecPartLabel.push_back("Omega");
}

void ReconstructedParticles::Print_RecPartLabel()
{
	int ii;
    cout<<"Reconstructed Particles"<<endl;
    cout<<"======================="<<endl;
	for(ii=0;ii<this->Get_nRecPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_RecPartLabel(ii) << endl;
	}
}

int ReconstructedParticles::Get_nRecPartLabel()
{
    return RecPartLabel.size();
}

string ReconstructedParticles::Get_RecPartLabel(int num)
{
    return RecPartLabel[num];
}
