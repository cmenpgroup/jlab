#include <vector>
#include <string>
#include "ParticleList.h"
#include "DetectedParticles.h"
#include "IntermediateParticles.h"
#include "ReconstructedParticles.h"
#include <iostream>
ParticleList::ParticleList()
{
    vector<string> temp;
    DetectedParticles DetList;
    IntermediateParticles IntList;
    ReconstructedParticles RecList;
 
    int ii;
	for(ii=0;ii<DetList.Get_nDetPartLabel();ii++){
        PartLabel.push_back(DetList.Get_DetPartLabel(ii));
    }
	for(ii=0;ii<IntList.Get_nIntPartLabel();ii++){
        PartLabel.push_back(IntList.Get_IntPartLabel(ii));
    }
	for(ii=0;ii<RecList.Get_nRecPartLabel();ii++){
        PartLabel.push_back(RecList.Get_RecPartLabel(ii));
    }
}

void ParticleList::Print_PartLabel()
{
	int ii;
    cout<<"All Particles in Analysis"<<endl;
    cout<<"========================="<<endl;
	for(ii=0;ii<this->Get_nPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_PartLabel(ii) << endl;
	}
}

int ParticleList::Get_nPartLabel()
{
    return PartLabel.size();
}

string ParticleList::Get_PartLabel(int num)
{
    return PartLabel[num];
}
