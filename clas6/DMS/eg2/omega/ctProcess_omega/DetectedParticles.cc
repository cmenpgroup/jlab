#include <vector>
#include <string>
#include "DetectedParticles.h"
#include <iostream>
DetectedParticles::DetectedParticles()
{
	DetPartLabel.push_back("Electron");
	DetPartLabel.push_back("Pi-");
	DetPartLabel.push_back("Pi+");
	DetPartLabel.push_back("Photon1");
	DetPartLabel.push_back("Photon2");
}

int DetectedParticles::Get_nDetPartLabel() {
    return DetPartLabel.size();
}

string DetectedParticles::Get_DetPartLabel(int num) {
    return DetPartLabel[num];
}

int DetectedParticles::Get_PartIndex(string label)
{
    int ret = -1;
    bool ifound = false;
    
    for(int ii=0; ii<this->Get_nDetPartLabel(); ii++){
        if (this->Get_DetPartLabel(ii).compare(label)==0) {
            ret = ii;
            ifound = true;
        }
    }
    if(!ifound) cout<<"DetectedParticle::Get_PartIndex(), Unknown label "<<label<<endl;
    return ret;
}

void DetectedParticles::Print_DetPartLabel()
{
	int ii;
    cout<<"Detected Particles"<<endl;
    cout<<"=================="<<endl;
	for(ii=0;ii<this->Get_nDetPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_DetPartLabel(ii) << endl;
	}
}
