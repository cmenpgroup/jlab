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

void DetectedParticles::Print_DetPartLabel()
{
	int ii;
    cout<<"Detected Particles"<<endl;
    cout<<"=================="<<endl;
	for(ii=0;ii<this->Get_nDetPartLabel();ii++){
		cout << ii+1 << "\t" << this->Get_DetPartLabel(ii) << endl;
	}
}
