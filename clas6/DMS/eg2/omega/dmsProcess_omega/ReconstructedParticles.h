#ifndef RECONSTRUCTEDPARTICLES_H
#define RECONSTRUCTEDPARTICLES_H
#include <vector>
#include <string>

using namespace std;

class ReconstructedParticles
{
public:
	vector<string> RecPartLabel;
	ReconstructedParticles();
	int Get_nRecPartLabel();
	string Get_RecPartLabel(int num);
	void Print_RecPartLabel();
};
#endif
