#ifndef DETECTEDPARTICLES_H
#define DETECTEDPARTICLES_H
#include <vector>
#include <string>

using namespace std;

class DetectedParticles
{
public:
	vector<string> DetPartLabel;
	DetectedParticles();
	int Get_nDetPartLabel();
	string Get_DetPartLabel(int num);
    int Get_PartIndex(string label);
	void Print_DetPartLabel();
};
#endif
