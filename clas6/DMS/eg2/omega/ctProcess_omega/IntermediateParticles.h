#ifndef INTERMEDIATEPARTICLES_H
#define INTERMEDIATEPARTICLES_H
#include <vector>
#include <string>

using namespace std;

class IntermediateParticles
{
public:
	vector<string> IntPartLabel;
	IntermediateParticles();
	int Get_nIntPartLabel();
	string Get_IntPartLabel(int num);
	void Print_IntPartLabel();
};

#endif
