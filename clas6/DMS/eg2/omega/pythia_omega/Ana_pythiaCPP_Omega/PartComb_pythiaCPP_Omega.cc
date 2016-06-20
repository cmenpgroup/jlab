#include "PartComb_pythiaCPP_Omega.h"

PartComb_pythiaCPP_Omega::PartComb_pythiaCPP_Omega()
{
    nCtr = 6; // combination counter
    // 1 : Number of pi+
    // 2 : Number of pi-
    // 3 : Number of pi0
    // 4 : Number of 2 photons
    // 5 : Number of omega mesons
    // 6 : Number of particle combinations listed below
    
    nCombPhoton = 8; // number of particle combination with 2 photons
    
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPi0>=1");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPi0==1");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPhoton>=2");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPhoton==2");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPi0>=1 && nOmega>=1");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPi0==1 && nOmega==1");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega>=1");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPhoton==2 && nOmega==1");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPi0==1 && nOmega==0");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPhoton==2 && nOmega==0");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPi0>=1 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPi0==1 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("nPip>=1 && nPim>=1 && nPhoton>=2 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("nPip==1 && nPim==1 && nPhoton==2 && nOmega==0 && nEta==0");
}


