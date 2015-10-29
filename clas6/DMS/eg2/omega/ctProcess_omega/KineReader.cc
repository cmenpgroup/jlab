#include <cmath>
#include <vector>
#include <string>
#include "KineReader.h"
#include <iostream>
KineReader::KineReader(TTree *tree)
{
    b_kine = tree->GetBranch("Kinematics");
    l_evtnum = b_kine->GetLeaf("EvtNum");
    l_q2 = b_kine->GetLeaf("Q2");
    l_nu = b_kine->GetLeaf("Nu");
    l_xb = b_kine->GetLeaf("Xb");
    l_w = b_kine->GetLeaf("W");
    l_xcorr = b_kine->GetLeaf("Xcorr");
    l_ycorr = b_kine->GetLeaf("Ycorr");
    l_zcorr = b_kine->GetLeaf("Zcorr");
    l_nelec = b_kine->GetLeaf("nElec");
    l_npip = b_kine->GetLeaf("nPip");
    l_npim = b_kine->GetLeaf("nPim");
    l_ngam = b_kine->GetLeaf("nGam");
}

void KineReader::ReadEntry(int num)
{
    b_kine->GetEntry(num);
}


