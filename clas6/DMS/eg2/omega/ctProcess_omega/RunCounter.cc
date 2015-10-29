#include <cmath>
#include <vector>
#include <string>
#include "RunCounter.h"
#include <iostream>

RunCounter::RunCounter()
{
    CtrLabel.push_back("Total Events");
    CtrLabel.push_back("Electron ID (All)");
    CtrLabel.push_back("Electron ID (Mom)");
    CtrLabel.push_back("Electron ID (CCnphe)");
    CtrLabel.push_back("Electron ID (ECPvsP)");
    CtrLabel.push_back("Electron ID (ECin)");
    CtrLabel.push_back("Electron ID (dtECSC)");
    CtrLabel.push_back("Electron ID (ECfid)");
    CtrLabel.push_back("Photon ID (All)");
    CtrLabel.push_back("Photon ID (Mom)");
    CtrLabel.push_back("Photon ID (Beta)");
    CtrLabel.push_back("Photon ID (Timing)");
    CtrLabel.push_back("Photon ID (ECfid)");
    CtrLabel.push_back("Photon ID (TOF Msq)");
    CtrLabel.push_back("Charged Pion ID (All)");
    CtrLabel.push_back("Pos. Pion ID (dBeta)");
    CtrLabel.push_back("Neg. Pion ID (dBeta)");
    CtrLabel.push_back("Omega ID (All)");
    CtrLabel.push_back("Omega ID (Mpi0)");
    CtrLabel.push_back("Omega ID (Q2)");
    CtrLabel.push_back("Omega ID (W)");
    CtrLabel.push_back("Omega ID (MPipPim)");
    CtrLabel.push_back("Omega ID (Zmatch)");
    CtrLabel.push_back("Omega ID (OpAngElectron)");
    CtrLabel.push_back("Omega ID (NumDetPart)");
    
    this->Init();
}

void RunCounter::Init()
{
    CtrStats.assign(this->Get_nCtrLabel(),0);
}

void RunCounter::Increment(string label)
{
    bool ifound = false;
    
    for(int ii=0; ii<this->Get_nCtrLabel(); ii++){
        if (this->Get_CtrLabel(ii).compare(label)==0) {
            CtrStats[ii]++;
            ifound = true;
        }
    }
    if(!ifound) cout<<"RunCounter::Increment(), Unknown label "<<label<<endl;
}

void RunCounter::Print()
{
    int ii;
    cout<<"Statistics Summary"<<endl;
    cout<<"=================="<<endl;
    for(ii=0;ii<this->Get_nCtrLabel();ii++){
        cout << this->Get_CtrLabel(ii) << "\t" << this->Get_CtrStats(ii) << "\t";
        cout << "("<<(float)this->Get_CtrStats(ii)/(float)this->Get_CtrStats(0)<<")"<< endl;
    }
}


