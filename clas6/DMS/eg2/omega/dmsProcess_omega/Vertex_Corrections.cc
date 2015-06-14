#include <vector>
#include <string>
#include "Vertex_Corrections.h"
#include <iostream>
#include <cmath>
Vertex_Corrections::Vertex_Corrections()
{
    // Initialize sector
    Sector = 0;
    
    // Initialize the target vertex
    TargetVert.SetX(0.0);
    TargetVert.SetY(0.0);
    TargetVert.SetZ(0.0);
    
    // Initialize the particle vertex
    ParticleVert.SetX(0.0);
    ParticleVert.SetY(0.0);
    ParticleVert.SetZ(0.0);
    
    // Initialize the particle directional vector
    ParticleDir.SetX(0.0);
    ParticleDir.SetY(0.0);
    ParticleDir.SetZ(0.0);

    // Initialize the particle's corrected vertex
    ParticleVertCorr.SetX(0.0);
    ParticleVertCorr.SetY(0.0);
    ParticleVertCorr.SetZ(0.0);
}

void Vertex_Corrections::Put_Target_Vertex(TVector3 V3){
    this->TargetVert = V3;
}

void Vertex_Corrections::Put_Particle_Vertex(TVector3 V3){
    this->ParticleVert = V3;
}

void Vertex_Corrections::Put_Particle_Dir(TVector3 V3){
    this->ParticleDir = V3;
}

void Vertex_Corrections::Put_Sector(int sector){
    this->Sector = sector;
}

void Vertex_Corrections::Correct_Vertex(){
    //Vertex correction (x,y,z)->(x_corr,y_corr,z_corr). Correction preserves (x,y) in coordinates of Sector1
    TVector3 RotatedVertPos = this->Get_Particle_Vertex();
    TVector3 RotatedVertDir = this->Get_Particle_Dir();
    TVector3 TargetPos = this->Get_Target_Vertex();
    int sect = this->Get_Sector();
    
    RotatedVertPos.RotateZ(-TMath::DegToRad()*60.*sect);
    RotatedVertDir.RotateZ(-TMath::DegToRad()*60.*sect);
    TargetPos.RotateZ(-TMath::DegToRad()*60.*sect);
    
    Float_t ShiftLength = (TargetPos.X()-RotatedVertPos.X())/RotatedVertDir.X();
    RotatedVertDir = ShiftLength*RotatedVertDir;
    RotatedVertPos = RotatedVertPos+RotatedVertDir;

    // set the particle's correct vertex
    ParticleVertCorr.SetXYZ((RotatedVertPos-TargetPos).X(),(RotatedVertPos-TargetPos).Y(),RotatedVertPos.Z());
}

// print the cut information
void Vertex_Corrections::Print_Vertex_Corrections()
{
    cout<<"Vertex Corrections"<<endl;
    cout<<"========================="<<endl;
    
    cout<<"Particle Vertex: "<<endl;
    TVector3 PartPos = this->Get_Particle_Vertex();
    PartPos.Print();
    
    cout<<"Particle Direction: "<<endl;
    TVector3 PartDir = this->Get_Particle_Dir();
    PartDir.Print();
    
    cout<<"Target: "<<endl;
    TVector3 Target = this->Get_Target_Vertex();
    Target.Print();
    
    cout<<"Corrected Vertex: "<<endl;
    TVector3 PartPosCorr = this->Get_Particle_Vertex_Corrected();
    PartPosCorr.Print();
    
}
