#ifndef VERTEXCORRECTIONS_H
#define VERTEXCORRECTIONS_H
#include <vector>
#include <string>
#include "TVector3.h"

using namespace std;

class Vertex_Corrections
{
    double ParticlePhi;
    TVector3 TargetVert;
    TVector3 ParticleVert;
    TVector3 ParticleDir;
    TVector3 ParticleVertCorr;
    
public:
    Vertex_Corrections();
    void Put_Particle_Phi(double phi);
    void Put_Target_Vertex(TVector3 V3);
    void Put_Particle_Vertex(TVector3 V3);
    void Put_Particle_Dir(TVector3 V3);
    double Get_Particle_Phi() {return ParticlePhi;};
    TVector3 Get_Target_Vertex() {return TargetVert;};
    TVector3 Get_Particle_Vertex() {return ParticleVert;};
    TVector3 Get_Particle_Dir() {return ParticleDir;};
    TVector3 Get_Particle_Vertex_Corrected() {return ParticleVertCorr;};
    void Correct_Vertex();
    void Print_Vertex_Corrections();
};
#endif
