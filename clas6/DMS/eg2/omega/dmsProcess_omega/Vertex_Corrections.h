#ifndef ECGEOMETRY_H
#define ECGEOMETRY_H
#include <vector>
#include <string>

using namespace std;

class Vertex_Corrections
{
    int Sector;
    TVector3 TargetVert;
    TVector3 ParticleVert;
    TVector3 ParticleDir;
    TVector3 ParticleVertCorr;
    
public:
    Vertex_Corrections();
    void Put_Sector(int sector);
    void Put_Target_Vertex(TVector3 V3);
    void Put_Particle_Vertex(TVector3 V3);
    void Put_Particle_Dir(TVector3 V3);
    int Get_Sector() {return Sector;};
    TVector3 Get_Target_Vertex() {return TargetVert;};
    TVector3 Get_Particle_Vertex() {return ParticleVert;};
    TVector3 Get_Particle_Dir() {return ParticleDir;};
    TVector3 Get_Particle_Vertex_Corrected() {return ParticleVertCorr;};
    void Correct_Vertex();
    void Print_Vertex_Corrections();
};
#endif
