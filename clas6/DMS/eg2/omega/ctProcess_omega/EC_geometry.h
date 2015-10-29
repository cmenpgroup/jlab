#ifndef ECGEOMETRY_H
#define ECGEOMETRY_H
#include <vector>
#include <string>

using namespace std;

class EC_geometry
{
    double EC_U;
    double EC_V;
    double EC_W;
    double EC_theta;
    double EC_phi;
    double ylow;
    double yhi;
    double rho;
public:
    EC_geometry();
    void Put_UVW(double u, double v, double w);
    double Get_EC_theta() {return EC_theta;};
    double Get_ylow() {return ylow;};
    double Get_yhi() {return yhi;};
    double Get_rho() {return rho;};
    double Get_U() {return EC_U;};
    double Get_V() {return EC_V;};
    double Get_W() {return EC_W;};
    double Get_Xlocal();
    double Get_Ylocal();
    bool Check_U();
    bool Check_V();
    bool Check_W();
    void Print_EC_geometry();
};
#endif
