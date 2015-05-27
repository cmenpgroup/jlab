#include <vector>
#include <string>
#include "EC_geometry.h"
#include <iostream>
#include <cmath>
EC_geometry::EC_geometry()
{
    EC_theta = 0.4363323; // radians
    ylow = -182.974;
    yhi = 189.956;
    rho = 1.097621; // radians
    
    EC_U = -1.0;
    EC_V = -1.0;
    EC_W = -1.0;
}

double EC_geometry::Get_Xlocal()
{
    double ret = 0.0;
    if(this->Check_U() && this->Check_V() && this->Check_W()) ret = this->Get_W()*cos(this->Get_rho()) - 0.5*this->Get_V();
    return ret;
}

double EC_geometry::Get_Ylocal()
{
    double ret = 0.0;
    if(this->Check_U() && this->Check_V() && this->Check_W()) ret = this->Get_ylow() + this->Get_U()*sin(this->Get_rho());
    return ret;
}

void EC_geometry::Put_UVW(double u, double v, double w){
    this->EC_U = u;
    this->EC_V = v;
    this->EC_W = w;
}

// check that the U-view is a positive number (ie has been filled)
bool EC_geometry::Check_U()
{
    bool ret = (this->Get_U() >= 0) ? true : false;
    return ret;
}

// check that the V-view is a positive number (ie has been filled)
bool EC_geometry::Check_V()
{
    bool ret = (this->Get_V() >= 0) ? true : false;
    return ret;
}

// check that the W-view is a positive number (ie has been filled)
bool EC_geometry::Check_W()
{
    bool ret = (this->Get_W() >= 0) ? true : false;
    return ret;
}

// print the cut information
void EC_geometry::Print_EC_geometry()
{
    cout<<"EC Geometry Constants"<<endl;
    cout<<"========================="<<endl;

    cout << "EC theta = " << this->Get_EC_theta()*180.0/3.14159 << endl;
    cout << "ylow = " << this->Get_ylow() << " (cm)" << endl;
    cout << "yhi = " << this->Get_yhi() << " (cm)" << endl;
    cout << "rho = " << this->Get_rho()*180.0/3.14159 << " (cm)" << endl;
    cout << endl;
}
