#include "MethodHeat.h"
#include <cstdio>
#include <cstring>


void MethodHeat::convertToParam(int i, Param& p)
{
    p.r = 0.0;
    p.p = 0.0;
    p.u = 0.0;
    p.v = 0.0;
    p.e = 0.0;
    p.T = T[i];
    p.M = 0.0;
}


void MethodHeat::init()
{
    mesh = new Mesh();
    mesh->initFromFiles((char*)"heat.1");

    T = new double[mesh->cCount];
    memset(T, 0, sizeof(double)*mesh->cCount);

    saveVTK(0);

}


void MethodHeat::run()
{

}

MethodHeat::~MethodHeat()
{
    delete[] T;
}
