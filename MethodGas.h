#pragma once

#include "method.h"


class MethodGas : public Method
{
public:
    virtual void convertToParam(int, Param&);
    virtual void init();
    virtual void run();

    ~MethodGas();
protected:

    void initValues();
    void flux(Param pl, Param pr, Vector n, double &fr, double &fu, double &fv, double &fe);
    void bnd(Edge *e, Param p1, Param &p2);

    double *ro, *ru, *rv, *re;
    double *int_ro, *int_ru, *int_rv, *int_re;

    double TMAX;
    double TAU;

    const double GAM = 1.4;
};



