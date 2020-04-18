#include "MethodGas.h"
#include <cstdio>
#include <cstring>
#include <cmath>


void MethodGas::convertToParam(int i, Param& p)
{
    p.r = ro[i];
    p.u = ru[i]/p.r;
    p.v = rv[i]/p.r;
    p.e = re[i]/p.r-p.magU()/2.;
    p.p = p.r*p.e*(GAM-1.0);
    p.T = 0.0;
    p.M = p.magU()/sqrt(GAM*p.p/p.r);
}


void MethodGas::init()
{
    mesh = new Mesh();
    mesh->initFromFiles((char*)"step.1");

    ro = new double[mesh->cCount];
    ru = new double[mesh->cCount];
    rv = new double[mesh->cCount];
    re = new double[mesh->cCount];

    initValues();
    saveVTK(0);

    int_ro = new double[mesh->cCount];
    int_ru = new double[mesh->cCount];
    int_rv = new double[mesh->cCount];
    int_re = new double[mesh->cCount];

    TMAX    = 4.0;
    TAU     = 1.e-4;
}


void MethodGas::initValues()
{
    for (int i = 0; i < mesh->cCount; i++) {
        ro[i] = 1.;
        ru[i] = 3.;
        rv[i] = 0.0;
        re[i] = (1./1.4)/(GAM-1.0)+0.5*(ru[i]*ru[i])/ro[i];
    }
}


void MethodGas::run()
{
    double t = 0.0;
    int step = 0;
    while (t < TMAX) {
        t += TAU;
        step++;

        memset(int_ro, 0, sizeof(double)*mesh->cCount);
        memset(int_ru, 0, sizeof(double)*mesh->cCount);
        memset(int_rv, 0, sizeof(double)*mesh->cCount);
        memset(int_re, 0, sizeof(double)*mesh->cCount);
        for (int ie = 0; ie < mesh->eCount; ie++) {
            Edge &e = mesh->edges[ie];
            int c1 = e.c1;
            Cell &cell1 = mesh->getCell(c1);

            Param p1, p2;
            convertToParam(c1, p1);
            if (e.c2 >= 0) {
                convertToParam(e.c2, p2);
            }
            else {
                bnd(&e, p1, p2);
            }

            double fr, fu, fv, fe;
            flux(p1, p2, e.n, fr, fu, fv, fe);
            fr *= e.l;
            fu *= e.l;
            fv *= e.l;
            fe *= e.l;
            int_ro[c1] -= fr;
            int_ru[c1] -= fu;
            int_rv[c1] -= fv;
            int_re[c1] -= fe;
            if (e.c2 >= 0) {
                int_ro[e.c2] += fr;
                int_ru[e.c2] += fu;
                int_rv[e.c2] += fv;
                int_re[e.c2] += fe;
            }
        }

        for (int i = 0; i < mesh->cCount; i++) {
            double CFL = TAU/mesh->cells[i].S;
            ro[i] += int_ro[i]*CFL;
            ru[i] += int_ru[i]*CFL;
            rv[i] += int_rv[i]*CFL;
            re[i] += int_re[i]*CFL;
        }

        if (step % 100 == 0)
        {
            try {
                saveVTK(step);
            } catch (int code) {
                if (code == 1){
                    printf("Exception have been encountered: value is not a number");
                    break;
                } else {
                    printf("Exception have been caught, but it's impossible to identify it.");
                }
            }
            printf("Calculation results for step %d are saved.\n", step);
        }
    }
}


void MethodGas::flux(Param pl, Param pr, Vector n, double &fr, double &fu, double &fv, double &fe)
{
    double rol, rul, rvl, rel;
    double ror, rur, rvr, rer;
    double frl, ful, fvl, fel;
    double frr, fur, fvr, fer;
    double alpha, unl, unr, q1,q2;

    unl = n.x*pl.u+n.y*pl.v;
    unr = n.x*pr.u+n.y*pr.v;

    rol = pl.r;
    rul = pl.r*pl.u;
    rvl = pl.r*pl.v;
    rel = pl.r*(pl.e+0.5*pl.magU());

    ror = pr.r;
    rur = pr.r*pr.u;
    rvr = pr.r*pr.v;
    rer = pr.r*(pr.e+0.5*pr.magU());

    frl = pl.r*unl;
    ful = frl*pl.u+pl.p*n.x;
    fvl = frl*pl.v+pl.p*n.y;
    fel = (rel+pl.p)*unl;

    frr = pr.r*unr;
    fur = frr*pr.u+pr.p*n.x;
    fvr = frr*pr.v+pr.p*n.y;
    fer = (rer+pr.p)*unr;

    q1 = sqrt(pl.p*GAM/pl.r)+fabs(pl.magU());
    q2 = sqrt(pr.p*GAM/pr.r)+fabs(pr.magU());
    alpha = (q1 > q2) ? q1 : q2;

    fr = 0.5*((frl+frr)-alpha*(ror-rol));
    fu = 0.5*((ful+fur)-alpha*(rur-rul));
    fv = 0.5*((fvl+fvr)-alpha*(rvr-rvl));
    fe = 0.5*((fel+fer)-alpha*(rer-rel));
}


void MethodGas::bnd(Edge *e, Param p1, Param &p2)
{
    switch (e->type) {
        case 1: // вытекание
            p2 = p1;
            break;
        case 2: // втекание
            p2.r = 1.;
            p2.u = 3.;
            p2.v = 0.0;
            p2.p = 1./1.4;
            p2.e = p2.p/p2.r/(GAM-1.0);
            p2.T = 0.0;
            p2.M = 0.0;
            break;
        case 3: // отражение
            p2 = p1;
            double Un = p1.u*e->n.x + p1.v*e->n.y;
            Vector V;
            V.x = e->n.x*Un*2.0;
            V.y = e->n.y*Un*2.0;
            p2.u = p1.u - V.x;
            p2.v = p1.v - V.y;
            break;
    }
}


MethodGas::~MethodGas()
{
    delete mesh;
    delete[] ro, ru, rv, re;
    delete[] int_ro, int_ru, int_rv, int_re;
}
