#ifndef POLARITON_H
#define POLARITON_H

#include "grid.h"
#include "param.h"
#include "wavefunction.h"
#include "field.h"
#include "rk.h"

class Polariton{
    private:
        Grid *_grid;
        Params *_params;
        Wavefunction *_wfc, *_wfx;
        RK *_rk;
        Field* _field;
    public:
        Polariton(Params *params);
        void set_param(Params *params);
        void evolve();
        void snapshot(const int ti);
};

#endif
