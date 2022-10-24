#ifndef GRID_H
#define GRID_H
#include "param.h"


class Grid{
    private:
        double *_x;
        double *_y;
        double *_t;
        int _nx, _ny,_nt;
        double _dx, _dy,_dt;
        Params *_params;
    public:
        Grid(Params *p);
        double *get_x();
        double *get_y();
        double *get_t();
        double x(int i);
        double y(int j);
        double t(int ti);
        ~Grid();
};


#endif
