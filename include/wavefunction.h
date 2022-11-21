#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <complex>
#include "param.h"
#include "grid.h"

class Wavefunction{
    private:
        cdouble **_wf;
        cdouble ***_wfbuf;
        cdouble *_x_row;
        cdouble *_y_row;
        Params *_params;
        Grid *_grid;
    public:
        Wavefunction(Params *params);
        void set_grid(Grid *grid);
        void init_zeros();
        cdouble** get_wf();
        cdouble* get_x_row(const int j);
        cdouble* get_y_row(const int i);
        void set_x_row(cdouble *arr, const int j);
        void set_y_row(cdouble *arr, const int i);
        cdouble norm();
        cdouble at(const int i, const int j);
        void copy_to_buf(const int ti);
        void operator/(const cdouble k);
        void save_wf(std::string filepath,const int ti);
        ~Wavefunction();
};

#endif
