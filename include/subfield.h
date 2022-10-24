#ifndef SUBFIELD_H
#define SUBFIELD_H

#include <vector>
#include <string>
#include "param.h"
#include "grid.h"
#include "utils.h"

class subField {
    private:
        double _fp, _tmin, _tmax;
        double _w0, _x0, _y0;
        int _l, _p;
        double _deltaw, _kpx, _kpy;

        std::string _definition;
        Grid *_grid;

        cdouble _gaussian(const double x, const double y, const cdouble tenv);
        cdouble _laguerre(const double x, const double y, const cdouble tenv);
        cdouble _tenv_sin2(const int ti);
        void _parse();

    public:
        subField (std::string def);
        void set_grid(Grid *grid);
        cdouble at(const int i, const int j, const int ti);
};

#endif
