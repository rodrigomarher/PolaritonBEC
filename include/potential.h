#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>
#include "param.h"
#include "utils.h"
#include "grid.h"

class Potential{
    private:
        Params *_params;
        Grid *_grid;
        std::string _definition;
        std::vector<std::string> _pot_params;
        cdouble **_pot;
        void _parse();

        cdouble _cylinder(const double x, const double y);
    public:
        Potential (Params *params);
        void set_def(std::string def);
        void set_grid(Grid *grid);
        void set_potential();
        cdouble at(const int i, const int j);
        ~Potential();
};
#endif
