#ifndef MASK_H
#define MASK_H

#include "param.h"
#include "grid.h"
#include "utils.h"

class Mask{
    private:
        Params *_param;
        Grid *_grid;
        cdouble *_xmask, *_ymask;

    public:
        Mask(Params *param);
        void set_grid(Grid *grid);
        void set_mask();
        cdouble x(const int i);
        cdouble y(const int j);
        void save_mask();
        ~Mask();
};


#endif
