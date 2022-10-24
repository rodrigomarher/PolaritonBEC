#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>
#include "utils.h"
#include "param.h"
#include "grid.h"
#include "subfield.h"

class Field{
    private:
        Params *_params;
        Grid *_grid;
        std::vector<subField*> _field_container;
        std::vector<std::string> _def_vec;
        cdouble **_field;

        void _parse();
        
    public:
        Field(Params *params);
        void set_grid(Grid *grid);
        void set_field();
        cdouble** get_field();
        void save_field(std::string filename);        
        void update(const int ti);
        cdouble at(const int ti, const int j);
        ~Field();
};

#endif
