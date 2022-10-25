#include <sstream>
#include <iostream>
#include "potential.h"

Potential::Potential(Params *params){
    _params = params;
    _pot = alloc2d<cdouble>(_params->nx,_params->ny);
}

void Potential::set_def(std::string def){
    _definition = def;
    _parse();
}

void Potential::set_grid(Grid *grid){
    _grid = grid;
}

void Potential::_parse(){
    char delim = ',';
    int counter = 0;
    std::string item;
    std::stringstream ss(_definition);
    while(getline(ss,item,delim)){
        _pot_params.push_back(item);
        counter++;
    }
}

cdouble Potential::_cylinder(const double x, const double y){
    double V0 = std::stod(_pot_params[1]);
    double r0 = std::stod(_pot_params[2]);
    double x0 = std::stod(_pot_params[3]);
    double y0 = std::stod(_pot_params[4]);
    double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0);

    if(r2<r0*r0){
        return cdouble(V0,0.0);
    }
    return cdouble(0.0,0.0);
}


void Potential::set_potential(){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx;i++){
        for(int j=0; j<_params->ny; j++){
            if(_pot_params[0]=="cylinder"){
                _pot[i][j] = _cylinder(_grid->x(i), _grid->y(j));
            }
            else{
                _pot[i][j] = cdouble(0.0,0.0);
            }
        }
    }
}

cdouble Potential::at(const int i, const int j){
    return _pot[i][j];
}

Potential::~Potential(){
    free2d<cdouble>(&_pot,_params->nx,_params->ny);
}
