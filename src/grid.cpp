#include "grid.h"
#include "utils.h"
#include <iostream>

Grid::Grid(Params *p){
    _params = p;
    _nx = _params->nx;
    _ny = _params->ny;
    _nt = _params->nt;
    _dx = _params->dx;
    _dy = _params->dy;
    _dt = _params->dt;
    
    _x = alloc1d<double>(_nx);
    _y = alloc1d<double>(_ny);
    _t = alloc1d<double>(_nt);
    linspace<double>(_x, _params->xmin, _params->xmax, _nx);
    linspace<double>(_y, _params->ymin, _params->ymax, _ny);
    std::cout<<_params->tmin<<" "<<_params->tmax<<std::endl;
    linspace<double>(_t, _params->tmin, _params->tmax, _nt);
    std::string path = _params->out_folder + "/t.dat";
    write_array1d<double>(_t,_nt, path);
    path = _params->out_folder + "/x.dat";
    write_array1d<double>(_x,_nx, path);
    path = _params->out_folder + "/y.dat";
    write_array1d<double>(_y,_ny, path);
}

double* Grid::get_x(){
    return _x;
}

double* Grid::get_y(){
    return _y;
}

double* Grid::get_t(){
    return _t;
}

double Grid::x(const int i){
    return _x[i];
}

double Grid::y(const int j){
    return _y[j];
}

double Grid::t(const int ti){
    return _t[ti];
}

Grid::~Grid(){
    free1d<double>(&_x, _nx);
    free1d<double>(&_y, _ny);
    free1d<double>(&_t, _nt);
}
