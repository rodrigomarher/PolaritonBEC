#include "wavefunction.h"
#include "utils.h"

Wavefunction::Wavefunction(Params *params){
    _params = params;
    _wf = alloc2d<cdouble>(_params->nx,_params->ny);
    _wfbuf = alloc3d<cdouble>(_params->nbuf, _params->nx, _params->ny);
    _x_row = alloc1d<cdouble>(_params->nx);
    _y_row = alloc1d<cdouble>(_params->ny);
    init_zeros();

}

void Wavefunction::set_grid(Grid *grid){
    _grid=grid;
}

void Wavefunction::init_zeros(){
    for(int i=0; i<_params->nx;i++){
        for(int j=0; j<_params->ny;j++){
            _wf[i][j] = std::complex<double>(0.0,0.0);
        }
    }
}

cdouble** Wavefunction::get_wf(){
    return _wf;
}

cdouble* Wavefunction::get_x_row(const int j){
    for(int i=0; i<_params->nx;i++){
        _x_row[i] = _wf[i][j];
    }
    return _x_row;
}

cdouble* Wavefunction::get_y_row(const int i){
    for(int j=0; j<_params->ny;j++){
        _y_row[j] = _wf[i][j];
    }
    return _y_row;
}

cdouble Wavefunction::at(const int i, const int j){
    return _wf[i][j];
}

cdouble Wavefunction::norm(){
    cdouble sum = 0.0;
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            sum += conj(_wf[i][j])*_wf[i][j];
        }
    }
    return sum*_params->dx*_params->dy;
}

void Wavefunction::set_x_row(cdouble *arr, const int j){
    for(int i=0; i<_params->nx; i++){
        _wf[i][j] = arr[i];
    } 
}

void Wavefunction::set_y_row(cdouble *arr, const int i){
    for(int j=0; j<_params->ny; j++){
        _wf[i][j] = arr[j];
    }
}

void Wavefunction::copy_to_buf(const int ti){
    for (int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            _wfbuf[ti][i][j] = _wf[i][j];
        }
    }

}

void Wavefunction::operator/(const cdouble k){
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            _wf[i][j] /= k;
        }
    }
}

void Wavefunction::save_wf(std::string filepath,const int ti){
    write_array2d_complex(_wfbuf[ti],_params->nx,_params->ny,filepath); 
}

Wavefunction::~Wavefunction(){
    free2d<cdouble>(&_wf, _params->nx, _params->ny);
    free3d<cdouble>(&_wfbuf, _params->nbuf, _params->nx, _params->ny);
    free1d<cdouble>(&_x_row, _params->nx);
    free1d<cdouble>(&_y_row, _params->ny);
}
