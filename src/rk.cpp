#include <iostream>
#include "rk.h"
#include "utils.h"

RK::RK(Params* params){
    _params = params;
    
    _k1 = alloc2d<cdouble>(_params->nx, _params->ny);
    _k2 = alloc2d<cdouble>(_params->nx, _params->ny);
    _k3 = alloc2d<cdouble>(_params->nx, _params->ny);
    _k4 = alloc2d<cdouble>(_params->nx, _params->ny);
}

cdouble RK::_rhs_f(const int i, const int j, const int ti, cdouble** kn, double pre){
    cdouble laplacian = _laplacian(i,j, _wf, kn, pre); 
    cdouble fp = _field->at(i,j);
    cdouble a = fp;
    cdouble b = (_pot->at(i,j) - 0.5*I*hbar*_params->gamma_c)*(_wf->at(i,j)+pre*kn[i][j]);
    cdouble c = _params->g*(_wf->at(i,j))*std::conj(_wf->at(i,j))*(_wf->at(i,j)+pre*kn[i][j]);
    cdouble d = _params->omega_c*(_wf->at(i,j)+pre*kn[i][j]) - 0.5*hbar*hbar/_params->mc*laplacian;
    return -I/hbar*(a+b+c+d);
}

cdouble RK::_laplacian(const int i, const int j,Wavefunction* wf, cdouble** kn, double pre){
    cdouble laplacian_x;
    cdouble laplacian_y;
    if(i == 0){
        laplacian_x = -2.0*(wf->at(i,j)+pre*kn[i][j]);
        laplacian_x += (wf->at(i+1,j)+pre*kn[i+1][j]);
    }
    else if(i==_params->nx-1){
        laplacian_x = (wf->at(i-1,j)+pre*kn[i-1][j]);
        laplacian_x +=  - 2.0*(wf->at(i,j)+pre*kn[i][j]);
    }
    else{
        laplacian_x = (wf->at(i-1,j)+pre*kn[i-1][j]);
        laplacian_x +=  - 2.0*(wf->at(i,j)+pre*kn[i][j]);
        laplacian_x += (wf->at(i+1,j)+pre*kn[i+1][j]);
    }
    laplacian_x /= _params->dx*_params->dx;

    if(j == 0){
        laplacian_y = -2.0*(wf->at(i,j)+pre*kn[i][j]);
        laplacian_y += (wf->at(i,j+1)+pre*kn[i][j+1]);
    }
    else if(j ==_params->ny-1){
        laplacian_y = (wf->at(i,j-1)+pre*kn[i][j-1]);
        laplacian_y +=  - 2.0*(wf->at(i,j)+pre*kn[i][j]);
    }
    else{
        laplacian_y = (wf->at(i,j-1)+pre*kn[i][j-1]);
        laplacian_y +=  - 2.0*(wf->at(i,j)+pre*kn[i][j]);
        laplacian_y += (wf->at(i,j+1)+pre*kn[i][j+1]);
    }
    laplacian_y /= _params->dy*_params->dy;
    return laplacian_x + laplacian_y;
}

void RK::_update_k(const int ti){
    _update_k1(ti);
    _update_k2(ti);
    _update_k3(ti);
    _update_k4(ti);
}

void RK::_update_k1(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            _k1[i][j] = _params->dt*_rhs_f(i,j,ti,_k1,0.0);
        }
    }
}

void RK::_update_k2(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            _k2[i][j] = _params->dt*_rhs_f(i,j,ti,_k1,0.5);
        }
    }
}

void RK::_update_k3(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            _k3[i][j] = _params->dt*_rhs_f(i,j,ti,_k2,0.5);
        }
    }
}

void RK::_update_k4(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx;i++){
        for(int j=0; j<_params->ny; j++){
            _k4[i][j] = _params->dt*_rhs_f(i,j,ti,_k3,1.0);
        }
    }
}

void RK::set_grid(Grid* grid){
    _grid = grid;
}

void RK::set_wf(Wavefunction* wf1){
    _wf = wf1; 
}

void RK::set_field(Field* field){
    _field = field;
}

void RK::set_mask(Mask* mask){
    _mask = mask;
}

void RK::set_pot(Potential *pot){
    _pot = pot;
}

void RK::step(const int ti){
    _update_k(ti);

    cdouble** wf = _wf->get_wf();
    cdouble u = 1.0/6.0;

    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            wf[i][j] += u*(_k1[i][j] + 2.0*_k2[i][j] + 2.0*_k3[i][j] + _k4[i][j])*_mask->x(i)*_mask->y(j);
        }
    }
}

RK::~RK(){
    free2d<cdouble>(&_k1,_params->nx,_params->ny);
    free2d<cdouble>(&_k2,_params->nx,_params->ny);
    free2d<cdouble>(&_k3,_params->nx,_params->ny);
    free2d<cdouble>(&_k4,_params->nx,_params->ny);
}
