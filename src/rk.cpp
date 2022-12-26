#include "rk.h"
#include "utils.h"

RK::RK(Params* params){
    _params = params;
    
    _k1x = alloc2d<cdouble>(_params->nx, _params->ny);
    _k2x = alloc2d<cdouble>(_params->nx, _params->ny);
    _k3x = alloc2d<cdouble>(_params->nx, _params->ny);
    _k4x = alloc2d<cdouble>(_params->nx, _params->ny);
    _k1c = alloc2d<cdouble>(_params->nx, _params->ny);
    _k2c = alloc2d<cdouble>(_params->nx, _params->ny);
    _k3c = alloc2d<cdouble>(_params->nx, _params->ny);
    _k4c = alloc2d<cdouble>(_params->nx, _params->ny);
}

cdouble RK::_rhs_fc(const int i, const int j, const int ti, cdouble** knc, cdouble** knx, double pre){
    cdouble laplacian = _laplacian(i,j, _wfc, knc, pre); 
    cdouble fp = _field->at(i,j);
    cdouble a = fp + _params->omega_r*(_wfx->at(i,j)+pre*knx[i][j]);
    cdouble b = (_ph_pot->at(i,j) - 0.5*I*hbar*_params->gamma_c)*(_wfc->at(i,j)+pre*knc[i][j]);
    cdouble c = _params->omega_c*(_wfc->at(i,j)+pre*knc[i][j]) - 0.5*hbar*hbar/_params->mc*laplacian;
    return -I/hbar*(a+b+c);
}

cdouble RK::_rhs_fx(const int i, const int j, const int ti, cdouble** knc, cdouble** knx, double pre){

    cdouble laplacian = _laplacian(i,j, _wfx, knx, pre); 
    cdouble a = _params->omega_r*(_wfc->at(i,j)+pre*knc[i][j]);
    cdouble b = (_ex_pot->at(i,j)- 0.5*I*hbar*_params->gamma_x)*(_wfx->at(i,j)+pre*knx[i][j]);
    cdouble c = _params->g*std::abs(_wfx->at(i,j))*std::abs(_wfx->at(i,j))*(_wfx->at(i,j)+pre*knx[i][j]);
    cdouble d = _params->omega_x*(_wfx->at(i,j)+pre*knx[i][j]) - 0.5*hbar*hbar/_params->mx*laplacian;
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
            _k1c[i][j] = _params->dt*_rhs_fc(i,j,ti,_k1c,_k1x,0.0);
            _k1x[i][j] = _params->dt*_rhs_fx(i,j,ti,_k1c,_k1x,0.0);
        }
    }
}

void RK::_update_k2(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            _k2c[i][j] = _params->dt*_rhs_fc(i,j,ti,_k1c,_k1x,0.5);
            _k2x[i][j] = _params->dt*_rhs_fx(i,j,ti,_k1c,_k1x,0.5);
        }
    }
}

void RK::_update_k3(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            _k3c[i][j] = _params->dt*_rhs_fc(i,j,ti,_k2c,_k2x,0.5);
            _k3x[i][j] = _params->dt*_rhs_fx(i,j,ti,_k2c,_k2x,0.5);
        }
    }
}

void RK::_update_k4(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx;i++){
        for(int j=0; j<_params->ny; j++){
            _k4c[i][j] = _params->dt*_rhs_fc(i,j,ti,_k3c,_k3x,1.0);
            _k4x[i][j] = _params->dt*_rhs_fx(i,j,ti,_k3c,_k3x,1.0);
        }
    }
}

void RK::set_grid(Grid* grid){
    _grid = grid;
}

void RK::set_wf(Wavefunction* wf1, Wavefunction* wf2){
    _wfc = wf1; _wfx = wf2;
}

void RK::set_field(Field* field){
    _field = field;
}

void RK::set_mask(Mask* mask){
    _mask = mask;
}

void RK::set_pot(Potential *ph_pot, Potential *ex_pot){
    _ph_pot = ph_pot;
    _ex_pot = ex_pot;
}

void RK::step(const int ti){
    _update_k(ti);

    cdouble** wfc = _wfc->get_wf();
    cdouble** wfx = _wfx->get_wf();
    cdouble u = 1.0/6.0;

    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            wfc[i][j] += u*(_k1c[i][j] + 2.0*_k2c[i][j] + 2.0*_k3c[i][j] + _k4c[i][j])*_mask->x(i)*_mask->y(j);
            wfx[i][j] += u*(_k1x[i][j] + 2.0*_k2x[i][j] + 2.0*_k3x[i][j] + _k4x[i][j])*_mask->x(i)*_mask->y(j);
        }
    }
}

RK::~RK(){
    free2d<cdouble>(&_k1x,_params->nx,_params->ny);
    free2d<cdouble>(&_k2x,_params->nx,_params->ny);
    free2d<cdouble>(&_k3x,_params->nx,_params->ny);
    free2d<cdouble>(&_k4x,_params->nx,_params->ny);
    free2d<cdouble>(&_k1c,_params->nx,_params->ny);
    free2d<cdouble>(&_k2c,_params->nx,_params->ny);
    free2d<cdouble>(&_k3c,_params->nx,_params->ny);
    free2d<cdouble>(&_k4c,_params->nx,_params->ny);
}
