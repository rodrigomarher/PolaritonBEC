#include "mask.h"

Mask::Mask(Params *param){
    _param = param;
    _xmask = alloc1d<cdouble>(_param->nx);
    _ymask = alloc1d<cdouble>(_param->ny);
}

void Mask::set_grid(Grid *grid){
    _grid = grid;
}

void Mask::set_mask(){
    double xb = _param->mask_xmax;
    double yb = _param->mask_ymax;
    double gamma = _param->mask_gamma;

    for(int i=0; i<_param->nx; i++){
        _xmask[i] = cdouble(1.0,0.0);
        if(_grid->x(i)<_grid->x(0)+xb){
            _xmask[i] = pow(cos(M_PI*(_grid->x(i)-(_grid->x(0)+xb))*gamma/(2.0*xb)),1.0/8.0);
            if(_xmask[i].real()<0.0){
                _xmask[i] = 0.0;
            }
        }
        else if( _grid->x(i)>_grid->x(_param->nx-1)-xb){
            _xmask[i] = pow(cos(M_PI*(_grid->x(i)-(_grid->x(_param->nx-1)-xb))*gamma/(2.0*xb)),1.0/8.0);
            if(_xmask[i].real()<0.0){
                _xmask[i] = 0.0;
            }
        }
    }

    for(int j=0; j<_param->ny; j++){
        _ymask[j] = cdouble(1.0,0.0);
        if(_grid->y(j)<_grid->y(0)+yb){
            _ymask[j] = pow(cos(M_PI*(_grid->y(j)-(_grid->y(0)+yb))*gamma/(2.0*xb)),1.0/8.0);
            if(_ymask[j].real()<0.0){
                _ymask[j] = 0.0;
            }
        }
        else if( _grid->y(j)>_grid->y(_param->ny-1)-yb){
            _ymask[j] = pow(cos(M_PI*(_grid->y(j)-(_grid->y(_param->ny-1)-xb))*gamma/(2.0*yb)),1.0/8.0);
            if(_ymask[j].real()<0.0){
                _ymask[j] = 0.0;
            }
        }
    }
}

cdouble Mask::x(const int i){
    return _xmask[i];
} 

cdouble Mask::y(const int j){
    return _ymask[j];
}

Mask::~Mask(){
    free1d<cdouble>(&_xmask,_param->nx);
    free1d<cdouble>(&_ymask,_param->ny);
}
