#include <iostream>
#include "subfield.h"

subField::subField(std::string def){
    _definition = def;
    _parse();
}

void subField::_parse(){
    char delim = ',';
    int counter = 0;
    std::stringstream ss(_definition);
    std::string item;
    std::vector<std::string> field_parameters;

    while(getline(ss,item,delim)){
        field_parameters.push_back(item);
        counter++;
    }
    if(counter!=12){
        std::cout<<"[ERROR] Invalid number of field parameters."<<std::endl;
    }
    _fp = std::stod(field_parameters[0]);
    _w0 = std::stod(field_parameters[1]);
    _x0 = std::stod(field_parameters[2]);
    _y0 = std::stod(field_parameters[3]);
    _tmin = std::stod(field_parameters[4]);
    _tmax = std::stod(field_parameters[5]);
    _omega = std::stod(field_parameters[6]);
    _kpx = std::stod(field_parameters[7]);
    _kpy = std::stod(field_parameters[8]);
    _l = std::stoi(field_parameters[9]);
    _p = std::stoi(field_parameters[10]);
    _cep = std::stod(field_parameters[11])*M_PI;
    std::cout<<"CEP = "<<_cep<<std::endl;

}

cdouble subField::_gaussian(const double x, const double y, const cdouble tenv){
    double r2 = (x-_x0)*(x-_x0) + (y-_y0)*(y-_y0);
    double w02 = _w0*_w0;

    return tenv*_fp*exp(-r2/w02);
}

cdouble subField::_laguerre(const double x, const double y, const cdouble tenv){
    double w02 = _w0*_w0;
    double sqrt2 = sqrt(2.0);
    double norm = sqrt(2*factorial(_p)/(M_PI*factorial(_l+_p))) ;
    double r2,theta,a,b,c;
    cdouble phase;

    r2 = (x-_x0)*(x-_x0)+(y-_y0)*(y-_y0);
    theta = atan2((y-_y0),(x-_x0));
    a = pow(sqrt(r2)*sqrt2/_w0,abs(_l));
    b = exp(-r2/w02);
    c = std::assoc_laguerre(_p,_l,2.0*r2/w02);
    phase = exp(-I*(double)_l*theta+I*(_kpx*x + _kpy*y));
    return tenv*_fp*(double)norm*a*b*c*phase;
}

cdouble subField::_tenv_sin2(const int ti){
    double deltat = _tmax-_tmin; 
    cdouble phase = exp(-I*_omega/hbar*_grid->t(ti) + _cep);
    if(_grid->t(ti) < _tmin){
        return cdouble(0.0,0.0);
    }
    if(_grid->t(ti) > _tmax){
        return cdouble(0.0,0.0);
        }
    cdouble env = cdouble(pow(sin(M_PI*(_grid->t(ti)-_tmin)/deltat),2),0.0)*phase;
    return env;
}

void subField::set_grid(Grid *grid){
    _grid = grid;
}

cdouble subField::at(const int i, const int j, const int ti){
    cdouble env = _tenv_sin2(ti);
    double x = _grid->x(i);
    double y = _grid->y(j);
    return _laguerre(x,y,env);
}


