#include <iostream>
#include <string>
#include <sstream>
#include "param.h"

Params::Params(){
    n_threads   = 8;
    g           = 2.0;
    omega_r     = 4.4;
    gamma_x     = 0.01;
    gamma_c     = 0.01;
    mc          = 2E-5*m0;
    mx          = 0.5*m0;
    xmin        = -60.0;
    xmax        =  60.0;
    ymin        = -60.0;
    ymax        =  60.0;
    nx          = 500;
    ny          = 500;
    dx          = (xmax-xmin)/(double)(nx-1);
    dy          = (ymax-ymin)/(double)(ny-1);
    tmin        = 0;
    tmax        = 20;
    nt          = 50000;
    dt          = (tmax-tmin)/(double)(nt-1);
    nfields     = 1;
    fields_def  = "2.0, 5.0, 10.0, 0.0, 0.0, 1.4, 0.0, -1.0, 0.0, 0, 0;";
                  //"2.0, 5.0,-10.0, 0.0, 0.0, 1.4, 0.0, 0.0, 0.0, 3, 0";

}

void Params::print_param(){
    std::string params;
    std::stringstream ss;
    ss<<"PolaritonBEC:\n"<<
        "----------------------\n"<<
        "g:         "<<g<<std::endl<<
        "omega_r:   "<<omega_r<<std::endl<<
        "gamma_x:   "<<gamma_x<<std::endl<<
        "gamma_c:   "<<gamma_c<<std::endl<<
        "mc:        "<<mc<<std::endl<<
        "xmin:      "<<xmin<<std::endl<<
        "xmax:      "<<xmax<<std::endl<<
        "ymin:      "<<ymin<<std::endl<<
        "ymax:      "<<ymax<<std::endl<<
        "nx:        "<<nx<<std::endl<<
        "ny:        "<<ny<<std::endl<<
        "dx:        "<<dx<<std::endl<<
        "dy:        "<<dy<<std::endl<<
        "tmin:      "<<tmin<<std::endl<<
        "tmax:      "<<tmax<<std::endl<<
        "nt:        "<<nt<<std::endl<<
        "nfields:   "<<nfields<<std::endl<<
        "field_def: "<<fields_def<<std::endl;

    std::cout<<ss.str();
}
