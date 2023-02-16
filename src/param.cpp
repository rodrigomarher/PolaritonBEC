#include <iostream>
#include <string>
#include <sstream>
#include "param.h"

Params::Params(){
    n_threads   = 16;
    g           = 0.01;
    omega_r     = 5.1;
    omega_c     = 1483.1;
    omega_x     = 1483.1;
    gamma_x     = 0.01;
    gamma_c     = 0.1;
    mc          = 2E-5*m0;
    mx          = 0.5*m0;
    xmin        = -60.0;
    xmax        =  60.0;
    ymin        = -60.0;
    ymax        =  60.0;
    nx          = 1024;
    ny          = 1024;
    dx          = (xmax-xmin)/(double)(nx-1);
    dy          = (ymax-ymin)/(double)(ny-1);
    tmin        = 0;
    tmax        = 2.0;
    nt          = 25000*2;
    dt          = (tmax-tmin)/(double)(nt-1);
    ndump       = 100;
    nbuf        = 250;
    mask_xmax   = 20.0;
    mask_ymax   = 20.0;
    mask_gamma  = 1.2;
    nfields     = 2;
    fields_def  = "10000.0, 2.0, 5.0, 0.0, 0.0, 0.0055731999999999995, 1483.1, -3.0, 0.0, 0, 0, 0.0;"
                  "10000.0, 2.0, 0.0, 0.0, 0.0, 0.0055731999999999995, 1483.1, 1.0, 0.0, 0, 0, 0.0;";
    ph_pot      = "cylinder, 0.0, 4.0, 0.0, 0.0;";
    ex_pot      = ",";
    out_folder  = "results_free_2pump";
}

void Params::print_param(){
    std::string params;
    std::stringstream ss;
    ss<<"PolaritonBEC:\n"<<
        "----------------------\n"<<
        "g:             "<<g<<std::endl<<
        "omega_r:       "<<omega_r<<std::endl<<
        "gamma_x:       "<<gamma_x<<std::endl<<
        "gamma_c:       "<<gamma_c<<std::endl<<
        "mc:            "<<mc<<std::endl<<
        "xmin:          "<<xmin<<std::endl<<
        "xmax:          "<<xmax<<std::endl<<
        "ymin:          "<<ymin<<std::endl<<
        "ymax:          "<<ymax<<std::endl<<
        "nx:            "<<nx<<std::endl<<
        "ny:            "<<ny<<std::endl<<
        "dx:            "<<dx<<std::endl<<
        "dy:            "<<dy<<std::endl<<
        "tmin:          "<<tmin<<std::endl<<
        "tmax:          "<<tmax<<std::endl<<
        "nt:            "<<nt<<std::endl<<
        "ndump:         "<<ndump<<std::endl<<
        "nbuf:          "<<nbuf<<std::endl<<
        "mask_xmax:     "<<mask_xmax<<std::endl<<
        "mask_ymax:     "<<mask_ymax<<std::endl<<
        "mask_gamma:    "<<mask_gamma<<std::endl<<
        "nfields:       "<<nfields<<std::endl<<
        "field_def:     "<<fields_def<<std::endl<<
        "ph_pot:        "<<ph_pot<<std::endl<<
        "ex_pot:        "<<ex_pot<<std::endl;

    std::cout<<ss.str();
}
