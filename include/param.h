#ifndef PARAM_H
#define PARAM_H

#include <cmath>
#include <string>

#define cdouble std::complex<double>
#define I 	std::complex<double>(0.0,1.0)

#define hbar    0.6582  //Reduced Planck constant [meV ps]
#define m0      5.11E3 //Electron's rest mass [meV um^-2 ps^2]

class Params{
    public:
        int n_threads;          //Number of threads for parallelization
        double g;               //Interaction coefficient [meV um^2]
        double omega_r;         //Rabi frequency [meV]
        double omega_c;         //Cavity photon frequency [meV]
        double omega_x;         //Cavity exciton frequency [meV] 
        double gamma_x;         //Exciton decay rate [1/ps] 
        double gamma_c;         //Photon decay rate [1/ps] 
        double mc;              //Polariton effective mass [m0]
        double mx;              //Polariton effective mass [m0]
        double xmin;            //Min x simulation box [um] 
        double xmax;            //Max x simulation box [um]
        double ymin;            //Min y simulation box [um]
        double ymax;            //Max y simulation box [um]
        int nx;                 //Number points along x
        int ny;                 //Number points along y
        double dx;              //Spatial step x [um]
        double dy;              //Spatial step y [um]
        double tmin;            //Min simulation time [ps]
        double tmax;            //Max simulation time [ps] 
        int nt;                 //Number temporal points
        double dt;              //Temporal step [ps]
        int ndump;
        int nbuf;
        double mask_xmax;       //Mask max x [um]
        double mask_ymax;       //Mask max y [um]
        double mask_gamma;      //Mask gamma value
        int nfields;            //Number of fields
        std::string fields_def; //Field defunition
        std::string ph_pot;     //Photonic potential
        std::string ex_pot;     //Excitonic potential
        
	Params();
        void print_param();
};

#endif
