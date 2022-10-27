#include <iostream>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include "polariton.h"

Polariton::Polariton(Params *params){
    _params = params;
    _grid = new Grid(_params); 

    _wfc = new Wavefunction(_params);
    _wfx = new Wavefunction(_params);
    _wfc->set_grid(_grid);
    _wfx->set_grid(_grid);
    _wfx->init_zeros();
    _wfc->init_zeros();

    _field = new Field(_params);
    _field->set_grid(_grid);
    _field->set_field();

    _ph_pot = new Potential(_params);
    _ex_pot = new Potential(_params);
    _ph_pot->set_grid(_grid);
    _ex_pot->set_grid(_grid);
    _ph_pot->set_def(_params->ph_pot);
    _ex_pot->set_def(_params->ex_pot);
    _ph_pot->set_potential();
    _ex_pot->set_potential();

    _mask = new Mask(_params);
    _mask->set_grid(_grid);
    _mask->set_mask();

    _rk = new RK(_params);
    _rk->set_grid(_grid);
    _rk->set_wf(_wfx, _wfc);
    _rk->set_field(_field);
    _rk->set_mask(_mask);
    _rk->set_pot(_ph_pot, _ex_pot);

    _params->print_param();
    
    std::cout<<"Photonic wf norm: "<<_wfc->norm()<<std::endl;
    std::cout<<"Excitonic wf norm: "<<_wfx->norm()<<std::endl;
     
    omp_set_num_threads(_params->n_threads);
    //_field->update(0);
    //_field->save_field("results/pump_field.dat");
};

void Polariton::set_param(Params *params){
    _params = params;
}

void Polariton::evolve(){
    for(int ti=0; ti<_params->nt;ti++){
        _field->update(ti);
        _rk->step(ti);
        if(ti%100==0){
            std::cout<<"ti: "<<ti;
            std::cout<<" Photonic wf norm: "<<_wfc->norm();
            std::cout<<" Excitonic wf norm: "<<_wfx->norm()<<std::endl;
            snapshot(ti);
        }
    }
}

void Polariton::snapshot(const int ti){
    std::stringstream ss;
    ss<<std::setw(8)<<std::setfill('0')<<ti;
    std::string name_c = "results/wfc_"+ss.str()+".dat";
    std::string name_x = "results/wfx_"+ss.str()+".dat";
    std::string name_p = "results/pump_"+ss.str()+".dat";
    _wfc->save_wf(name_c); 
    _wfx->save_wf(name_x);
    _field->save_field(name_p);
}
