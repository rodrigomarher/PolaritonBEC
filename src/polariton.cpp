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
    _mask->save_mask();

    _rk = new RK(_params);
    _rk->set_grid(_grid);
    _rk->set_wf(_wfc, _wfx);
    _rk->set_field(_field);
    _rk->set_mask(_mask);
    _rk->set_pot(_ph_pot, _ex_pot);

    _params->print_param();
    
     
    omp_set_num_threads(_params->n_threads);
    //_field->update(0);
    //_field->save_field("results/pump_field.dat");
};

void Polariton::set_param(Params *params){
    _params = params;
}

void Polariton::evolve(){
    int counter_buffer = 0;
    cdouble **wf_ph_avg, **wf_ex_avg;
    wf_ph_avg = alloc2d<cdouble>(_params->nx, _params->ny);
    wf_ex_avg = alloc2d<cdouble>(_params->nx, _params->ny);
    for(int i=0;i<_params->nx;i++){
        for(int j=0; j>_params->ny;j++){
            wf_ph_avg[i][j] = cdouble(0.0,0.0); 
            wf_ex_avg[i][j] = cdouble(0.0,0.0); 
        }
    }
    for(int ti=0; ti<_params->nt;ti++){
        _field->update(ti);
        _rk->step(ti);
    	

	for(int i=0;i<_params->nx;i++){
	    for(int j=0; j<_params->ny;j++){
	        wf_ph_avg[i][j] += _wfc->at(i,j)*_params->dt; 
	        wf_ex_avg[i][j] += _wfx->at(i,j)*_params->dt; 
	    }
	}


        if(ti%_params->ndump==0){
            std::cout<<(counter_buffer)%_params->nbuf<<" "<<counter_buffer<<std::endl;
            _wfc->copy_to_buf((counter_buffer)%_params->nbuf);
            _wfx->copy_to_buf((counter_buffer)%_params->nbuf);
            _field->copy_to_buf((counter_buffer)%_params->nbuf);
            std::cout<<"ti: "<<ti;
            std::cout<<" Photonic wf norm: "<<_wfc->norm();
            std::cout<<" Excitonic wf norm: "<<_wfx->norm()<<std::endl;
            //snapshot(ti);
            counter_buffer++;

            if(counter_buffer%_params->nbuf==0){
                snapshot(ti);
            }
        }
    }
    std::string path = "results/wfc_avg.dat";
    write_array2d_complex(wf_ph_avg,_params->nx,_params->ny,path);
    path = "results/wfx_avg.dat";
    write_array2d_complex(wf_ex_avg,_params->nx,_params->ny,path);
}

void Polariton::snapshot(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nbuf; i++){
        int number = ti - _params->ndump*(_params->nbuf-i-1);
        std::stringstream ss;
        ss<<std::setw(8)<<std::setfill('0')<<number;
        std::string name_c = "results/wfc_"+ss.str()+".dat";
        std::string name_x = "results/wfx_"+ss.str()+".dat";
        std::string name_p = "results/pump_"+ss.str()+".dat";
        _wfc->save_wf(name_c,i); 
        _wfx->save_wf(name_x,i);
        _field->save_field(name_p,i);
    }
}
