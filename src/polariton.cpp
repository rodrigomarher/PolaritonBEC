#include <iostream>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <filesystem>
#include "polariton.h"

Polariton::Polariton(Params *params){
    _params = params;
    _grid = new Grid(_params); 

    _wf = new Wavefunction(_params);
    _wf->set_grid(_grid);
    _wf->init_zeros();

    _field = new Field(_params);
    _field->set_grid(_grid);
    _field->set_field();

    _pot = new Potential(_params);
    _pot->set_grid(_grid);
    _pot->set_def(_params->ph_pot);
    _pot->set_potential();

    _mask = new Mask(_params);
    _mask->set_grid(_grid);
    _mask->set_mask();
    _mask->save_mask();

    _rk = new RK(_params);
    _rk->set_grid(_grid);
    _rk->set_wf(_wf);
    _rk->set_field(_field);
    _rk->set_mask(_mask);
    _rk->set_pot(_pot);

    _params->print_param();
    if(!std::filesystem::exists(_params->out_folder)){
        std::filesystem::create_directory(_params->out_folder);
    }
     
    omp_set_num_threads(_params->n_threads);
    //_field->update(0);
    //_field->save_field("results/pump_field.dat");
};

void Polariton::set_param(Params *params){
    _params = params;
}

void Polariton::evolve(){
    int counter_buffer = 0;
    cdouble **wf_avg;
    wf_avg = alloc2d<cdouble>(_params->nx, _params->ny);
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            wf_avg[i][j] = cdouble(0.0,0.0);
        }
    }

    for(int ti=0; ti<_params->nt;ti++){
        _field->update(ti);
        _rk->step(ti);
    
        #pragma omp parallel for collapse(1) schedule(dynamic)
        for(int i=0; i<_params->nx; i++){
            for(int j=0; j<_params->ny; j++){
                wf_avg[i][j] += _wf->at(i,j)*_params->dt;
            }
        }

        if(ti%_params->ndump==0){
            std::cout<<(counter_buffer)%_params->nbuf<<" "<<counter_buffer<<std::endl;
            _wf->copy_to_buf((counter_buffer)%_params->nbuf);
            _field->copy_to_buf((counter_buffer)%_params->nbuf);
            std::cout<<"ti: "<<ti;
            std::cout<<" WF norm: "<<_wf->norm();
            //snapshot(ti);
            counter_buffer++;

            if(counter_buffer%_params->nbuf==0){
                std::cout<<"Doing snapshot..."<<std::endl;
                snapshot(ti);
            }
        }
    }
    std::string path = _params->out_folder+"/wf_avg.dat";
    write_array2d_complex(wf_avg, _params->nx, _params->ny, path);
    path= _params->out_folder+"/wf_lst.dat";
    write_array2d_complex(_wf->get_wf(),_params->nx, _params->ny, path);
}

void Polariton::snapshot(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nbuf; i++){
        int number = ti - _params->ndump*(_params->nbuf-i-1);
        std::stringstream ss;
        ss<<std::setw(8)<<std::setfill('0')<<number;
        std::string name_c = _params->out_folder + "/wf_" + ss.str() + ".dat";
        std::string name_p = _params->out_folder + "/pump_" + ss.str() + ".dat";
        _wf->save_wf(name_c,i); 
        _field->save_field(name_p,i);
    }
}
