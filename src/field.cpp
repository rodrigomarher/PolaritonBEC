#include <cmath>
#include <iostream>
#include "field.h"


Field::Field(Params *params){
    _params = params;
    _field = alloc2d<cdouble>(_params->nx, _params->ny);
    _fieldbuf = alloc3d<cdouble>(_params->nbuf,_params->nx, _params->ny);
    _parse();
}

void Field::_parse(){
    char delim= ';';
    int counter = 0;
    std::stringstream ss(_params->fields_def);
    std::string item;
    while(getline(ss, item, delim)){
        _def_vec.push_back(item);
        counter++;
    }
    if(counter != _params->nfields){
        std::cout<<"[ERROR] n_fields not equal to number of definitions."<<std::endl;
    }
}

void Field::set_grid(Grid *grid){
    _grid = grid;
}

void Field::set_field(){
    for(int k=0; k<_params->nfields;k++){
        _field_container.push_back(new subField(_def_vec[k]));
        _field_container[k]->set_grid(_grid);
    }
}

cdouble** Field::get_field(){
    return _field;
}

cdouble Field::at(const int i, const int j){
    return _field[i][j];
}

void Field::update(const int ti){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_params->nx; i++){
        for(int j=0; j<_params->ny; j++){
            cdouble field_sum = cdouble(0.0,0.0);
            for(int k=0; k<_params->nfields;k++){
               field_sum += _field_container[k]->at(i, j, ti);
            }
            _field[i][j] = field_sum;
        }
    }
}

void Field::save_field(std::string filename, const int ti){
   write_array2d_complex(_fieldbuf[ti], _params->nx, _params->ny, filename);
}

void Field::copy_to_buf(const int ti){
    for (int i=0;i<_params->nx;i++){
        for(int j=0; j<_params->ny;j++){
            _fieldbuf[ti][i][j] = _field[i][j];
        }
    }
}
Field::~Field(){
    free2d<cdouble>(&_field, _params->nx, _params->ny);
}





