#ifndef RK_H
#define RK_H

#include "param.h"
#include "grid.h"
#include "wavefunction.h"
#include "field.h"
#include "mask.h"
#include "potential.h"


class RK{
    private:
        Params* _params;
        Grid* _grid;
        Field* _field;
        Wavefunction* _wf;
        Mask* _mask;
        Potential *_pot;
        
        //TODO: Store kn matrices in vector containers
        cdouble **_k1, **_k2, **_k3, **_k4;

        cdouble _rhs_f(const int i, const int j, const int ti, cdouble** kn, double pre);
        
        void _update_k(const int ti);
        void _update_k1(const int ti);
        void _update_k2(const int ti);
        void _update_k3(const int ti);
        void _update_k4(const int ti);
        cdouble _laplacian(const int i, const int j, Wavefunction* wf, cdouble** kn, double pre);
	
    public:
        RK(Params* params);
        void set_grid(Grid* grid);
        void set_wf(Wavefunction* wfc);
        void set_field(Field* field);
        void set_mask(Mask *mask);
        void set_pot(Potential *pot);
        void step(const int ti);
        ~RK();


};

#endif
