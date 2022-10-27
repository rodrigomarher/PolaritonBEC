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
        Wavefunction* _wfx;
        Wavefunction* _wfc;
        Mask* _mask;
        Potential *_ph_pot, *_ex_pot;
        
        //TODO: Store kn matrices in vector containers
        cdouble **_k1x, **_k2x, **_k3x, **_k4x;
        cdouble **_k1c, **_k2c, **_k3c, **_k4c;

        cdouble _rhs_fc(const int i, const int j, const int ti, cdouble** knc, cdouble** knx, double pre);
        cdouble _rhs_fx(const int i, const int j, const int ti, cdouble** knc, cdouble** knx, double pre);
        
        void _update_k(const int ti);
        void _update_k1(const int ti);
        void _update_k2(const int ti);
        void _update_k3(const int ti);
        void _update_k4(const int ti);
        cdouble _laplacian(const int i, const int j, Wavefunction* wf, cdouble** kn, double pre);
	
    public:
        RK(Params* params);
        void set_grid(Grid* grid);
        void set_wf(Wavefunction* wfx, Wavefunction* wfc);
        void set_field(Field* field);
        void set_mask(Mask *mask);
        void set_pot(Potential *ph_pot, Potential *ex_pot);
        void step(const int ti);
        ~RK();


};

#endif
