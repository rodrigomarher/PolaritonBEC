#include "param.h"
#include "polariton.h"

int main(){
    Params params;
    Polariton polariton(&params);
    polariton.evolve();
    return 0;
}
