#ifndef PARAMETRS_MODELING_H
#define PARAMETRS_MODELING_H
#include "header.h"
struct parametrs_modeling
{
    //constants
    const double C1D2 = 1.0 / 2.0;
    const double C1D3 = 1.0 / 3.0;
    const double C2D3 = 2.0 / 3.0;
    const double C4D3 = 1.0 / 2.0;

    const double pi = 3.14159265358979323846264;
    const double DNA = 6.022045000e23;
    const double kB = 1.380622e-23;
    const double hP = 6.6260755e-34;
    const double nuP = 0.38e0; //0.38d0 !0.38d0 !Poisson's coefficient  AL==0.34 ;

    std::string dat_file = {"SPH_CU_R10.dat"};
    std::string dump_file = {"SPH_CU_R10.dump"};
    
    double  dbp = 0.1e-3; //distance_between_particles
    int npaa = 15; //number of particles along axis

    /* data */
};


#endif // !PARAMETRS_MODELING_H
