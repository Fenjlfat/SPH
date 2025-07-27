#ifndef SPH_PARTICLE_H
#define SPH_PARTICLE_H
#include "header.h"
struct particles
{
    //kordinate
    double IX   = 0.0;
    double IY   = 0.0;
    double IZ   = 0.0;
    //velocity
    double IVX  = 0.0;
    double IVY  = 0.0;
    double IVZ  = 0.0;

    std::vector<double*> IXX_Ptr{&IX, &IY, &IZ};    //integer, dimension(3), parameter :: IXX=[IX,IY,IZ]
    std::vector<double*> IVV_Ptr{&IVX, &IVY, &IVZ};    //integer, dimension(3), parameter :: IVV=[IVX,IVY,IVZ]

    double IMAS = 0.0;
    double IDNS = 0.0;   //density
    double IDDNS = 0.0;   //time derivative of density
    double IU   = 0.0;   //internal energy
    double IDU   = 0.0;   //time derivative of internal energy
    double IP   = 0.0;   //pressure
    double ISUB = 0.0;   //substance identifier
    double ICS  = 0.0;   //sound velocity
    double IT    = 0.0;
    double IG  = 0.0;
    double IKK = 0.0;     
    double ICP  = 0.0;

    double IHS   = 0.0;   //kernel scale 

    std::vector<double> IACS{0.0, 0.0, 0.0};    //integer, dimension(3), parameter :: IACS=[13,14,15]
    std::vector<double> IDX{0.0, 0.0, 0.0};    //integer, dimension(3), parameter :: IDX =[19,20,21]

    double IUXY  = 0.0;   
    double IUXZ  = 0.0;   
    double IUYZ  = 0.0;   
    double IUXX  = 0.0;   
    double IUYY  = 0.0;   
    double IUZZ  = 0.0;    
    
    double IULL  = 0.0;     
    //stress
    double ISXY = 0.0;
    double ISXZ = 0.0;
    double ISYZ = 0.0;
    double ISXX = 0.0;
    double ISYY = 0.0;
    double ISZZ = 0.0;

    double IYY = 0.0;
    double RHOI = 0.0;
    double RHOD = 0.0;
    double RHOO = 0.0;

    double IWXX = 0.0;
    double IWYY = 0.0;
    double IWZZ = 0.0;
    double IWXY = 0.0;
    double IWXZ = 0.0;
    double IWYZ = 0.0;
    double IRXY = 0.0;
    double IRXZ = 0.0;
    double IRYZ = 0.0;

    double ICV = 0.0;
    double IETA1 = 0.0;
    double IMisSTR = 0.0;
    double IVREALX  = 0.0;   //Vx  physical
    double IVREALY  = 0.0;   //Vx  physical
    double IVREALZ  = 0.0;   //Vx  physical
    std::vector<double*> IVREALV_Ptr{&IVREALX, &IVREALY, &IVREALZ};   //integer, dimension(3), parameter :: IVREALV=[IVREALX,IVREALY,IVREALZ]
};
#endif // !SPH_PARTICLE_H

