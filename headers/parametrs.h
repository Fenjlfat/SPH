#ifndef PARAMETRS_H
#define PARAMETRS_H
#include "header.h"
struct parametrs
{
    //constants
    const double C1D2 = 1.0 / 2.0;
    const double C1D3 = 1.0 / 3.0;
    const double C2D3 = 2.0 / 3.0;
    const double C4D3 = 4.0 / 3.0;
    const double CC12 = 24.0 * sqrt(3.0);

    // Константы и глобальные параметры
    const double pi = 3.14159265358979323846264;
    const double DNA = 6.022045000e23;  //avogadro
    const double kB = 1.380622e-23;     //boltzman
    const double hP = 6.6260755e-34;
    const double nuP = 0.38e0; //0.38d0 !0.38d0 !Poisson's coefficient  AL==0.34 ;
    
    const double AMOL = 64.685;       // кг/моль (для Cu)
    const double EpsL = 8.0 * 1.6e-19; // эВ/б
    const double VIMB = 0.0642;       // скорость иммобилизации
    const double RHOD0 = 1.0e11;      // 1/m²
    const double DKA = 3.392;         // коэффициент аннигиляции
    const double C_FCC = 0.5 * sqrt(2.0) * pow(4.0, 1.0/3.0);
    
    std::string dat_file = {"/mnt/disk1/LINUX/SPH/PROGRAMMS/C++/SPH/SPH_CU_R10.dat"};
    std::string dump_file = {"/mnt/disk1/LINUX/SPH/PROGRAMMS/C++/SPH/SPH_CU_R10.dump"};
    
    //======================================================system
    double scale = 1.00;
    double  dbp = 1.0e-3; //distance_between_particles
    int npaax = 5; //number of particles along axis
    int npaay = 5; //number of particles along axis
    int npaaz = 5; //number of particles along axis

    //======================================================dislocations
    int DP0  =  9;  //number of variables
    int NSS0  =  12;  //number of groups of dislocations
    //===============variables
    int JRHOD = 0; 
    int JRHOI = 1;
    int JVD = 2; 
    int DNBXX = 3;
    int DNBYY = 4;
    int DNBZZ = 5;
    int DNBXY = 6;
    int DNBXZ = 7;
    int DNBYZ = 8;
    //======================================================
};


#endif // !PARAMETRS_MODELING_H
