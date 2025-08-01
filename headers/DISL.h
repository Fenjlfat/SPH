#include "header.h"
struct dislocations
{
    int DP0  =  9;  //number of variables
    int NSS0  =  12;  //number of groups of dislocations
    //======================================================
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
    //std::vector<std::vector<std::vector<double>>> DISL(NPT, std::vector<std::vector<double>>(NSS0, std::vector<double>(размер)))
    std::vector<std::vector<double>> DISL;
    dislocations(int NSS0, int DP0, 0.0)
    //std::vector<std::vector<double*>> DISL{{},{},};
};
