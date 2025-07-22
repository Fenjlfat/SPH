#include "header.h"
#include "parametrs_modeling.h"
#include "particles.h";
void SPOUT(std::vector<particles> &particle)
{
    parametrs_modeling PAR;
    std::fstream F;
    F.open(PAR.dat_file);
    for (int i = 0; i < particle.size(); i++)
    {
        std::cout << PAR.dat_file << std::endl;
        F << particle[i].IX << particle[i].IY <<std::endl;
        /* code */
    }
    F.close();
}