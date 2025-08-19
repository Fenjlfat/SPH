#include <iostream>
#include "header.h"
#include "particles.h"
#include "parametrs.h"


void INIT(std::vector<particles> &particle, parametrs &parametr);
void NLIST(std::vector<particles> &particle, parametrs &parametr); 
void MOVE(std::vector<particles> &particle, parametrs &parametr);
void DISLOC(std::vector<particles> &particle, parametrs &parametr);
void SPOUT(int NTS, int NPT, std::vector<particles> &particle);


int main() 
{
    parametrs parametr;
    int NPT = parametr.npaax * parametr.npaay * parametr.npaaz; // number of particle

    //create vector of structure that contains value of particles
    std::vector<particles> particle(NPT);

    //initialization of vector particle
    INIT(particle, parametr);
    
    for (int NTS = 0; NTS < 1500; NTS++)
    {
        if (NTS == 0 || NTS % 10 == 0) NLIST(particle, parametr);
        MOVE(particle, parametr);
        DISLOC(particle, parametr);
        //if (NTS == 0 || NTS % 10 == 0) SPOUT(NTS, NPT, particle);
        SPOUT(NTS, NPT, particle);
        std::cout << "NTS=" << NTS << "\n";
    }
    
    std::cout << "well done!" << std::endl;
    return 0;
}
