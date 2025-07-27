#include <iostream>
#include "header.h"
#include "particles.h"
#include "parametrs.h"

void INIT(std::vector<particles> &particle);
void MOVE(std::vector<particles> &particle);
void SPOUT(int NTS, int NPT, std::vector<particles> &particle);


int main() 
{
    parametrs par;
    int NPT = par.npaax * par.npaay * par.npaaz; // number of particle
    
    //create vector of structure that contains value of particles
    std::vector<particles> particle(NPT);
    
    //initialization of vector particle
    INIT(particle);
    
    for (int NTS = 0; NTS < 1; NTS++)
    {
        MOVE(particle);
        SPOUT(NTS, NPT, particle);
    }
    
    std::cout << "well done!" << std::endl;
    return 0;
}
