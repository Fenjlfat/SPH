#include <iostream>
#include "header.h"
#include "particles.h"
#include "parametrs.h"

void EOS(double &E_EOS, double &P_EOS, double &T_EOS, double &CV_EOS, double &CS_EOS, double &DNS_EOS);
void INIT(double &E_EOS, double &P_EOS, double &T_EOS, double &CV_EOS, double &CS_EOS, double &DNS_EOS, std::vector<particles> &particle);
void NLIST(std::vector<particles> &particle, parametrs &parametr); 
void MOVE(double &TIME, double &DTAU, std::vector<particles> &particle, parametrs &parametr);
void DISLOC(double DTAU, std::vector<particles>& particle, parametrs &parametr);
void SPOUT(int NTS, int NPT, std::vector<particles> &particle);


int main() 
{
    parametrs parametr;
    int NPT = parametr.npaax * parametr.npaay * parametr.npaaz; // number of particle

    //create vector of structure that contains value of particles
    std::vector<particles> particle(NPT);

    double E_EOS = 0.0, P_EOS = 0.0, T_EOS = 0.0, CV_EOS = 0.0, CS_EOS = 0.0, DNS_EOS = 0.0;
    EOS(E_EOS, P_EOS, T_EOS, CV_EOS, CS_EOS, DNS_EOS); 
    //initialization of vector particle
    INIT(E_EOS, P_EOS, T_EOS, CV_EOS, CS_EOS, DNS_EOS, particle);
    
    for (int NTS = 0; NTS < 1; NTS++)
    {
        NLIST(particle, parametr);
        MOVE(TIME, DTAU, particle, parametr);
        DISLOC(DTAU, particle, parametr);
        SPOUT(NTS, NPT, particle);
    }
    
    std::cout << "well done!" << std::endl;
    return 0;
}
