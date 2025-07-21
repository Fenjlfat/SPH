#include <iostream>
#include "particles.h"
int main() 
{
    //particle PART;
    std::vector<particles> particle(10); 
    for (int i = 0; i < 10; i++)
    {
        particle[i].IX = i;
        particle[i].IY = static_cast<double> (i);
        particle[i].IZ = static_cast<double> (i);
        particle[i].IVX = static_cast<double> (i);
        particle[i].IVY = static_cast<double> (i);
        particle[i].IVZ = static_cast<double> (i);
    }
    for (int i = 0; i < 10; i++)
    {
        std::cout << "координата X=" << particle[i].IX << std::endl;
        std::cout << "ссылка на X=" << *particle[i].IXX_Ptr[0] << std::endl;
    }
    
    *particle[0].IXX_Ptr[0] = 345;
    *particle[2].IXX_Ptr[0] = 555;
    std::cout << "0_IX = " << particle[0].IX << std::endl;
    std::cout << "2_IX = " << particle[2].IX << std::endl;
    std::cout << "SPH" << std::endl;
    return 0;
}
