#include <iostream>
#include "particle.h"
int main() 
{
    //particle PART;
    std::vector<particle> PART(10); 
    for (int i = 0; i < 10; i++)
    {
        PART[i].IX = i;
        PART[i].IY = static_cast<double> (i);
        PART[i].IZ = static_cast<double> (i);
        PART[i].IVX = static_cast<double> (i);
        PART[i].IVY = static_cast<double> (i);
        PART[i].IVZ = static_cast<double> (i);
    }
    for (int i = 0; i < 10; i++)
    {
        std::cout << "координата X=" << PART[i].IX << std::endl;
        std::cout << "ссылка на X=" << *PART[i].IXX_Ptr[0] << std::endl;
    }
    
    *PART[0].IXX_Ptr[0] = 345;
    *PART[2].IXX_Ptr[0] = 555;
    std::cout << "0_IX = " << PART[0].IX << std::endl;
    std::cout << "2_IX = " << PART[2].IX << std::endl;
    std::cout << "SPH" << std::endl;
    return 0;
}
