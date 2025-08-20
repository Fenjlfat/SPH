#include "header.h"
#include "parametrs.h"
#include "particles.h"
void SPOUT(int NTS, int NPT, std::vector<particles> &particle)
{
    parametrs par;
    
    double XMAX = particle[0].IX;
    double XMIN = particle[0].IX;
    double YMAX = particle[0].IY;
    double YMIN = particle[0].IY;
    double ZMAX = particle[0].IZ;
    double ZMIN = particle[0].IZ;

    for (const auto &p : particle) 
    {
        if (XMIN > p.IX) XMIN = p.IX;
        if (XMAX < p.IX) XMAX = p.IX;
        if (YMIN > p.IY) YMIN = p.IY;
        if (YMAX < p.IY) YMAX = p.IY;
        if (ZMIN > p.IZ) ZMIN = p.IZ;
        if (ZMAX < p.IZ) ZMAX = p.IZ;
    }
    
    std::fstream FILE;
    FILE.open(par.dump_file, std::ios::app); //dozapis
    //FILE.open(par.dump_file);
    if (!FILE.is_open())
    {
        std::cout << "file is not open" << std::endl;
    }
    else
    {
        FILE << "ITEM: TIMESTEP" << std::endl;
        FILE << "   " << NTS << std::endl;
        FILE << "ITEM: NUMBER OF ATOMS" << std::endl;
        FILE << "   " << NPT << std::endl;
        FILE << "ITEM: BOX BOUNDS pp pp pp" << std::endl;
        FILE << "   " << XMIN * par.scale << "   " << XMAX * par.scale << std::endl;
        FILE << "   " << YMIN * par.scale << "   " << YMAX * par.scale << std::endl;
        FILE << "   " << ZMIN * par.scale << "   " << ZMAX * par.scale << std::endl;
        //FILE << "ITEM: ATOMS id type x y z "<< std::endl;
        FILE << "ITEM: ATOMS id type x y z vx vy vz c_dns c_pres c_t Sxx Syy Szz Sxy Sxz Syz uxx uyy uzz uxy uxz uyz wxx wyy wzz wxy wxz wyz RhoD RhoI YY Si Wi Ui" << std::endl;
        //Title="ITEM: ATOMS id type x y z vx vy vz c_dns c_pres c_t Sxx Syy Szz Sxy Sxz Syz"//&
        // " uxx uyy uzz uxy uxz uyz wxx wyy wzz wxy wxz wyz RhoD RhoI YY Si Wi Ui"
        //realization otput to file
        int counter = 0;
        for (const auto &p : particle)
        {
            double STRI = sqrt(pow((p.ISXX - p.ISYY), 2) + pow((p.ISYY - p.ISZZ), 2) + pow((p.ISZZ - p.ISXX), 2) +
                6.e0 * (pow(p.ISXY, 2) + pow(p.ISXZ, 2) + pow(p.ISYZ,2))) * 1.e-9 / sqrt(2.e0);

            double WPLI = sqrt(pow((p.IWXX - p.IWYY), 2) + pow((p.IWYY - p.IWZZ), 2) + pow((p.IWZZ - p.IWXX), 2) +
                1.5e0 * (pow(p.IWXY, 2) + pow(p.IWXZ, 2) + pow(p.IWYZ, 2))) * 1.e0 / sqrt(2.e0);

            double UMCI = sqrt(pow((p.IUXX - p.IUYY), 2) + pow((p.IUYY - p.IUZZ), 2) + pow((p.IUZZ - p.IUXX), 2) +
                1.5e0 * (pow(p.IUXY, 2) + pow(p.IUXZ, 2) + pow(p.IUYZ, 2))) * 1.e0 / sqrt(2.e0);

            FILE << "   " << counter << "   " << 1 
                << "    " << p.IX * par.scale << "   " << p.IY * par.scale << "   " << p.IZ * par.scale 
                << "    " << p.IVX << "    " << p.IVY << "    " << p.IVZ
                << "    " << p.IDNS << "    " << p.IP * 1.e-9 << "    " << p.IT 
                << "    " << p.ISXX * 1.e-9 << "    " << p.ISYY * 1.e-9 << "    " << p.ISZZ * 1.e-9
                << "    " << p.ISXY * 1.e-9 << "    " << p.ISXZ * 1.e-9 << "    " << p.ISYZ * 1.e-9 
                << "    " << p.IUXX << "    " << p.IUYY << "    " << p.IUZZ 
                << "    " << p.IUXY << "    " << p.IUXZ << "    " << p.IUYZ 
                << "    " << p.IWXX << "    " << p.IWYY << "    " << p.IWZZ 
                << "    " << p.IWXY << "    " << p.IWXZ << "    " << p.IWYZ 
                << "    " << p.RHOD << "    " << p.RHOI << "    " << p.IYY * 1.e-9 
                << "    " << p.IMisSTR*1.e-9 << "    " << WPLI << "    " << UMCI <<std::endl;
            
            counter++;
       }    
    } 
    FILE.close();
}