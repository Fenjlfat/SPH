#include "header.h"
#include "particles.h"


int NMES0=50;		//number of meshs
int NNEB0=500;		//maximal number of neighbours
int NPAT[NMES0][NMES0][NMES0][NNEB0];
int MESH [NPT0][3];


void NLIST(std::vector<particles> &particle) 
{
    int NMESX = 0;
    int NMESY = 0;
    int NMESZ = 0;

    double XMAX = particle[0].IX;
    double XMIN = XMAX;
    double YMAX = particle[0].IY;
    double YMIN = YMAX;
    double ZMAX = particle[0].IZ;
    double ZMIN = ZMAX;
    double HMAX = particle[0].IHS;

    // Находим границы области и максимальный масштаб ядра
    for (int I = 1; I < particle.size(); I++) 
    {
        double XX = particle[I].IX;
        double YY = particle[I].IY;
        double ZZ = particle[I].IZ;

        XMIN = std::min(XMIN, XX);
        XMAX = std::max(XMAX, XX);
        YMIN = std::min(YMIN, YY);
        YMAX = std::max(YMAX, YY);
        ZMIN = std::min(ZMIN, ZZ);
        ZMAX = std::max(ZMAX, ZZ);

        HMAX = std::max(HMAX, particle[I].IHS);
    }

    // Проверка на нулевой масштаб ядра
    if (HMAX <= 0.0) 
    {
        std::cerr << "ERROR message from SPHEP_NLIST: ZERO Kernel Scale HMAX" << std::endl;
        return;
    }

    // Вычисление размеров сетки
    NMESX = std::min(static_cast<int>((XMAX - XMIN) / (2.0 * HMAX)), NMES0);
    NMESY = std::min(static_cast<int>((YMAX - YMIN) / (2.0 * HMAX)), NMES0);
    NMESZ = std::min(static_cast<int>((ZMAX - ZMIN) / (2.0 * HMAX)), NMES0);

    if (std::max({NMESX, NMESY, NMESZ}) > NMES0) 
    {
        std::cerr << "ERROR message from SPHEP_NLIST: Increase upto NMES0=" 
                  << std::max({NMESX, NMESY, NMESZ}) << std::endl;
        return;
    }

    double DXX = 1.0 / (XMAX - XMIN);
    double DYY = 1.0 / (YMAX - YMIN);
    double DZZ = 1.0 / (ZMAX - ZMIN);

    // Инициализация NPAT
    NPAT.assign(NMESX + 1, std::vector<std::vector<std::vector<int>>>(NMESY + 1, std::vector<std::vector<int>>(NMESZ + 1, std::vector<int>(NNEB0 + 1, 0))));

    // Заполнение сетки частицами
    for (int I = 0; I < particle.size(); ++I) {
        int NX = static_cast<int>((particle[I].IX - XMIN) * DXX * (NMESX - 1)) + 1;
        int NY = static_cast<int>((particle[I].IY - YMIN) * DYY * (NMESY - 1)) + 1;
        int NZ = static_cast<int>((particle[I].IZ - ZMIN) * DZZ * (NMESZ - 1)) + 1;

        NX = std::min(NX, NMESX);
        NY = std::min(NY, NMESY);
        NZ = std::min(NZ, NMESZ);

        MESH[I][0] = NX; //MESH[I][IX] = NX;
        MESH[I][1] = NY; //MESH[I][IY] = NY;
        MESH[I][2] = NZ; //MESH[I][IZ] = NZ;

        int NUM = NPAT[NX][NY][NZ][0] + 1;
        if (NUM > NNEB0) 
        {
            std::cerr << "ERROR(NLIST): Increase NNEB0: " << NNEB0 << " is not enough!!!" << std::endl;
            std::cerr << I << " " << NX << " " << NY << " " << NZ << std::endl;
            std::cerr << particle[I].IX << " " << particle[I].IY << " " << particle[I].IZ << std::endl;
            return;
        }

        NPAT[NX][NY][NZ][NUM] = I;
        NPAT[NX][NY][NZ][0] = NUM;
    }
}