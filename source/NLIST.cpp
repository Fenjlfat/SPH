#include "header.h"
#include "particles.h"
#include "parametrs.h"

void NLIST(std::vector<particles> &particle, parametrs &parametr) 
{
    //parametrs param;
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

    //Вычисление размеров сетки. 
    //Размер сетки выбирается так, чтобы размер ячейки был не менее 2*HMAX 
    //чтобы гарантировать, что все соседи частицы находятся в соседних ячейках
    NMESX = std::min(static_cast<int>((XMAX - XMIN) / (2.0 * HMAX)), parametr.NMES0);
    NMESY = std::min(static_cast<int>((YMAX - YMIN) / (2.0 * HMAX)), parametr.NMES0);
    NMESZ = std::min(static_cast<int>((ZMAX - ZMIN) / (2.0 * HMAX)), parametr.NMES0);
    if (std::max({NMESX, NMESY, NMESZ}) > parametr.NMES0) 
    {
        std::cerr << "ERROR message from SPHEP_NLIST: Increase upto NMES0=" 
                  << std::max({NMESX, NMESY, NMESZ}) << std::endl;
        return;
    }

    //Вычисляются коэффициенты для преобразования координат в индексы ячеек
    double DXX = 1.0 / (XMAX - XMIN);
    double DYY = 1.0 / (YMAX - YMIN);
    double DZZ = 1.0 / (ZMAX - ZMIN);

    // Заполнение сетки частицами
    for (int I = 0; I < particle.size(); ++I) 
    {
        //Для каждой частицы вычисляются индексы ячейки, в которую она попадает
        int NX = static_cast<int>((particle[I].IX - XMIN) * DXX * (NMESX - 1));
        int NY = static_cast<int>((particle[I].IY - YMIN) * DYY * (NMESY - 1));
        int NZ = static_cast<int>((particle[I].IZ - ZMIN) * DZZ * (NMESZ - 1));
        //NX = std::min(NX, NMESX);
        //NY = std::min(NY, NMESY);
        //NZ = std::min(NZ, NMESZ);

        //Индексы ячейки сохраняются в структуре частицы
        particle[I].MESH[0] = NX; //MESH[I][IX] = NX;
        particle[I].MESH[1] = NY; //MESH[I][IY] = NY;
        particle[I].MESH[2] = NZ; //MESH[I][IZ] = NZ;

        int NUM = parametr.NPAT[NX][NY][NZ][0] + 1;
        if (NUM > parametr.NNEB0) 
        {
            std::cerr << "ERROR(NLIST): Increase NNEB0: " << parametr.NNEB0 << " is not enough!!!" << std::endl;
            std::cerr << I << " " << NX << " " << NY << " " << NZ << std::endl;
            std::cerr << particle[I].IX << " " << particle[I].IY << " " << particle[I].IZ << std::endl;
            return;
        }
        //Частица добавляется в соответствующую ячейку сетки
        parametr.NPAT[NX][NY][NZ][NUM] = I;
        //Записывается общее количество частиц находящихся в ячейке 
        parametr.NPAT[NX][NY][NZ][0] = NUM;
    }
}