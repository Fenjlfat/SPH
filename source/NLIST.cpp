#include "HEADER.h"   //use SPHEP_MEM
#include "SPH_NLIST.h"

SPH_NLIST::SPH_NLIST(int NPT_get)
{
    NPT = NPT_get;
    std::vector<std::vector<std::vector<std::vector<int>>>> 
        NPAT_m(NMES0, std::vector < std::vector < std::vector<int>>>(NMES0, std::vector<std::vector<int>>(NMES0, std::vector<int>(NNEB0, 0))));
    NPAT = NPAT_m;
    std::vector<std::vector<int>> MESH_m(NPT_get, std::vector<int>(3, 0));
    MESH = MESH_m;
}

SPH_NLIST::~SPH_NLIST()
{
}

void ZERO(int NMES0, int NNEB0, std::vector<std::vector<std::vector<std::vector<int>>>> &NPAT)
{
    for (int i = 0; i < NMES0; i++)
    {
        for (int j = 0; j < NMES0; j++)
        {
            for (int k = 0; k < NMES0; k++)
            {
                for (int n = 0; n < NNEB0; n++)
                {
                    NPAT[i][j][k][n] = 0;
                }
            }
            
        }
    }
}


int MAX(int a, int b, int c)
{
    if (a > b && a > c)
    {
        return a;
    }
    else if (b > a && b > c)
    {
        return b;
    }
    else
    {
        return c;
    }
}


void SPH_NLIST::SPHEP_NLIST(std::vector<std::vector<long double>> &FS)
{
    // Assuming SPHEP_MEM is a header file that defines necessary variables and functions
    // and FS is a function that returns a double value based on the input indices.

    VARIABLE V;
    PARAMETRS_MODELING PAR;

    int NX, NY, NZ, NUM{0};
    double DXX, DYY, DZZ;
    
    double HMAX = FS[0][V.IHS];
    double XX, YY, ZZ;
    double XMAX = FS[0][V.IX];
    double XMIN = XMAX;
    double YMAX = FS[0][V.IY];
    double YMIN = YMAX;
    double ZMAX = FS[0][V.IZ];
    double ZMIN = ZMAX;

    for (int I = 1; I < NPT; I++) //do I = 2, NPT
    {     
        XX = FS[I][V.IX];
        YY = FS[I][V.IY];
        ZZ = FS[I][V.IZ];

        if (XMIN > XX) XMIN = XX;
        if (XMAX < XX) XMAX = XX;
        if (YMIN > YY) YMIN = YY;
        if (YMAX < YY) YMAX = YY;
        if (ZMIN > ZZ) ZMIN = ZZ;
        if (ZMAX < ZZ) ZMAX = ZZ;

        if (HMAX < FS[I][V.IHS]) HMAX = FS[I][V.IHS];
    }

    if (HMAX <= 0.0) 
    {
        std::cout << "ERROR message from SPHEP_NLIST: ZERO Kernel Scale HMAX" << std::endl;
        NLIST_EXIT = 1;
    }

    NMESX = std::min(static_cast<int>((XMAX - XMIN) / (2.0 * HMAX)), NMES0);
    NMESY = std::min(static_cast<int>((YMAX - YMIN) / (2.0 * HMAX)), NMES0);
    NMESZ = std::min(static_cast<int>((ZMAX - ZMIN) / (2.0 * HMAX)), NMES0);
    //std::cout << "NMESX=" << NMESX << "  NMESY=" << NMESY << "  NMESZ=" << NMESZ << std::endl;
    //if (NMESX > NMES0 || NMESY > NMES0 || NMESZ > NMES0)
    if (MAX( NMESX, NMESY, NMESZ ) > NMES0) 
    {
        std::cout << "ERROR message from SPHEP_NLIST: Increase upto NMES0=" << MAX( NMESX, NMESY, NMESZ ) << std::endl;
        NLIST_EXIT = 1;
    }

    DXX = 1.0 / (XMAX - XMIN);
    DYY = 1.0 / (YMAX - YMIN);
    DZZ = 1.0 / (ZMAX - ZMIN);
    //std::cout << "DXX=" << DXX << "  DYY=" << DYY << "  DZZ=" << DZZ << std::endl;
    ZERO(NMES0, NNEB0, NPAT); // ôóíêöèÿ îáíóëåíèÿ âåêòîðà
    for (int I = 1; I <= NPT; I++) {
        NX = static_cast<int>((FS[I-1][V.IX] - XMIN) * DXX * (NMESX - 1)) + 1;
        NY = static_cast<int>((FS[I-1][V.IY] - YMIN) * DYY * (NMESY - 1)) + 1;
        NZ = static_cast<int>((FS[I-1][V.IZ] - ZMIN) * DZZ * (NMESZ - 1)) + 1;
        //std::cout << "NUM=" << NUM << "  NX=" << NX << "  NY=" << NY << "  NZ=" << NZ << std::endl;
        //std::cout << "NUM=" << "MESH=" << MESH[I][V.IX] << "  NPAT=" << NPAT[NX][NY][NZ][NUM] << std::endl;
        if (NX > NMESX) NX = NMESX;
        if (NY > NMESY) NY = NMESY;
        if (NZ > NMESZ) NZ = NMESZ;

        MESH[I-1][V.IX] = NX;
        MESH[I-1][V.IY] = NY;
        MESH[I-1][V.IZ] = NZ;

        NUM = NPAT[NX][NY][NZ][0] + 1;
        if (NUM > NNEB0) {
            std::cout << "ERROR(NLIST): Increase NNEB0: " << NNEB0 << " is not enough!!!" << std::endl;
            std::cout << I-1 << " " << NX << " " << NY << " " << NZ << std::endl;
            std::cout << FS[I-1][V.IX] << " " << FS[I-1][V.IY] << " " << FS[I-1][V.IZ] << std::endl;
            NLIST_EXIT = 1;
            break;
        }
        NPAT[NX][NY][NZ][NUM] = I+1;
        
        NPAT[NX][NY][NZ][0] = NUM;
        /*std::cout << "NPT=" << NPT << std::endl;
        std::cout << "NUM=" << NUM << "  NX=" << NX << "  NY=" << NY << "  NZ=" << NZ <<  std::endl;
        std::cout << "NUM=" << "MESH=" << MESH[I][V.IX] << "  NPAT=" << NPAT[NX][NY][NZ][NUM] << std::endl;
        getchar();*/

    }
    /*
    for (int i = 1; i <= 50; i++)
    {
        for (int j = 1; j <= 50; j++)
        {
            for (int k = 1; k <= 50; k++)
            {
                for (int n = 1; n <= 500; n++)
                {
                    std::cout << "NPAT[" << i << "][" << j << "][" << k << "][" << n << "] = " << NPAT[NX][NY][NZ][NUM] << std::endl;
                    getchar();
                }
            }
        }

    }
    getchar();
    */
    //std::cout << "hello world!" << std::endl;
    
}