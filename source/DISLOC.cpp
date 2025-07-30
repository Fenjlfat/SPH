// Константы и глобальные параметры
constexpr double AMOL = 64.685;       // кг/моль (для Cu)
constexpr double EpsL = 8.0 * 1.6e-19; // эВ/бор
constexpr double VIMB = 0.0642;       // скорость иммобилизации
constexpr double RHOD0 = 1.0e11;      // 1/m²
constexpr double DKA = 3.392;         // коэффициент аннигиляции


void SPHEP_DISLOC(std::vector<MaterialFields>& FS, std::vector<std::vector<DislocationData>>& DISL, int NPT, int NSS0, double DTAU, double DNA) 
{
    const double CC12 = 24.0 * sqrt(3.0);
    const double C_FCC = 0.5 * sqrt(2.0) * pow(4.0, 1.0/3.0);

    for (int I = 0; I < NPT; ++I) 
    {
        double CT = sqrt(FS[I].IG / FS[I].IDNS);
        double BTR = (0.0046 * FS[I].IT + 0.281) * 1.0e-5;
        
        double m1 = AMOL * 1.0e-3 / DNA;
        double DNC = FS[I].IDNS / m1;
        double d1 = pow(DNC, C1D3);
        double BURG = C_FCC / d1;
        
        double LOC_RHOD = 0.0;
        double LOC_RHOI = 0.0;
        
        // Расчет напряжения (для Cu)
        double YY = (80.0e6 + 5.075 * FS[I].IG * BURG * sqrt(FS[I].RHOI)) * 
                   (1.0 - (FS[I].IT - 300.0) / 1060.0);
        YY = (YY < 1.0e6) ? 1.0e6 : YY;
        FS[I].IYY = YY;

        for (int J = 0; J < NSS0; ++J) 
        {
            double FORCE = BURG * (FS[I].ISXX * DISL[I][J].DNBXX + 
                                   FS[I].ISYY * DISL[I][J].DNBYY + 
                                   FS[I].ISZZ * DISL[I][J].DNBZZ + 
                                  2.0 * FS[I].ISXY * DISL[I][J].DNBXY + 
                                  2.0 * FS[I].ISXZ * DISL[I][J].DNBXZ + 
                                  2.0 * FS[I].ISYZ * DISL[I][J].DNBYZ);
            
            double YKF = fabs(BURG * YY * 0.5);
            if (FORCE > YKF) 
            {
                FORCE -= YKF;
            } else if (FORCE < -YKF) 
            {
                FORCE += YKF;
            } else 
            {
                FORCE = 0.0;
            }

            double ff = FORCE / (CT * BTR);
            double aff = fabs(ff);
            double DMN = pow(108.0 * aff + CC12 * sqrt(1.0 + 6.75 * pow(aff, 2)), C1D3);
            
            DISL[I][J].JVD = 0.0;
            if (aff > 0.0) 
            {
                DISL[I][J].JVD = CT * ff * pow(sqrt(fabs(DMN / (6.0 * aff) - 2.0/(aff * DMN))), 3);
            }

            // Обновление компонент мегатензора при вращении
            DISL[I][J].DNBXX += (-2.0 * DISL[I][J].DNBXY * FS[I].IRXY - 
                                 2.0 * DISL[I][J].DNBXZ * FS[I].IRXZ) * DTAU;
            DISL[I][J].DNBYY += (2.0 * DISL[I][J].DNBXY * FS[I].IRXY - 
                                 2.0 * DISL[I][J].DNBYZ * FS[I].IRYZ) * DTAU;
            DISL[I][J].DNBZZ += (2.0 * DISL[I][J].DNBXZ * FS[I].IRXZ + 
                                 2.0 * DISL[I][J].DNBYZ * FS[I].IRYZ) * DTAU;
            DISL[I][J].DNBXY += ((DISL[I][J].DNBXX - DISL[I][J].DNBYY) * FS[I].IRXY - 
                                 DISL[I][J].DNBXZ * FS[I].IRYZ - 
                                 DISL[I][J].DNBYZ * FS[I].IRXZ) * DTAU;
            DISL[I][J].DNBXZ += ((DISL[I][J].DNBXX - DISL[I][J].DNBZZ) * FS[I].IRXZ + 
                                DISL[I][J].DNBXY * FS[I].IRYZ - 
                                DISL[I][J].DNBYZ * FS[I].IRXY) * DTAU;
            DISL[I][J].DNBYZ += ((DISL[I][J].DNBYY - DISL[I][J].DNBZZ) * FS[I].IRYZ + 
                                 DISL[I][J].DNBXY * FS[I].IRXZ - 
                                 DISL[I][J].DNBXZ * FS[I].IRXY) * DTAU;
            
            // Пластическая релаксация
            double VRDT = DISL[I][J].JVD * DISL[I][J].JRHOD * BURG * DTAU;
            FS[I].IWXX += DISL[I][J].DNBXX * VRDT;
            FS[I].IWYY += DISL[I][J].DNBYY * VRDT;
            FS[I].IWZZ += DISL[I][J].DNBZZ * VRDT;
            FS[I].IWXY += DISL[I][J].DNBXY * VRDT;
            FS[I].IWXZ += DISL[I][J].DNBXZ * VRDT;
            FS[I].IWYZ += DISL[I][J].DNBYZ * VRDT;

            // Кинетические уравнения
            double DRHOD = DISL[I][J].JRHOD * fabs(FORCE * DISL[I][J].JVD) * 0.133 * BURG / EpsL;
            double DRHOI = 0.0;
            if (DISL[I][J].JRHOD > RHOD0) {
                DRHOI = VIMB * (DISL[I][J].JRHOD - RHOD0) * sqrt(DISL[I][J].JRHOI);
            }

            double DRHOA = DKA * BURG * fabs(DISL[I][J].JVD) * DISL[I][J].JRHOD * 
                          (2.0 * DISL[I][J].JRHOD + DISL[I][J].JRHOI);
            double DRHOAI = DKA * BURG * fabs(DISL[I][J].JVD) * DISL[I][J].JRHOD * DISL[I][J].JRHOI;

            DISL[I][J].JRHOD += DTAU * (DRHOD - DRHOI - DRHOA);
            DISL[I][J].JRHOI += DTAU * (DRHOI - DRHOAI);

            LOC_RHOD += DISL[I][J].JRHOD;
            LOC_RHOI += DISL[I][J].JRHOI;
        }

        FS[I].RHOD = LOC_RHOD;
        FS[I].RHOI = LOC_RHOI;
        FS[I].IETA1 = 3.0 * pow(BURG, 2) * FS[I].RHOD / (16.0 * BTR);
    }
}



#include <cmath>
#include <vector>

// Предполагаем, что DISL - это 3D массив (или вектор)
// Например: std::vector<std::vector<std::vector<double>>> DISL(NPT, std::vector<std::vector<double>>(NSS0, std::vector<double>(размер)));

for (int I = 0; I < NPT; I++) {
    for (int J = 0; J < NSS0; J++) {
        // Вычисляем B_COUNT и N_COUNT (аналогично Fortran-коду)
        int B_COUNT = (J) % 3 + 1;  // Fortran использует индексацию с 1
        int N_COUNT = (J - (B_COUNT - 1)) / 3 + 1;  // Корректировка для индексации с 0
        
        // Нормализующие константы
        double VN = 1.0 / sqrt(3.0);  // Для нормального вектора
        double VB = 1.0 / sqrt(2.0);  // Для вектора Бюргерса
        
        std::vector<double> VN1(3);  // Нормальный вектор
        std::vector<double> VB1(3);  // Вектор Бюргерса
        
        // Выбор направления нормального вектора (аналог Fortran SELECT CASE)
        switch (N_COUNT) {
            case 1:
                VN1 = {VN, VN, VN};  // Все компоненты положительные
                // Выбор направления вектора Бюргерса для case 1
                switch (B_COUNT) {
                    case 1:
                        VB1 = {-VB, 0.0, VB};
                        break;
                    case 2:
                        VB1 = {0.0, -VB, VB};
                        break;
                    case 3:
                        VB1 = {-VB, VB, 0.0};
                        break;
                }
                break;
            case 2:
                VN1 = {-VN, VN, VN};  // Первая компонента отрицательная
                switch (B_COUNT) {
                    case 1:
                        VB1 = {VB, 0.0, VB};
                        break;
                    case 2:
                        VB1 = {0.0, -VB, VB};
                        break;
                    case 3:
                        VB1 = {VB, VB, 0.0};
                        break;
                }
                break;
            case 3:
                VN1 = {VN, -VN, VN};  // Вторая компонента отрицательная
                switch (B_COUNT) {
                    case 1:
                        VB1 = {-VB, 0.0, VB};
                        break;
                    case 2:
                        VB1 = {0.0, VB, VB};
                        break;
                    case 3:
                        VB1 = {VB, VB, 0.0};
                        break;
                }
                break;
            case 4:
                VN1 = {VN, VN, -VN};  // Третья компонента отрицательная
                switch (B_COUNT) {
                    case 1:
                        VB1 = {VB, 0.0, VB};
                        break;
                    case 2:
                        VB1 = {0.0, VB, VB};
                        break;
                    case 3:
                        VB1 = {-VB, VB, 0.0};
                        break;
                }
                break;
        }
        
        // Заполнение массива DISL (предполагаем, что JRHOD, JRHOI и т.д. - это константы)
        DISL[I][J][JRHOD] = 1.0e11;                     // Плотность дислокаций
        DISL[I][J][JRHOI] = 6.648e12;                    // Исходная плотность 10^12
        DISL[I][J][JVD] = 0.0;                           // Скорость дислокаций
        DISL[I][J][DNBXX] = VN1[0] * VB1[0];             // Компоненты тензора
        DISL[I][J][DNBYY] = VN1[1] * VB1[1];
        DISL[I][J][DNBZZ] = VN1[2] * VB1[2];
        DISL[I][J][DNBXY] = 0.5 * (VN1[0] * VB1[1] + VN1[1] * VB1[0]);  // Смешанные компоненты
        DISL[I][J][DNBXZ] = 0.5 * (VN1[0] * VB1[2] + VN1[2] * VB1[0]);
        DISL[I][J][DNBYZ] = 0.5 * (VN1[1] * VB1[2] + VN1[2] * VB1[1]);
    }
}



#include "HEADER.h"
#include "SPH_DISLOC.h"


SPH_DISLOC::SPH_DISLOC(int NPT_get)
{
    
    NPT = NPT_get;
    //std::vector<std::vector<double>> DISL_m(NPT_m, std::vector<double>(DP, 0.00));
    std::vector<std::vector<std::vector<long double>>> DISL_m(NPT_get, std::vector<std::vector<long double>>(NSS0, std::vector<long double>(DP, 0.00)));
    DISL = DISL_m;
}

SPH_DISLOC::~SPH_DISLOC()
{
}


//!=======================================================================BEG_SPHEP_DISLOC============================================================================== =
double SPH_DISLOC::SPHEP_DISLOC(
    double DTAU,
    std::vector<std::vector<long double>>& FS) //subroutine SPHEP_DISLOC
{
    //use SPHEP_MEM
    //implicit none;
    VARIABLE V;
    PARAMETRS_MODELING PAR;
    
    double CT, FORCE, ff, aff, BTR, d1, DNC, BURG, YKF, YY, VRDT, DMN;
    double LOC_RHOD, LOC_RHOI, DRHOD, DRHOI, DRHOA, DRHOAI;
    double C1D3 = 1.e0 / 3.e0;
    double AMOL = 64.685e0;  //!AL !kg / mol;               //!!!!CU;
    double EpsL = 8.e0 * 1.6e-19; //!eV / b;          //!ostaetsya dly alumin takoyje
    double VIMB = 0.0642e0;                          //   !ñêîðîñòü èììîáèëèçàöèè  !posmotret v statye v jap
    double RHOD0 = 1.e11;			//!1 / m2                                  !posmotret v statye v jap
    double DKA = 3.392e0;// !10.d0			!anihilation coefficient       !posmotret v statye v jap
    double CC12 = 24.e0 * sqrt(3.e0);
    double C_FCC = 0.5e0 * sqrt(2.e0) * pow((4.e0), C1D3);
    double m1 = AMOL * 1.e-3 / PAR.DNA;
    for (int I = 0; I < NPT; I++) //do I = 1, NPT
    {
        CT = sqrt(FS[I][V.IG] / FS[I][V.IDNS]);
        BTR = (0.0046e0 * FS[I][V.IT] + 0.281e0) * 1.0e-5;        //!1.2 * 10 ^ 5   !CU
        //!!!!!!!!!!!!!!!!BTR = (0.0046d0 * FS(I, IT) + 0.281d0) * 1.0d - 5; //        !1.2 * 10 ^ 5   !AL    !posmotret v statye s KVS
        //!Burgers vector;        
        DNC = FS[I][V.IDNS] / m1;
        d1 = pow(DNC,C1D3);
        BURG = C_FCC / d1;
        LOC_RHOD = 0.e0;
        LOC_RHOI = 0.e0;
        YY = (80.e6 + 5.075e0 * FS[I][V.IG] * BURG * sqrt(FS[I][V.RHOI])) * (1.e0 - (FS[I][V.IT] - 300.e0) / 1060.e0);   //!CU;
        //!YY = (20.d6 + 5.075d0 !êîýôôèöèåíò óïðî÷íåíèÿ !posmotret v statye v jap;
        //!*FS(I, IG)* BURG* dsqrt(FS(I, RHOI)))* (1.d0 - (FS(I, IT) - 300.d0) / 633.d0)   !AL;
        if (YY < 1.e6) YY = 1.e6;
        FS[I][V.IYY] = YY;
        for (int J = 0; J < NSS0; J++) //do J = 1, NSS0
        {
            FORCE = BURG * (FS[I][V.ISXX] * DISL[I][J][DNBXX] + FS[I][V.ISYY] * DISL[I][J][DNBYY] + FS[I][V.ISZZ] * DISL[I][J][DNBZZ] +
                FS[I][V.ISXY] * 2.e0 * DISL[I][J][DNBXY] +
                FS[I][V.ISXZ] * 2.e0 * DISL[I][J][DNBXZ] +
                FS[I][V.ISYZ] * 2.e0 * DISL[I][J][DNBYZ]);
            YKF = abs(BURG * YY * 0.5e0);
            if (FORCE > YKF) //then
            {
                FORCE = FORCE - YKF;
            }
            else if (FORCE < -YKF) //then
            {
                FORCE = FORCE + YKF;
            }
            else
            {
                FORCE = 0.e0;
            }//endif

            ff = FORCE / (CT * BTR);
            aff = abs(ff);
            DMN = pow((108.e0 * aff + CC12 * sqrt(1.e0 + 6.75e0 * pow(aff,2))),C1D3);     //!AL and CU
            if (aff > 0.e0)
            {
                DISL[I][J][JVD] = CT * ff * pow((sqrt(abs(DMN / (6.e0 * aff) - 2.e0 / (aff * DMN)))), 3);
            }
            //!change of mega_ij components at rotation;
            DISL[I][J][DNBXX] = DISL[I][J][DNBXX] + (-2.e0 * DISL[I][J][DNBXY] * FS[I][V.IRXY] - 2.e0 * DISL[I][J][DNBXZ] * FS[I][V.IRXZ]) * DTAU;
            
            DISL[I][J][DNBYY] = DISL[I][J][DNBYY] + (2.e0 * DISL[I][J][DNBXY] * FS[I][V.IRXY] - 2.e0 * DISL[I][J][DNBYZ] * FS[I][V.IRYZ]) * DTAU;
            
            DISL[I][J][DNBZZ] = DISL[I][J][DNBZZ] + (2.e0 * DISL[I][J][DNBXZ] * FS[I][V.IRXZ] + 2.e0 * DISL[I][J][DNBYZ] * FS[I][V.IRYZ]) * DTAU;
            
            DISL[I][J][DNBXY] = DISL[I][J][DNBXY] +
                ((DISL[I][J][DNBXX] - DISL[I][J][DNBYY]) * FS[I][V.IRXY] - DISL[I][J][DNBXZ] *
                    FS[I][V.IRYZ] - DISL[I][J][DNBYZ] * FS[I][V.IRXZ]) * DTAU;
            
            DISL[I][J][DNBXZ] = DISL[I][J][DNBXZ] +
                ((DISL[I][J][DNBXX] - DISL[I][J][DNBZZ]) * FS[I][V.IRXZ] + DISL[I][J][DNBXY] *
                    FS[I][V.IRYZ] - DISL[I][J][DNBYZ] * FS[I][V.IRXY]) * DTAU;
            
            DISL[I][J][DNBYZ] = DISL[I][J][DNBYZ] +
                ((DISL[I][J][DNBYY] - DISL[I][J][DNBZZ]) * FS[I][V.IRYZ] + DISL[I][J][DNBXY] *
                    FS[I][V.IRXZ] - DISL[I][J][DNBXZ] * FS[I][V.IRXY]) * DTAU;
            
            if (DISL[I][J][DNBXY]!= DISL[I][J][DNBXY])
            {
                std::cout << "ERROR DNBXY";
                getchar();
            }
            if (DISL[I][J][DNBXZ] != DISL[I][J][DNBXZ])
            {
                std::cout << "ERROR DNBXZ";
                getchar();
            }
            if (DISL[I][J][DNBYZ] != DISL[I][J][DNBYZ])
            {
                std::cout << "ERROR DNBYZ";
                getchar();
            }
            if (DISL[I][J][DNBXX] != DISL[I][J][DNBXX])
            {
                std::cout << "ERROR DNBXX";
                getchar();
            }
            if (DISL[I][J][DNBYY] != DISL[I][J][DNBYY])
            {
                std::cout << "ERROR DNBYY";
                getchar();
            }
            if (DISL[I][J][DNBZZ] != DISL[I][J][DNBZZ])
            {
                std::cout << "ERROR DNBZZ";
                getchar();
            }
            

            //!plastic relaxation;
            VRDT = DISL[I][J][JVD] * DISL[I][J][JRHOD] * BURG * DTAU;
            if (VRDT != VRDT)
            {
                std::cout << "ERROR VRDT";
                getchar();
            }
            FS[I][V.IWXX] = FS[I][V.IWXX] + DISL[I][J][DNBXX] * VRDT;
            FS[I][V.IWYY] = FS[I][V.IWYY] + DISL[I][J][DNBYY] * VRDT;
            FS[I][V.IWZZ] = FS[I][V.IWZZ] + DISL[I][J][DNBZZ] * VRDT;
            FS[I][V.IWXY] = FS[I][V.IWXY] + DISL[I][J][DNBXY] * VRDT;
            FS[I][V.IWXZ] = FS[I][V.IWXZ] + DISL[I][J][DNBXZ] * VRDT;
            FS[I][V.IWYZ] = FS[I][V.IWYZ] + DISL[I][J][DNBYZ] * VRDT;
            /*std::cout << "I=" << I << "  J=" << J << "   " << " DISL[I][J][DNBXX]=" << DISL[I][J][DNBXX] << " DISL[I][J][DNBYY]=" << DISL[I][J][DNBYY]
                << " DISL[I][J][DNBXY]=" << DISL[I][J][DNBXY];
            getchar();*/
            //!kinetic equations;
            DRHOD = DISL[I][J][JRHOD] * abs(FORCE * DISL[I][J][JVD]) * 0.133e0 * BURG / EpsL;           //!? ? ? AL and CU;
            if (DRHOD != DRHOD)
            {
                std::cout << "ERROR DRHOD";
                getchar();
            }
            DRHOI = 0.e0;
            //double KOEF_JROI = abs(DISL[I][J][JRHOI]);
            if (DISL[I][J][JRHOD] > RHOD0) DRHOI = VIMB * (DISL[I][J][JRHOD] - RHOD0) * sqrt(DISL[I][J][JRHOI]);
            if (DRHOI != DRHOI)
            {
                std::cout << "ERROR DRHOI";
                getchar();
            }
            
            //if (DISL[I][J][JRHOD] > RHOD0) DRHOI = VIMB * (DISL[I][J][JRHOD] - RHOD0) * sqrt(KOEF_JROI);
            //!;
            DRHOA = DKA * BURG * abs(DISL[I][J][JVD]) * DISL[I][J][JRHOD] * (2.e0 * DISL[I][J][JRHOD] + DISL[I][J][JRHOI]);
            DRHOAI = DKA * BURG * abs(DISL[I][J][JVD]) * DISL[I][J][JRHOD] * DISL[I][J][JRHOI];;
            //!;
            DISL[I][J][JRHOD] = DISL[I][J][JRHOD] + DTAU * (DRHOD - DRHOI - DRHOA);
            /*
            std::cout
                
                << "   DRHOD=" << DRHOD << "   DISL[I][J][JRHOD]=" << DISL[I][J][JRHOD] << "   DISL[I][J][JVD]=" << DISL[I][J][JVD] << "\n"
                << "   DRHOI=" << DRHOI << "   DISL[I][J][JRHOI]=" << DISL[I][J][JRHOI] << "\n"
                << "   DRHOA=" << DRHOA << "\n"
                << "   DRHOAI=" << DRHOAI << "\n"
                << "   aff=" << aff << "   FORCE=" << FORCE << "   CT=" << CT << "   BTR=" << BTR
                << "   IG=" << FS[I][V.IG]
                << "   IDNS=" << FS[I][V.IDNS] << "\n"
                << std::endl;
                */

            DISL[I][J][JRHOI] = DISL[I][J][JRHOI] + DTAU * (DRHOI - DRHOAI);
            if (DISL[I][J][JRHOI] <= 0.00)
            {

                std::cout << "ERROR DRHOI";
                std::cout << "DISL[I][J][JRHOI]=" << DISL[I][J][JRHOI];
                std::cout << "  DRHOI=" << DRHOI;
                std::cout << "  DRHOAI=" << DRHOAI;
                getchar();
            }
            /*if (DISL[I][J][JRHOI] < 0.00)
            {
                DISL[I][J][JRHOI] = 0.00;
            }*/
            LOC_RHOD = LOC_RHOD + DISL[I][J][JRHOD];
            LOC_RHOI = LOC_RHOI + DISL[I][J][JRHOI];
        }
        FS[I][V.RHOD] = LOC_RHOD;
        FS[I][V.RHOI] = LOC_RHOI;
        FS[I][V.IETA1] = 3.e0 * pow(BURG,2) * FS[I][V.RHOD] / (16.e0 * BTR);
    }//enddo
    //6 format(1pe16.8, 20(e16.8));
    return 0;
}//end subroutine SPHEP_DISLOC
//!=======================================================================END_SPHEP_DISLOC============================================================================== =



void SPH_DISLOC::SPHEP_DISL()
{
    for (int I = 0; I < NPT; I++) {
        for (int J = 1; J <= NSS0; J++) {
            int B_COUNT = (J - 1) % 3 + 1;
            int N_COUNT = static_cast<int>((J - B_COUNT) / 3 + 1);
            //std::cout << "B_COUNT=" << B_COUNT << std::endl;
            //std::cout << "N_COUNT=" << N_COUNT << std::endl;
            //getchar();
            double VN = 1.0 / std::sqrt(3.0);
            double VB = 1.0 / std::sqrt(2.0);
            std::vector<double> VN1(3);
            std::vector<double> VB1(3);

            switch (N_COUNT) {
            case 1:
                VN1 = { VN, VN, VN };
                switch (B_COUNT) {
                case 1: VB1 = { -VB, 0.0, VB }; break;
                case 2: VB1 = { 0.0, -VB, VB }; break;
                case 3: VB1 = { -VB, VB, 0.0 }; break;
                }
                break;
            case 2:
                VN1 = { -VN, VN, VN };
                switch (B_COUNT) {
                case 1: VB1 = { VB, 0.0, VB }; break;
                case 2: VB1 = { 0.0, -VB, VB }; break;
                case 3: VB1 = { VB, VB, 0.0 }; break;
                }
                break;
            case 3:
                VN1 = { VN, -VN, VN };
                switch (B_COUNT) {
                case 1: VB1 = { -VB, 0.0, VB }; break;
                case 2: VB1 = { 0.0, VB, VB }; break;
                case 3: VB1 = { VB, VB, 0.0 }; break;
                }
                break;
            case 4:
                VN1 = { VN, VN, -VN };
                switch (B_COUNT) {
                case 1: VB1 = { VB, 0.0, VB }; break;
                case 2: VB1 = { 0.0, VB, VB }; break;
                case 3: VB1 = { -VB, VB, 0.0 }; break;
                }
                break;
            }

            /*std::cout << "I=" << I << "  J=" << J - 1 << "   ";
            std::cout << "VN1{";
            for (auto var : VN1)
            {
                std::cout << var << "; ";
            }

            std::cout << "       VB1{";
            for (auto var : VB1)
            {
                std::cout << var << "; ";
            }
            getchar();*/
            DISL[I][J - 1][JRHOD] = 1.0e11;
            DISL[I][J - 1][JRHOI] = 6.648e12; // 6.3e12 ORIGIN 10^12
            DISL[I][J - 1][JVD] = 0.0;
            DISL[I][J - 1][DNBXX] = VN1[0] * VB1[0];
            DISL[I][J - 1][DNBYY] = VN1[1] * VB1[1];
            DISL[I][J - 1][DNBZZ] = VN1[2] * VB1[2];
            DISL[I][J - 1][DNBXY] = 0.5 * (VN1[0] * VB1[1] + VN1[1] * VB1[0]);
            DISL[I][J - 1][DNBXZ] = 0.5 * (VN1[0] * VB1[2] + VN1[2] * VB1[0]);
            DISL[I][J - 1][DNBYZ] = 0.5 * (VN1[1] * VB1[2] + VN1[2] * VB1[1]);
            /*std::cout << "I=" << I << "  J=" << J - 1 << "   " << " DISL[I][J][DNBXX]=" << DISL[I][J - 1][DNBXX] << " DISL[I][J][DNBYY]=" << DISL[I][J - 1][DNBYY]
                << " DISL[I][J][DNBXY]=" << DISL[I][J - 1][DNBXY];
            getchar();*/
        }
    }

}