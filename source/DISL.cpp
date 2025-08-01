#include "header.h"
void DISL()
{
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

}

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