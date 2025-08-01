#include "header.h"
#include "parametrs.h"
#include "particles.h"

//initialization vector of structure particle
void INIT(double &E_EOS, double &P_EOS, double &T_EOS, double &CV_EOS, double &CS_EOS, double &DNS_EOS, std::vector<particles> &particle)
{
    parametrs param;
    
    int idParticle = 0; //number particle
    //double *xx;
    //double *yy;
    //double *zz;
    

    for (int i = 0; i < param.npaax; i++)
    {
        for (int j = 0; j < param.npaay; j++)
        {
            for (int k = 0; k < param.npaaz; k++)
            {
                double xx = i * param.dbp;
                double yy = j * param.dbp;
                double zz = k * param.dbp;
                //koordinate
                particle[idParticle].IX = xx;
                particle[idParticle].IY = yy;
                particle[idParticle].IZ = zz;
                //velocity
                particle[idParticle].IVX = 0.00;
                particle[idParticle].IVY = 0.00;
                particle[idParticle].IVZ = 0.00;
                
                /*for (int ALF = 0; ALF < 3; ALF++)
				{
					*particle[idParticle].IVV_Ptr[ALF] = 0.e0; //- 98.d0;
					*particle[idParticle].IVREALV_Ptr[ALF] = *particle[idParticle].IVV_Ptr[ALF];
				}*/

                *particle[idParticle].IVV_Ptr[0] = particle[idParticle].IVX;
                *particle[idParticle].IVV_Ptr[1] = particle[idParticle].IVY;
                *particle[idParticle].IVV_Ptr[2] = particle[idParticle].IVZ; 
                
                *particle[idParticle].IVREALV_Ptr[0] = *particle[idParticle].IVV_Ptr[0];
                *particle[idParticle].IVREALV_Ptr[1] = *particle[idParticle].IVV_Ptr[1];
                *particle[idParticle].IVREALV_Ptr[2] = *particle[idParticle].IVV_Ptr[2];

                //thermodynamic parameters
				particle[idParticle].IU = E_EOS;
                particle[idParticle].IP = P_EOS;
				particle[idParticle].IT = T_EOS;
                particle[idParticle].ICV = CV_EOS;
                particle[idParticle].ICS = CS_EOS;
                particle[idParticle].IDNS = DNS_EOS;
				particle[idParticle].IMAS = particle[idParticle].IDNS * pow(param.dbp, 3);
				particle[idParticle].IKK = particle[idParticle].IDNS * pow(particle[idParticle].ICS, 2);
				particle[idParticle].ICP = param.nuP;
				particle[idParticle].IG = 1.5e0 * particle[idParticle].IKK * 
                    (1.e0 - 2.e0 * particle[idParticle].ICP) / (1.e0 + particle[idParticle].ICP);
				particle[idParticle].IHS = 1.e0 * param.dbp;
				//std::cout << " IHS=" << FS[id_Particle][VAR.IHS] << std::endl;
				//getchar();

                idParticle += 1;
            }
            
        }
        
    }

    // initialization vector disl
    // std::vector<std::vector<std::vector<double>>> DISL(NPT, std::vector<std::vector<double>>(NSS0, std::vector<double>(размер)));
    for (int I = 0; I < idParticle; I++) {
        for (int J = 0; J < param.NSS0; J++) {
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
            
            // Заполнение массива DISL 
            //DISL[I][J][JRHOD] = 1.0e11;                     // Плотность дислокаций
            //DISL[I][J][JRHOI] = 6.648e12;                    // Исходная плотность 10^12
            //DISL[I][J][JVD] = 0.0;                           // Скорость дислокаций
            //DISL[I][J][DNBXX] = VN1[0] * VB1[0];             // Компоненты тензора
            //DISL[I][J][DNBYY] = VN1[1] * VB1[1];
            //DISL[I][J][DNBZZ] = VN1[2] * VB1[2];
            //DISL[I][J][DNBXY] = 0.5 * (VN1[0] * VB1[1] + VN1[1] * VB1[0]);  // Смешанные компоненты
            //DISL[I][J][DNBXZ] = 0.5 * (VN1[0] * VB1[2] + VN1[2] * VB1[0]);
            //DISL[I][J][DNBYZ] = 0.5 * (VN1[1] * VB1[2] + VN1[2] * VB1[1]);

            particle[I].dislocation[J][param.JRHOD] = 1.0e11;          // Плотность дислокаций
            particle[I].dislocation[J][param.JRHOI] = 6.648e12;        // Исходная плотность 10^12
            particle[I].dislocation[J][param.JVD] = 0.0;               // Скорость дислокаций
            particle[I].dislocation[J][param.DNBXX]= VN1[0] * VB1[0];  // Компоненты тензора
            particle[I].dislocation[J][param.DNBYY] = VN1[1] * VB1[1];
            particle[I].dislocation[J][param.DNBZZ] = VN1[2] * VB1[2];
            particle[I].dislocation[J][param.DNBXY] = 0.5 * (VN1[0] * VB1[1] + VN1[1] * VB1[0]);   // Смешанные компоненты
            particle[I].dislocation[J][param.DNBXZ] = 0.5 * (VN1[0] * VB1[2] + VN1[2] * VB1[0]);
            particle[I].dislocation[J][param.DNBYZ] = 0.5 * (VN1[1] * VB1[2] + VN1[2] * VB1[1]);
        }
    }
}
