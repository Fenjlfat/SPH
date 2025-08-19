#include "header.h"
#include "parametrs.h"
#include "particles.h"

void SLRelax(double &RO, double &U, double &T,  double &P, double &CS, double &CV, double SUBS);

//initialization vector of structure particle
void INIT(std::vector<particles> &particle, parametrs &parametr)
{
    //parametrs param;
    
    int idParticle = 0; //number particle
    //double *xx;
    //double *yy;
    //double *zz;
    
    double E_EOS = 0.0, P_EOS = 0.0, T_EOS = 300.0, CV_EOS = 0.0, CS_EOS = 0.0, DNS_EOS = 8300.0;
    //EOS(E_EOS, P_EOS, T_EOS, CV_EOS, CS_EOS, DNS_EOS);
    SLRelax(DNS_EOS, E_EOS, T_EOS, P_EOS, CS_EOS, CV_EOS, parametr.SUBZ); 

    for (int i = 0; i < parametr.npaax; i++)
    {
        for (int j = 0; j < parametr.npaay; j++)
        {
            for (int k = 0; k < parametr.npaaz; k++)
            {
                double xx = i * parametr.dbp;
                double yy = j * parametr.dbp;
                double zz = k * parametr.dbp;
                //koordinate
                particle[idParticle].IX = xx;
                particle[idParticle].IY = yy;
                particle[idParticle].IZ = zz;
                //velocity
                particle[idParticle].IVX = 0.00;
                particle[idParticle].IVY = 0.00;
                particle[idParticle].IVZ = 0.00;
                
                //if (xx > (parametr.npaax - 2) * parametr.dbp)
                if(xx < 0.001)
                {
                    particle[idParticle].IVX = 100.00;
                }
                
                

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
				particle[idParticle].IMAS = particle[idParticle].IDNS * pow(parametr.dbp, 3);
				particle[idParticle].IKK = particle[idParticle].IDNS * pow(particle[idParticle].ICS, 2);
				particle[idParticle].ICP = parametr.nuP;
				particle[idParticle].IG = 1.5e0 * particle[idParticle].IKK * 
                    (1.e0 - 2.e0 * particle[idParticle].ICP) / (1.e0 + particle[idParticle].ICP);
				particle[idParticle].IHS = 1.e0 * parametr.dbp;
				//std::cout << " IHS=" << FS[id_Particle][VAR.IHS] << std::endl;
				//getchar();

                idParticle += 1;
            }
            
        }
        
    }

    // initialization vector disl
    // std::vector<std::vector<std::vector<double>>> DISL(NPT, std::vector<std::vector<double>>(NSS0, std::vector<double>(размер)));
    for (int I = 0; I < idParticle; I++) {
        for (int J = 0; J < parametr.NSS0; J++) {
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

            particle[I].dislocation[J][parametr.JRHOD] = 1.0e11;          // Плотность дислокаций
            particle[I].dislocation[J][parametr.JRHOI] = 6.648e12;        // Исходная плотность 10^12
            particle[I].dislocation[J][parametr.JVD] = 0.0;               // Скорость дислокаций
            particle[I].dislocation[J][parametr.DNBXX] = VN1[0] * VB1[0];  // Компоненты тензора
            particle[I].dislocation[J][parametr.DNBYY] = VN1[1] * VB1[1];
            particle[I].dislocation[J][parametr.DNBZZ] = VN1[2] * VB1[2];
            particle[I].dislocation[J][parametr.DNBXY] = 0.5 * (VN1[0] * VB1[1] + VN1[1] * VB1[0]);   // Смешанные компоненты
            particle[I].dislocation[J][parametr.DNBXZ] = 0.5 * (VN1[0] * VB1[2] + VN1[2] * VB1[0]);
            particle[I].dislocation[J][parametr.DNBYZ] = 0.5 * (VN1[1] * VB1[2] + VN1[2] * VB1[1]);
        }
    }
}
