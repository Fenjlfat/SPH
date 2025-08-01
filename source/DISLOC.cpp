#include "header.h"
#include "parametrs.h"
#include "particles.h"


void SPHEP_DISLOC(std::vector<particles>& particle, double DTAU) 
{
    parametrs param;

    //for (int I = 0; I < NPT; ++I) 
    for(auto &p: particle)
    {
        double CT = sqrt(p.IG / p.IDNS);
        double BTR = (0.0046 * p.IT + 0.281) * 1.0e-5;
        
        double m1 = param.AMOL * 1.0e-3 / param.DNA;
        double DNC = p.IDNS / m1;
        double d1 = pow(DNC, param.C1D3);
        double BURG = param.C_FCC / d1;
        
        double LOC_RHOD = 0.0;
        double LOC_RHOI = 0.0;
        
        // Расчет напряжения (для Cu)
        double YY = (80.0e6 + 5.075 * p.IG * BURG * sqrt(p.RHOI)) * (1.0 - (p.IT - 300.0) / 1060.0);
        YY = (YY < 1.0e6) ? 1.0e6 : YY;
        p.IYY = YY;

        for (int J = 0; J < param.NSS0; ++J) 
        {
            double FORCE = BURG * (p.ISXX * p.dislocation[J][param.DNBXX] + 
                                   p.ISYY * p.dislocation[J][param.DNBYY] + 
                                   p.ISZZ * p.dislocation[J][param.DNBZZ] + 
                                  2.0 * p.ISXY * p.dislocation[J][param.DNBXY] + 
                                  2.0 * p.ISXZ * p.dislocation[J][param.DNBXZ] + 
                                  2.0 * p.ISYZ * p.dislocation[J][param.DNBYZ]);
            
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
            double DMN = pow(108.0 * aff + param.CC12 * sqrt(1.0 + 6.75 * pow(aff, 2)), param.C1D3);
            
            p.dislocation[J][param.JVD] = 0.0;
            if (aff > 0.0) 
            {
                p.dislocation[J][param.JVD] = CT * ff * pow(sqrt(fabs(DMN / (6.0 * aff) - 2.0/(aff * DMN))), 3);
            }

            // Обновление компонент тензора omega при вращении
            p.dislocation[J][param.DNBXX] += (-2.0 * p.dislocation[J][param.DNBXY] * p.IRXY - 2.0 * p.dislocation[J][param.DNBXZ] * p.IRXZ) * DTAU;
            p.dislocation[J][param.DNBYY] +=  (2.0 * p.dislocation[J][param.DNBXY] * p.IRXY - 2.0 * p.dislocation[J][param.DNBYZ] * p.IRYZ) * DTAU;
            p.dislocation[J][param.DNBZZ] +=  (2.0 * p.dislocation[J][param.DNBXZ] * p.IRXZ + 2.0 * p.dislocation[J][param.DNBYZ] * p.IRYZ) * DTAU;
            p.dislocation[J][param.DNBXY] += ((p.dislocation[J][param.DNBXX] - p.dislocation[J][param.DNBYY]) * p.IRXY - 
                                               p.dislocation[J][param.DNBXZ] * p.IRYZ - p.dislocation[J][param.DNBYZ] * p.IRXZ) * DTAU;
            p.dislocation[J][param.DNBXZ] += ((p.dislocation[J][param.DNBXX] - p.dislocation[J][param.DNBZZ]) * p.IRXZ + 
                                               p.dislocation[J][param.DNBXY] * p.IRYZ - p.dislocation[J][param.DNBYZ] * p.IRXY) * DTAU;
            p.dislocation[J][param.DNBYZ] += ((p.dislocation[J][param.DNBYY] - p.dislocation[J][param.DNBZZ]) * p.IRYZ + 
                                               p.dislocation[J][param.DNBXY] * p.IRXZ - p.dislocation[J][param.DNBXZ] * p.IRXY) * DTAU;
            
            // Пластическая релаксация
            double VRDT = p.dislocation[J][param.JVD] * p.dislocation[J][param.JRHOD] * BURG * DTAU;
            p.IWXX += p.dislocation[J][param.DNBXX] * VRDT;
            p.IWYY += p.dislocation[J][param.DNBYY] * VRDT;
            p.IWZZ += p.dislocation[J][param.DNBZZ] * VRDT;
            p.IWXY += p.dislocation[J][param.DNBXY] * VRDT;
            p.IWXZ += p.dislocation[J][param.DNBXZ] * VRDT;
            p.IWYZ += p.dislocation[J][param.DNBYZ] * VRDT;

            // Кинетические уравнения
            double DRHOD = p.dislocation[J][param.JRHOD] * fabs(FORCE * p.dislocation[J][param.JVD]) * 0.133 * BURG / param.EpsL;
            double DRHOI = 0.0;
            if (p.dislocation[J][param.JRHOD] > param.RHOD0) {
                DRHOI = param.VIMB * (p.dislocation[J][param.JRHOD] - param.RHOD0) * sqrt(p.dislocation[J][param.JRHOI]);
            }

            double DRHOA = param.DKA * BURG * fabs(p.dislocation[J][param.JVD]) * p.dislocation[J][param.JRHOD] * 
                          (2.0 * p.dislocation[J][param.JRHOD] + p.dislocation[J][param.JRHOI]);
            double DRHOAI = param.DKA * BURG * fabs(p.dislocation[J][param.JVD]) * p.dislocation[J][param.JRHOD] * p.dislocation[J][param.JRHOI];

            p.dislocation[J][param.JRHOD] += DTAU * (DRHOD - DRHOI - DRHOA);
            p.dislocation[J][param.JRHOI] += DTAU * (DRHOI - DRHOAI);

            LOC_RHOD += p.dislocation[J][param.JRHOD];
            LOC_RHOI += p.dislocation[J][param.JRHOI];
        }

        p.RHOD = LOC_RHOD;
        p.RHOI = LOC_RHOI;
        p.IETA1 = 3.0 * pow(BURG, 2) * p.RHOD / (16.0 * BTR);
    }
}