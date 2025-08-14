#include "header.h"
#include "parametrs.h"
#include "particles.h"


void SPHEP_DISLOC(std::vector<particles>& particle, parametrs &parametr) 
{
    //parametrs param;

    //for (int I = 0; I < NPT; ++I) 
    for(auto &p: particle)
    {
        double CT = sqrt(p.IG / p.IDNS);
        double BTR = (0.0046 * p.IT + 0.281) * 1.0e-5;
        
        double m1 = parametr.AMOL * 1.0e-3 / parametr.DNA;
        double DNC = p.IDNS / m1;
        double d1 = pow(DNC, parametr.C1D3);
        double BURG = parametr.C_FCC / d1;
        
        double LOC_RHOD = 0.0;
        double LOC_RHOI = 0.0;
        
        // Расчет напряжения (для Cu)
        double YY = (80.0e6 + 5.075 * p.IG * BURG * sqrt(p.RHOI)) * (1.0 - (p.IT - 300.0) / 1060.0);
        YY = (YY < 1.0e6) ? 1.0e6 : YY;
        p.IYY = YY;

        for (int J = 0; J < parametr.NSS0; ++J) 
        {
            double FORCE = BURG * (p.ISXX * p.dislocation[J][parametr.DNBXX] + 
                                   p.ISYY * p.dislocation[J][parametr.DNBYY] + 
                                   p.ISZZ * p.dislocation[J][parametr.DNBZZ] + 
                                  2.0 * p.ISXY * p.dislocation[J][parametr.DNBXY] + 
                                  2.0 * p.ISXZ * p.dislocation[J][parametr.DNBXZ] + 
                                  2.0 * p.ISYZ * p.dislocation[J][parametr.DNBYZ]);
            
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
            double DMN = pow(108.0 * aff + parametr.CC12 * sqrt(1.0 + 6.75 * pow(aff, 2)), parametr.C1D3);
            
            p.dislocation[J][parametr.JVD] = 0.0;
            if (aff > 0.0) 
            {
                p.dislocation[J][parametr.JVD] = CT * ff * pow(sqrt(fabs(DMN / (6.0 * aff) - 2.0/(aff * DMN))), 3);
            }

            // Обновление компонент тензора omega при вращении
            p.dislocation[J][parametr.DNBXX] += (-2.0 * p.dislocation[J][parametr.DNBXY] * p.IRXY - 2.0 * p.dislocation[J][parametr.DNBXZ] * p.IRXZ) * parametr.DTAU;
            p.dislocation[J][parametr.DNBYY] +=  (2.0 * p.dislocation[J][parametr.DNBXY] * p.IRXY - 2.0 * p.dislocation[J][parametr.DNBYZ] * p.IRYZ) * parametr.DTAU;
            p.dislocation[J][parametr.DNBZZ] +=  (2.0 * p.dislocation[J][parametr.DNBXZ] * p.IRXZ + 2.0 * p.dislocation[J][parametr.DNBYZ] * p.IRYZ) * parametr.DTAU;
            p.dislocation[J][parametr.DNBXY] += ((p.dislocation[J][parametr.DNBXX] - p.dislocation[J][parametr.DNBYY]) * p.IRXY - 
                                               p.dislocation[J][parametr.DNBXZ] * p.IRYZ - p.dislocation[J][parametr.DNBYZ] * p.IRXZ) * parametr.DTAU;
            p.dislocation[J][parametr.DNBXZ] += ((p.dislocation[J][parametr.DNBXX] - p.dislocation[J][parametr.DNBZZ]) * p.IRXZ + 
                                               p.dislocation[J][parametr.DNBXY] * p.IRYZ - p.dislocation[J][parametr.DNBYZ] * p.IRXY) * parametr.DTAU;
            p.dislocation[J][parametr.DNBYZ] += ((p.dislocation[J][parametr.DNBYY] - p.dislocation[J][parametr.DNBZZ]) * p.IRYZ + 
                                               p.dislocation[J][parametr.DNBXY] * p.IRXZ - p.dislocation[J][parametr.DNBXZ] * p.IRXY) * parametr.DTAU;
            
            // Пластическая релаксация
            double VRDT = p.dislocation[J][parametr.JVD] * p.dislocation[J][parametr.JRHOD] * BURG * parametr.DTAU;
            p.IWXX += p.dislocation[J][parametr.DNBXX] * VRDT;
            p.IWYY += p.dislocation[J][parametr.DNBYY] * VRDT;
            p.IWZZ += p.dislocation[J][parametr.DNBZZ] * VRDT;
            p.IWXY += p.dislocation[J][parametr.DNBXY] * VRDT;
            p.IWXZ += p.dislocation[J][parametr.DNBXZ] * VRDT;
            p.IWYZ += p.dislocation[J][parametr.DNBYZ] * VRDT;

            // Кинетические уравнения
            double DRHOD = p.dislocation[J][parametr.JRHOD] * fabs(FORCE * p.dislocation[J][parametr.JVD]) * 0.133 * BURG / parametr.EpsL;
            double DRHOI = 0.0;
            if (p.dislocation[J][parametr.JRHOD] > parametr.RHOD0) {
                DRHOI = parametr.VIMB * (p.dislocation[J][parametr.JRHOD] - parametr.RHOD0) * sqrt(p.dislocation[J][parametr.JRHOI]);
            }

            double DRHOA = parametr.DKA * BURG * fabs(p.dislocation[J][parametr.JVD]) * p.dislocation[J][parametr.JRHOD] * 
                          (2.0 * p.dislocation[J][parametr.JRHOD] + p.dislocation[J][parametr.JRHOI]);
            double DRHOAI = parametr.DKA * BURG * fabs(p.dislocation[J][parametr.JVD]) * p.dislocation[J][parametr.JRHOD] * p.dislocation[J][parametr.JRHOI];

            p.dislocation[J][parametr.JRHOD] += parametr.DTAU * (DRHOD - DRHOI - DRHOA);
            p.dislocation[J][parametr.JRHOI] += parametr.DTAU * (DRHOI - DRHOAI);

            LOC_RHOD += p.dislocation[J][parametr.JRHOD];
            LOC_RHOI += p.dislocation[J][parametr.JRHOI];
        }

        p.RHOD = LOC_RHOD;
        p.RHOI = LOC_RHOI;
        p.IETA1 = 3.0 * pow(BURG, 2) * p.RHOD / (16.0 * BTR);
    }
}