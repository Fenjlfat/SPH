#include "header.h"
#include "parametrs.h"
#include "particles.h"

// Реализация вспомогательных функций
void FWQ(double QQ, double HS, double& WQ) {
    if (QQ <= 1.0) {
        WQ = 1.0 - 1.5 * QQ * QQ + 0.75 * QQ * QQ * QQ;
    } else if (QQ <= 2.0) {
        double QM2 = 2.0 - QQ;
        WQ = 0.25 * QM2 * QM2 * QM2;
    } else {
        WQ = 0.0;
    }
    // WQ *= DKOEFW[NDIM] / pow(HS, NDIM); // DKOEFW должен быть определен
}

void FDWDQ(double QQ, double HS, double& DWDQ) {
    if (QQ <= 1.0) {
        DWDQ = -3.0 * QQ + 2.25 * QQ * QQ;
    } else if (QQ <= 2.0) {
        double QM2 = 2.0 - QQ;
        DWDQ = -0.75 * QM2 * QM2;
    } else {
        DWDQ = 0.0;
    }
    // DWDQ *= DKOEFW[NDIM] / (QQ * pow(HS, NDIM+2)); // DKOEFW должен быть определен
}

void FD2WDQ2(double QQ, double HS, double& D2WDQ2) {
    if (QQ <= 1.0) {
        D2WDQ2 = -3.0 + 4.5 * QQ;
    } else if (QQ <= 2.0) {
        double QM2 = 2.0 - QQ;
        D2WDQ2 = 1.5 * QM2;
    } else {
        D2WDQ2 = 0.0;
    }
    // D2WDQ2 *= DKOEFW[NDIM] / pow(HS, NDIM+2); // DKOEFW должен быть определен
}
//!===================================================== =
double FWQ(double QQ, double HS)//subroutine FWQ(QQ, HS, WQ)
{
    //use SPHEP_MEM
    //implicit none
    VARIABLE V;
    double QM2;
    double WQ;
    if (QQ <= 1.e0) //then
    {
        WQ = 1.e0 - 1.5e0 * QQ * QQ + 0.75e0 * QQ * QQ * QQ;
    }
    else if (QQ <= 2.e0) //then
    {
        QM2 = 2.e0 - QQ;
        WQ = 0.25e0 * QM2 * QM2 * QM2;
    }
    else
    {
        WQ = 0.e0;

    }//endif
    
    WQ = WQ * V.DKOEFW[V.NDIM-1] / (pow(HS, V.NDIM));
    //std::cout << "QQ=" << QQ << "   WQ=" << WQ << "   HS=" << HS << std::endl;
    //getchar();
    //!DKOEFW = (DKOEFW0 / HS * *NDIM) / HS * *2;
    return WQ;
}//end
//!===================================================== =

//!===================================================== =
double FDWDQ(double QQ, double HS) //subroutine FDWDQ(QQ, HS, DWDQ)
{
    //use SPHEP_MEM
    //implicit none
    VARIABLE V;
    double DWDQ;
    double QM2{1.00};
    if (QQ <= 1.e0)//then
    {
        DWDQ = -3.e0 * QQ + 2.25e0 * QQ * QQ;
    }
    else if (QQ <= 2.e0) //then
    {
        QM2 = 2.e0 - QQ;
        DWDQ = -0.75e0 * QM2 * QM2;
    }
    else
    {
        DWDQ = 0.e0;
    }//endif
    DWDQ = DWDQ * V.DKOEFW[V.NDIM-1] / (QQ * pow(HS, (V.NDIM + 2)));
    //std::cout << "QQ=" << QQ << "   QM=" << QM2 << "   HS=" << HS << "   DWDQ=" << DWDQ << std::endl;
    //getchar();
    return DWDQ;

}//end;!===================================================== =

//!===================================================== =
double FD2WDQ2(double QQ, double HS)//subroutine FD2WDQ2(QQ, HS, DWDQ)
{
    //use SPHEP_MEM
    //implicit none;
    VARIABLE V;
    double DWDQ;
    double QM2;
    //!
    if (QQ <= 1.e0)//then
    {
        DWDQ = -3.e0 + 4.5e0 * QQ;
    }
    else if (QQ <= 2.e0) //then
    {
        QM2 = 2.e0 - QQ;
        DWDQ = 1.5e0 * QM2;
    }
    else
    {
        DWDQ = 0.e0;
    }//endif
    DWDQ = DWDQ * V.DKOEFW[V.NDIM-1] / pow(HS, (V.NDIM + 2));
    return DWDQ;
    //end;
}//!===================================================== =

void acceleration(int &I, int &J, std::vector<particles> &p, double &MNO1, double &ART, double &DWDQ)
{
    //acceleration
    for (int ALF = 0; ALF < 3; ALF++) //do ALF = 1, NDIM
    {
        for (int BET = 0; BET < 3; BET++) //do BET = 1, NDIM
        {
            MNO1 = 0.e0;
            
            if (ALF == 0 && BET == 0) MNO1 = p[J].IMAS * ((p[I].IP - p[I].ISXX + p[J].IP - p[J].ISXX) / (p[I].IDNS * p[J].IDNS) + ART);
            
            if (ALF == 1 && BET == 1) MNO1 = p[J].IMAS * ((p[I].IP - p[I].ISYY + p[J].IP - p[J].ISYY) / (p[I].IDNS * p[J].IDNS) + ART);
            
            if (ALF == 2 && BET == 2) MNO1 = p[J].IMAS * ((p[I].IP - p[I].ISZZ + p[J].IP - p[J].ISZZ) / (p[I].IDNS * p[J].IDNS) + ART);
            
            if (ALF == 0 && BET == 1) MNO1 = p[J].IMAS * ((-p[I].ISXY - p[J].ISXY) / (p[I].IDNS * p[J].IDNS));
            
            if (ALF == 1 && BET == 0) MNO1 = p[J].IMAS * ((-p[I].ISXY - p[J].ISXY) / (p[I].IDNS * p[J].IDNS));
            
            if (ALF == 0 && BET == 2) MNO1 = p[J].IMAS * ((-p[I].ISXZ - p[J].ISXZ) / (p[I].IDNS * p[J].IDNS));
            
            if (ALF == 2 && BET == 0) MNO1 = p[J].IMAS * ((-p[I].ISXZ - p[J].ISXZ) / (p[I].IDNS * p[J].IDNS));
            
            if (ALF == 1 && BET == 2) MNO1 = p[J].IMAS * ((-p[I].ISYZ - p[J].ISYZ) / (p[I].IDNS * p[J].IDNS));
            
            if (ALF == 2 && BET == 1) MNO1 = p[J].IMAS * ((-p[I].ISYZ - p[J].ISYZ) / (p[I].IDNS * p[J].IDNS));
            
            p[I].IACS[ALF] -= MNO1 * DWDQ * (*p[I].IXX_Ptr[BET] - *p[J].IXX_Ptr[BET]);
            
            p[I].IDU += 0.5e0 * MNO1 * DWDQ * (*p[I].IVV_Ptr[ALF] - *p[J].IVV_Ptr[ALF]) * (*p[I].IXX_Ptr[BET] - *p[J].IXX_Ptr[BET]);
        }
    }
}

void rotation_rates(int &I, int &J, std::vector<particles> &p, double &DWDQ)
{
    //rotation rates
    p[I].IRXY += 0.5e0 * ((*p[I].IVV_Ptr[0] - *p[J].IVV_Ptr[0]) * (*p[I].IXX_Ptr[1] - *p[J].IXX_Ptr[1]) -
        (*p[I].IVV_Ptr[1] - *p[J].IVV_Ptr[1]) * (*p[I].IXX_Ptr[0] - *p[J].IXX_Ptr[0])) * DWDQ * p[J].IMAS / (p[J].IDNS);

    p[I].IRXZ += 0.5e0 * ((*p[I].IVV_Ptr[0] - *p[J].IVV_Ptr[0]) * (*p[I].IXX_Ptr[2] - *p[J].IXX_Ptr[2]) -
        (*p[I].IVV_Ptr[2] - *p[J].IVV_Ptr[2]) * (*p[I].IXX_Ptr[0] - *p[J].IXX_Ptr[0])) * DWDQ * p[J].IMAS / (p[J].IDNS);

    p[I].IRYZ += 0.5e0 * ((*p[I].IVV_Ptr[1] - *p[J].IVV_Ptr[1]) * (*p[I].IXX_Ptr[2] - *p[J].IXX_Ptr[2]) -
        (*p[I].IVV_Ptr[2] - *p[J].IVV_Ptr[2]) * (*p[I].IXX_Ptr[1] - *p[J].IXX_Ptr[1])) * DWDQ * p[J].IMAS / (p[J].IDNS);
}

void macroscopic_deformation(int &I, int &J, std::vector<particles> &p, double &DWDQ, double &DTAU)
{
    //macroscopic deformation
    p[I].IUXY -= 0.5e0 * ((*p[I].IVV_Ptr[0] - *p[J].IVV_Ptr[0]) *
        (*p[I].IXX_Ptr[1] - *p[J].IXX_Ptr[1]) + (*p[I].IVV_Ptr[1]) - *p[J].IVV_Ptr[1] *
        (*p[I].IXX_Ptr[0] - *p[J].IXX_Ptr[0])) * DWDQ * DTAU * p[J].IMAS / (p[J].IDNS);

    p[I].IUXZ -= 0.5e0 * ((*p[I].IVV_Ptr[0] - *p[J].IVV_Ptr[0]) *
        (*p[I].IXX_Ptr[2] - *p[J].IXX_Ptr[2]) + (*p[I].IVV_Ptr[2] - *p[J].IVV_Ptr[2]) *
        (*p[I].IXX_Ptr[0] - *p[J].IXX_Ptr[0])) * DWDQ * DTAU * p[J].IMAS / (p[J].IDNS);

    p[I].IUYZ -= 0.5e0 * ((*p[I].IVV_Ptr[1]) - *p[J].IVV_Ptr[1] *
        (*p[I].IXX_Ptr[2] - *p[J].IXX_Ptr[2]) + (*p[I].IVV_Ptr[2] - *p[J].IVV_Ptr[2]) *
        (*p[I].IXX_Ptr[1] - *p[J].IXX_Ptr[1])) * DWDQ * DTAU * p[J].IMAS / (p[J].IDNS);

    p[I].IUXX -= 0.5e0 * ((*p[I].IVV_Ptr[0] - *p[J].IVV_Ptr[0]) *
        (*p[I].IXX_Ptr[0] - *p[J].IXX_Ptr[0]) + (*p[I].IVV_Ptr[0] - *p[J].IVV_Ptr[0]) *
        (*p[I].IXX_Ptr[0] - *p[J].IXX_Ptr[0])) * DWDQ * DTAU * p[J].IMAS / (p[J].IDNS);

    p[I].IUYY -= 0.5e0 * ((*p[I].IVV_Ptr[1] - *p[J].IVV_Ptr[1]) *
        (*p[I].IXX_Ptr[1] - *p[J].IXX_Ptr[1]) + (*p[I].IVV_Ptr[1] - *p[J].IVV_Ptr[1]) *
        (*p[I].IXX_Ptr[1] - *p[J].IXX_Ptr[1])) * DWDQ * DTAU * p[J].IMAS / (p[J].IDNS);

    p[I].IUZZ -= 0.5e0 * ((*p[I].IVV_Ptr[2] - *p[J].IVV_Ptr[2]) *
        (*p[I].IXX_Ptr[2] - *p[J].IXX_Ptr[2]) + (*p[I].IVV_Ptr[2] - *p[J].IVV_Ptr[2]) *
        (*p[I].IXX_Ptr[2] - *p[J].IXX_Ptr[2])) * DWDQ * DTAU * p[J].IMAS / (p[J].IDNS);
}

int MAX(int a, int b)
{
    if (a == b || a > b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

int MIN(int a, int b)
{
    if (a < b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

//!=======================================================================BEG_SPHEP_MOVE============================================================================== =
void MOVE(double TIME, std::vector<particles> &particle) 
{
    parametrs param;
    double RR2, RR, QQ, HS;
    double MNO1,  SUM,  DWDQ, WQ;
    int NUM, NN;
    double MIU, ART;
    //double DX_RESCALE, DX_SHIFT;

    double DTAU = abs(particle[0].IHS / particle[0].ICS);
    for (const auto &p: particle) 
    {
        DTAU = (abs(p.IHS / p.ICS) < DTAU) ? abs(p.IHS < p.ICS) : DTAU; //DTAU1 = abs(p.IHS / p.ICS); //if (DTAU1 < DTAU) DTAU = DTAU1;
    }
    DTAU *= param.COEF_DTAU;

    for (int I = 0; I < particle.size(); I++) //do I = 1, NPT
    {
        particle[I].IRXY = 0.e0;
        particle[I].IRXZ = 0.e0;
        particle[I].IRYZ = 0.e0;
        for (int ALF = 0; ALF < 3; ALF++) particle[I].IACS[ALF] = 0.e0;
        particle[I].IDDNS = 0.e0;
        particle[I].IDU = 0.e0;
        particle[I].RHOO = 0.e0;
        //FS[I][IWSQ] = 0.e0;
        
        // Получаем координаты ячейки для текущей частицы
        int NX = particle[I].MESH[0];   //MESH[I][IX];
        int NY = particle[I].MESH[1];   //MESH[I][IY];
        int NZ = particle[I].MESH[2];   //MESH[I][IZ];
        int NQ = 0;

        /*// Цикл по соседним ячейкам (3x3x3 область)
        for (int NX1 = std::max(NX-1, 1); NX1 <= std::min(NX+1, NMES0); ++NX1) {
            for (int NY1 = std::max(NY-1, 1); NY1 <= std::min(NY+1, NMES0); ++NY1) {
                for (int NZ1 = std::max(NZ-1, 1); NZ1 <= std::min(NZ+1, NMES0); ++NZ1) {
                    int NUM = NPAT[NX1-1][NY1-1][NZ1-1][0]; // -1 из-за индексации с 0 в C++
                    
                    // Цикл по частицам в текущей ячейке
                    for (int NN = 1; NN <= NUM; ++NN) {
                        int J = NPAT[NX1-1][NY1-1][NZ1-1][NN];*/
        //for (int NX1 = MAX((NX - 1), 0); NX1 < MIN((NX + 1), value.NMES0); NX1++) //do NX1 = max(NX - 1, 1), min(NX + 1, NMES0)
        for (int NX1 = MAX((NX - 1), 1); NX1 < MIN((NX + 1), param.NMES0); NX1++) //do NX1 = max(NX - 1, 1), min(NX + 1, NMES0)
        {
            //for (int NY1 = MAX((NY - 1), 0); NY1 < MIN((NY + 1), value.NMES0); NY1++) //do NY1 = max(NY - 1, 1), min(NY + 1, NMES0)
            for (int NY1 = MAX((NY - 1), 1); NY1 < MIN((NY + 1), param.NMES0); NY1++) //do NY1 = max(NY - 1, 1), min(NY + 1, NMES0)
            {
                //for (int NZ1 = MAX((NZ - 1), 0); NZ1 < MIN((NZ + 1), value.NMES0); NZ1++) //do NZ1 = max(NZ - 1, 1), min(NZ + 1, NMES0)
                for (int NZ1 = MAX((NZ - 1), 1); NZ1 < MIN((NZ + 1), param.NMES0); NZ1++) //do NZ1 = max(NZ - 1, 1), min(NZ + 1, NMES0)
                {
                    NUM = NPAT[NX1][NY1][NZ1][0];
                    for (int NN = 1; NN <= NUM; NN++)   //do NN = 1, NUM
                    {
                        int J = NPAT[NX1][NY1][NZ1][NN];
                        if (J != I)//then
                        {
                            RR2 = pow((particle[I].IX - particle[J].IX), 2) + 
                                        pow((particle[I].IY - particle[J].IY), 2) + 
                                            pow((particle[I].IZ - particle[J].IZ), 2);
                            HS = particle[J].IHS;
                            if (RR2 < 4.e0 * HS * HS) //then
                            {
                                NQ = NQ + 1;
                                RR = sqrt(RR2);
                                QQ = RR / HS; //!*HS_Inv;
                                WQ = FWQ(QQ, HS);
                                DWDQ = FDWDQ(QQ, HS);
                                particle[I].RHOO += particle[J].IMAS * WQ;

                                //velocity divergence;
                                SUM = 0.e0;
                                for (int ALF = 0; ALF < 3; ALF++) //do ALF = 1, NDIM
                                {
                                    SUM += (*particle[I].IVV_Ptr[ALF] - *particle[J].IVV_Ptr[ALF]) * (*particle[I].IXX_Ptr[ALF] - *particle[J].IXX_Ptr[ALF]);
                                }
                                
                                //artificial viscosity;
                                if (SUM < 0.e0)
                                {
                                    MIU = 3.e0 * SUM / sqrt(RR2 + 0.01e0 * HS * HS) * 2.e0;
                                    ART = MIU * (-(particle[I].ICS + particle[J].ICS) + 4.e0 * MIU) / (particle[I].IDNS + particle[J].IDNS);
                                }
                                else
                                {
                                    ART = 0.e0;
                                }

                                particle[I].IDDNS += particle[J].IMAS * SUM * DWDQ;

                                acceleration(I, J, particle, MNO1, ART, DWDQ);

                                rotation_rates(I, J, particle, DWDQ);
                                
                                macroscopic_deformation(I, J, particle, DWDQ, DTAU);
                            }
                        }
                    }
                }
            }
        }
    }


    //!FS = FSN обновление компонент
    for (auto &p: particle) //do I = 1, NPT
    {
        //!change of Wij components at rotation
        p.IWXX += (-2.e0 * p.IWXY * p.IRXY - 2.e0 * p.IWXZ * p.IRXZ) * DTAU;
        p.IWYY +=  (2.e0 * p.IWXY * p.IRXY - 2.e0 * p.IWYZ * p.IRYZ) * DTAU;
        p.IWZZ +=  (2.e0 * p.IWXZ * p.IRXZ + 2.e0 * p.IWYZ * p.IRYZ) * DTAU;
        p.IWXY += ((p.IWXX - p.IWYY) * p.IRXY - p.IWXZ * p.IRYZ - p.IWYZ * p.IRXZ) * DTAU;
        p.IWXZ += ((p.IWXX - p.IWZZ) * p.IRXZ + p.IWXY * p.IRYZ - p.IWYZ * p.IRXY) * DTAU;
        p.IWYZ += ((p.IWYY - p.IWZZ) * p.IRYZ + p.IWXY * p.IRXZ - p.IWXZ * p.IRXY) * DTAU;
        //!change of Uij components at rotation
        p.IUXX += (-2.e0 * p.IUXY * p.IRXY - 2.e0 * p.IUXZ * p.IRXZ) * DTAU;
        p.IUYY +=  (2.e0 * p.IUXY * p.IRXY - 2.e0 * p.IUYZ * p.IRYZ) * DTAU;
        p.IUZZ +=  (2.e0 * p.IUXZ * p.IRXZ + 2.e0 * p.IUYZ * p.IRYZ) * DTAU;
        p.IUXY = p.IUXY + ((p.IUXX - p.IUYY) * p.IRXY - p.IUXZ * p.IRYZ - p.IUYZ * p.IRXZ) * DTAU; 
        p.IUXZ = p.IUXZ + ((p.IUXX - p.IUZZ) * p.IRXZ + p.IUXY * p.IRYZ - p.IUYZ * p.IRXY) * DTAU;
        p.IUYZ = p.IUYZ + ((p.IUYY - p.IUZZ) * p.IRYZ + p.IUXY * p.IRXZ - p.IUXZ * p.IRXY) * DTAU;
        //!stress deviators
        p.IULL = p.IUXX + p.IUYY + p.IUZZ;
        p.ISXY = 2.e0 * p.IG * (p.IUXY - p.IWXY);
        p.ISXZ = 2.e0 * p.IG * (p.IUXZ - p.IWXZ);
        p.ISYZ = 2.e0 * p.IG * (p.IUYZ - p.IWYZ);
        p.ISXX = 2.e0 * p.IG * (p.IUXX - param.C1D3 * p.IULL - p.IWXX);
        p.ISYY = 2.e0 * p.IG * (p.IUYY - param.C1D3 * p.IULL - p.IWYY);
        p.ISZZ = 2.e0 * p.IG * (p.IUZZ - param.C1D3 * p.IULL - p.IWZZ);
        p.IMisSTR = sqrt(pow((p.ISXX - p.ISYY), 2) + pow((p.ISYY - p.ISZZ), 2) + pow((p.ISZZ - p.ISXX), 2) +
                                6.e0 * (pow(p.ISXY, 2) + pow(p.ISXZ, 2) + pow(p.ISYZ, 2))) / sqrt(2.e0);
        p.IDNS = p.IDNS + p.IDDNS * DTAU;
        p.IU = p.IU + p.IDU * DTAU;    //!internal energy

        //EOS BRASS
        
        EOS_KH(p.IDNS, p.IU, p.IT, p.IP, p.ICS, p.ICV, param.SUBZ);
        
        p.IKK = p.IDNS * pow(p.ICS, 2);
        p.IG = 1.5e0 * p.IKK * (1.e0 - 2.e0 * p.ICP) / (1.e0 + p.ICP);
        
        for (int ALF = 0; ALF < 3; ALF++)
        {
            p.IDX[ALF] += *p.IVV_Ptr[ALF] * DTAU;
            *p.IXX_Ptr[ALF] += *p.IVV_Ptr[ALF] * DTAU;
            *p.IVREALV_Ptr[ALF] += p.IACS[ALF] * DTAU;
        }
        //!Rescale
        /*for (int ALF = 0; ALF < 3; ALF++) //do ALF = 1, 3 //num_axi
        {
            DX_RESCALE = -FS[I][V.IXX[ALF]] * V.vepsXYZ[ALF] * DTAU;
            FS[I][V.IDX[ALF]] = FS[I][V.IDX[ALF]] + DX_RESCALE;
            FS[I][V.IVV[ALF]] = FS[I][V.IVREALV[ALF]] + DX_RESCALE / DTAU;
            FS[I][V.IXX[ALF]] = FS[I][V.IXX[ALF]] + DX_RESCALE;
        }*/
        /*
        //!Boundary condition
        for (int ALF = 0; ALF < 3; ALF++) //do ALF = 1, 3 //num_axi   TRI - Axial
        {
            if (FS[I][V.IXX[ALF]] < CMINMAX[ALF][0]) 
            {
                DX_SHIFT = CMINMAX[ALF][0] - FS[I][V.IXX[ALF]];
                FS[I][V.IDX[ALF]] = FS[I][V.IDX[ALF]] + DX_SHIFT;
                FS[I][V.IVV[ALF]] = FS[I][V.IVV[ALF]] + DX_SHIFT / DTAU;
                FS[I][V.IXX[ALF]] = CMINMAX[ALF][0];
            }
            if (FS[I][V.IXX[ALF]] > CMINMAX[ALF][1]) 
            {
                DX_SHIFT = CMINMAX[ALF][1] - FS[I][V.IXX[ALF]];
                FS[I][V.IDX[ALF]] = FS[I][V.IDX[ALF]] + DX_SHIFT;
                FS[I][V.IVV[ALF]] = FS[I][V.IVV[ALF]] + DX_SHIFT / DTAU;
                FS[I][V.IXX[ALF]] = CMINMAX[ALF][1];
            }
        }
        */
    }

    TIME = TIME + DTAU;
    return TIME;
}//end subroutine SPHEP_MOVE
//!=====================================END_SPHEP_MOVE================================================ =