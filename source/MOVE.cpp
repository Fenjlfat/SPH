#include "header.h"
#include "parametrs.h"
#include "particles.h"


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
    double SUBZ = 0.e0;
    double DNS_EOS = 0.00;
    double T_EOS = 0.00;
    double E_EOS = 0.00;
    double P_EOS = 0.00;
    double CS_EOS = 0.00;
    double CV_EOS = 0.00;

    double RR2, RR, QQ, HS;
    double MNO1,  SUM,  DWDQ, SUM1, WQ, MNO2;
    double NX, NY, NZ;
    int NUM, NN, NQ;
    double MIU, ART;
    double DX_RESCALE, DX_SHIFT;

    double DTAU = abs(particle[0].IHS / particle[0].ICS);
    for (const auto &p: particle) 
    {
        DTAU = (abs(p.IHS / p.ICS) < DTAU) ? abs(p.IHS < p.ICS) : DTAU; //DTAU1 = abs(p.IHS / p.ICS); //if (DTAU1 < DTAU) DTAU = DTAU1;
    }
    DTAU *= param.COEF_DTAU;

    for (int I = 0; I < NPT; I++) //do I = 1, NPT
    {
        FS[I][V.IRXY] = 0.e0;
        FS[I][V.IRXZ] = 0.e0;
        FS[I][V.IRYZ] = 0.e0;
        for (int ALF = 0; ALF < V.NDIM; ALF++) //do ALF = 1, NDIM;
        {
            FS[I][V.IACS[ALF]] = 0.e0;
        }//enddo
        FS[I][V.IDDNS] = 0.e0;
        FS[I][V.IDU] = 0.e0;
        FS[I][V.RHOO] = 0.e0;
        //FS[I][IWSQ] = 0.e0;
        //!Ã¶ÃšÃªÃ« Ã¯Ã® Ã±Ã®Ã±Ã¥Ã€Ã¿Ã¬
        NX = value.MESH[I][V.IX];
        NY = value.MESH[I][V.IY];
        NZ = value.MESH[I][V.IZ];
        NQ = 0;
        
        //for (int NX1 = MAX((NX - 1), 0); NX1 < MIN((NX + 1), value.NMES0); NX1++) //do NX1 = max(NX - 1, 1), min(NX + 1, NMES0)
        for (int NX1 = MAX((NX - 1), 1); NX1 < MIN((NX + 1), value.NMES0); NX1++) //do NX1 = max(NX - 1, 1), min(NX + 1, NMES0)
        {
            //for (int NY1 = MAX((NY - 1), 0); NY1 < MIN((NY + 1), value.NMES0); NY1++) //do NY1 = max(NY - 1, 1), min(NY + 1, NMES0)
            for (int NY1 = MAX((NY - 1), 1); NY1 < MIN((NY + 1), value.NMES0); NY1++) //do NY1 = max(NY - 1, 1), min(NY + 1, NMES0)
            {
                //for (int NZ1 = MAX((NZ - 1), 0); NZ1 < MIN((NZ + 1), value.NMES0); NZ1++) //do NZ1 = max(NZ - 1, 1), min(NZ + 1, NMES0)
                for (int NZ1 = MAX((NZ - 1), 1); NZ1 < MIN((NZ + 1), value.NMES0); NZ1++) //do NZ1 = max(NZ - 1, 1), min(NZ + 1, NMES0)
                {
                    NUM = value.NPAT[NX1][NY1][NZ1][0];
                    //std::cout << "NUM=" << NUM << std::endl;
                    //getchar();
                    for (int NN = 1; NN <= NUM; NN++)   //do NN = 1, NUM
                    {
                        int J = value.NPAT[NX1][NY1][NZ1][NN];
                        if (J != I)//then
                        {
                            RR2 = pow((FS[I][V.IX] - FS[J][V.IX]), 2) + pow((FS[I][V.IY] - FS[J][V.IY]), 2) + pow((FS[I][V.IZ] - FS[J][V.IZ]), 2);
                            IVIX = FS[I][V.IX];
                            JVIX = FS[J][V.IX];
                            HS = FS[J][V.IHS];
                            if (RR2 < 4.e0 * HS * HS) //then
                            {
                                NQ = NQ + 1;
                                RR = sqrt(RR2);
                                QQ = RR / HS; //!*HS_Inv;
                                WQ = FWQ(QQ, HS);
                                //std::cout << "WQ=" << WQ << "   QQ=" << QQ << "   HS=" << HS << std::endl;
                                //getchar();
                                DWDQ = FDWDQ(QQ, HS);
                                //std::cout << "DWDQ=" << DWDQ << std::endl;
                                //getchar();
                                //DWDQ = FD2WDQ2(QQ, HS);   /////????????????????????????
                                //std::cout << "DWDQ=" << DWDQ << std::endl;
                                //getchar();
                                FS[I][V.RHOO] = FS[I][V.RHOO] + FS[J][V.IMAS] * WQ;

                                //!velocity divergence;
                                SUM = 0.e0;
                                for (int ALF = 0; ALF < V.NDIM; ALF++) //do ALF = 1, NDIM
                                {
                                    SUM = SUM + (FS[I][V.IVV[ALF]] - FS[J][V.IVV[ALF]]) * (FS[I][V.IXX[ALF]] - FS[J][V.IXX[ALF]]);
                                }//enddo;
                                
                                //!artificial viscosity;
                                if (SUM < 0.e0) //then
                                {
                                    MIU = 3.e0 * SUM / sqrt(RR2 + 0.01e0 * HS * HS) * 2.e0;
                                    ART = MIU * (-(FS[I][V.ICS] + FS[J][V.ICS]) + 4.e0 * MIU) / (FS[I][V.IDNS] + FS[J][V.IDNS]);
                                }
                                else
                                {
                                    ART = 0.e0;
                                }//endif

                                //!acceleration
                                for (int ALF = 0; ALF < V.NDIM; ALF++) //do ALF = 1, NDIM
                                {
                                    for (int BET = 0; BET < V.NDIM; BET++) //do BET = 1, NDIM
                                    {
                                        MNO1 = 0.e0;
                                        //std::cout << "MNO1=" << MNO1 << std::endl;
                                        if (ALF == 0 && BET == 0) 
                                            MNO1 = FS[J][V.IMAS] * ((FS[I][V.IP] - FS[I][V.ISXX] + FS[J][V.IP] - FS[J][V.ISXX]) / (FS[I][V.IDNS] * FS[J][V.IDNS]) + ART);
                                        
                                        if (ALF == 1 && BET == 1) 
                                            MNO1 = FS[J][V.IMAS] * ((FS[I][V.IP] - FS[I][V.ISYY] + FS[J][V.IP] - FS[J][V.ISYY]) / (FS[I][V.IDNS] * FS[J][V.IDNS]) + ART);
                                        
                                        if (ALF == 2 && BET == 2) 
                                            MNO1 = FS[J][V.IMAS] * ((FS[I][V.IP] - FS[I][V.ISZZ] + FS[J][V.IP] - FS[J][V.ISZZ]) / (FS[I][V.IDNS] * FS[J][V.IDNS]) + ART);
                                        
                                        if (ALF == 0 && BET == 1) 
                                            MNO1 = FS[J][V.IMAS] * ((-FS[I][V.ISXY] - FS[J][V.ISXY]) / (FS[I][V.IDNS] * FS[J][V.IDNS]));
                                        
                                        if (ALF == 1 && BET == 0) 
                                            MNO1 = FS[J][V.IMAS] * ((-FS[I][V.ISXY] - FS[J][V.ISXY]) / (FS[I][V.IDNS] * FS[J][V.IDNS]));
                                        
                                        if (ALF == 0 && BET == 2) 
                                            MNO1 = FS[J][V.IMAS] * ((-FS[I][V.ISXZ] - FS[J][V.ISXZ]) / (FS[I][V.IDNS] * FS[J][V.IDNS]));
                                        
                                        if (ALF == 2 && BET == 0) 
                                            MNO1 = FS[J][V.IMAS] * ((-FS[I][V.ISXZ] - FS[J][V.ISXZ]) / (FS[I][V.IDNS] * FS[J][V.IDNS]));
                                        
                                        if (ALF == 1 && BET == 2) 
                                            MNO1 = FS[J][V.IMAS] * ((-FS[I][V.ISYZ] - FS[J][V.ISYZ]) / (FS[I][V.IDNS] * FS[J][V.IDNS]));
                                        
                                        if (ALF == 2 && BET == 1) 
                                            MNO1 = FS[J][V.IMAS] * ((-FS[I][V.ISYZ] - FS[J][V.ISYZ]) / (FS[I][V.IDNS] * FS[J][V.IDNS]));
                                        
                                        
                                        FS[I][V.IACS[ALF]] = FS[I][V.IACS[ALF]] - MNO1 * DWDQ * (FS[I][V.IXX[BET]] - FS[J][V.IXX[BET]]);
                                        /*if (I == 0 && ALF == 1)
                                        {
                                            std::cout << "ALF=" << ALF << "   BET=" << BET << std::endl;
                                            std::cout << "MNO1=" << MNO1 << "   DWDQ=" << DWDQ << std::endl;
                                            std::cout << "FS["<<I<<"][V.IXX[BET]] = " << FS[I][V.IXX[BET]] << "    FS[" << J << "][V.IXX[BET]] = " << FS[J][V.IXX[BET]] << std::endl;
                                            std::cout << "J=" << J << "   NPAT["<<NX1<<"]["<<NY1<<"]["<<NZ1<<"]["<<NN<<"]" << std::endl;
                                            getchar();
                                            //V.MNO1_0_fix = MNO1;
                                            //V.DWDQ_FIX = DWDQ;
                                            //IVIX=FS[I][V.IXX[BET]];
                                            //JVIX=FS[J][V.IXX[BET]];

                                        }*/
                                        /*if (I == 0 && ALF == 1)
                                        {
                                            std::cout << "FS[I][V.IACS[" << ALF << "]]=" << FS[I][V.IACS[ALF]] << std::endl;
                                            std::cout << "FS[I][V.IXX[" << BET << "]]=" << FS[I][V.IXX[BET]] << std::endl;
                                            std::cout << "FS[J][V.IXX[" << BET << "]]=" << FS[J][V.IXX[BET]] << std::endl;
                                            std::cout << "DWDQ=" << DWDQ << std::endl;
                                            std::cout << "MNO1=" << MNO1 << std::endl;
                                            //std::cout << "MNO1=" << std::setprecision(25) << MNO1 << std::endl;
                                            std::cout << "I=" << I << std::endl;
                                            std::cout << "J=" << J<< std::endl;
                                            std::cout << "ALF=" << ALF<< std::endl;
                                            getchar();
                                        }*/
                                        
                                        FS[I][V.IDU] = FS[I][V.IDU] + 0.5e0 * MNO1 * DWDQ *
                                            (FS[I][V.IVV[ALF]] - FS[J][V.IVV[ALF]]) * (FS[I][V.IXX[BET]] - FS[J][V.IXX[BET]]);
                                    }//enddo
                                }//enddo

                                FS[I][V.IDDNS] = FS[I][V.IDDNS] + FS[J][V.IMAS] * SUM * DWDQ;
                                //std::cout << DWDQ << std::endl;
                                //getchar();
                                // 
                                //!rotation rates
                                FS[I][V.IRXY] = FS[I][V.IRXY] + 0.5e0 * ((FS[I][V.IVV[0]] - FS[J][V.IVV[0]]) * (FS[I][V.IXX[1]] - FS[J][V.IXX[1]]) -
                                    (FS[I][V.IVV[1]] - FS[J][V.IVV[1]]) * (FS[I][V.IXX[0]] - FS[J][V.IXX[0]])) * DWDQ * FS[J][V.IMAS] / (FS[J][V.IDNS]);

                                FS[I][V.IRXZ] = FS[I][V.IRXZ] + 0.5e0 * ((FS[I][V.IVV[0]] - FS[J][V.IVV[0]]) * (FS[I][V.IXX[2]] - FS[J][V.IXX[2]]) -
                                    (FS[I][V.IVV[2]] - FS[J][V.IVV[2]]) * (FS[I][V.IXX[0]] - FS[J][V.IXX[0]])) * DWDQ * FS[J][V.IMAS] / (FS[J][V.IDNS]);

                                FS[I][V.IRYZ] = FS[I][V.IRYZ] + 0.5e0 * ((FS[I][V.IVV[1]] - FS[J][V.IVV[1]]) * (FS[I][V.IXX[2]] - FS[J][V.IXX[2]]) -
                                    (FS[I][V.IVV[2]] - FS[J][V.IVV[2]]) * (FS[I][V.IXX[1]] - FS[J][V.IXX[1]])) * DWDQ * FS[J][V.IMAS] / (FS[J][V.IDNS]);
                                //!macroscopic deformation
                                FS[I][V.IUXY] = FS[I][V.IUXY] - 0.5e0 * ((FS[I][V.IVV[0]] - FS[J][V.IVV[0]]) *
                                    (FS[I][V.IXX[1]] - FS[J][V.IXX[1]]) + (FS[I][V.IVV[1]]) - FS[J][V.IVV[1]] *
                                    (FS[I][V.IXX[0]] - FS[J][V.IXX[0]])) * DWDQ * DTAU * FS[J][V.IMAS] / (FS[J][V.IDNS]);

                                FS[I][V.IUXZ] = FS[I][V.IUXZ] - 0.5e0 * ((FS[I][V.IVV[0]] - FS[J][V.IVV[0]]) *
                                    (FS[I][V.IXX[2]] - FS[J][V.IXX[2]]) + (FS[I][V.IVV[2]] - FS[J][V.IVV[2]]) *
                                    (FS[I][V.IXX[0]] - FS[J][V.IXX[0]])) * DWDQ * DTAU * FS[J][V.IMAS] / (FS[J][V.IDNS]);

                                FS[I][V.IUYZ] = FS[I][V.IUYZ] - 0.5e0 * ((FS[I][V.IVV[1]]) - FS[J][V.IVV[1]] *
                                    (FS[I][V.IXX[2]] - FS[J][V.IXX[2]]) + (FS[I][V.IVV[2]] - FS[J][V.IVV[2]]) *
                                    (FS[I][V.IXX[1]] - FS[J][V.IXX[1]])) * DWDQ * DTAU * FS[J][V.IMAS] / (FS[J][V.IDNS]);

                                FS[I][V.IUXX] = FS[I][V.IUXX] - 0.5e0 * ((FS[I][V.IVV[0]] - FS[J][V.IVV[0]]) *
                                    (FS[I][V.IXX[0]] - FS[J][V.IXX[0]]) + (FS[I][V.IVV[0]] - FS[J][V.IVV[0]]) *
                                    (FS[I][V.IXX[0]] - FS[J][V.IXX[0]])) * DWDQ * DTAU * FS[J][V.IMAS] / (FS[J][V.IDNS]);

                                FS[I][V.IUYY] = FS[I][V.IUYY] - 0.5e0 * ((FS[I][V.IVV[1]] - FS[J][V.IVV[1]]) *
                                    (FS[I][V.IXX[1]] - FS[J][V.IXX[1]]) + (FS[I][V.IVV[1]] - FS[J][V.IVV[1]]) *
                                    (FS[I][V.IXX[1]] - FS[J][V.IXX[1]])) * DWDQ * DTAU * FS[J][V.IMAS] / (FS[J][V.IDNS]);

                                FS[I][V.IUZZ] = FS[I][V.IUZZ] - 0.5e0 * ((FS[I][V.IVV[2]] - FS[J][V.IVV[2]]) *
                                    (FS[I][V.IXX[2]] - FS[J][V.IXX[2]]) + (FS[I][V.IVV[2]] - FS[J][V.IVV[2]]) *
                                    (FS[I][V.IXX[2]] - FS[J][V.IXX[2]])) * DWDQ * DTAU * FS[J][V.IMAS] / (FS[J][V.IDNS]);

                            }//endif;
                        }//endif;
                    }//enddo;
                }//enddo;
            }//enddo;
        }//enddo
    }//enddo //!(I)


    //!FS = FSN
    for (int I = 0; I < NPT; I++) //do I = 1, NPT
    {
        //!change of Wij components at rotation
        FS[I][V.IWXX] = FS[I][V.IWXX] + (-2.e0 * FS[I][V.IWXY] * FS[I][V.IRXY] - 2.e0 * FS[I][V.IWXZ] * FS[I][V.IRXZ]) * DTAU;
        FS[I][V.IWYY] = FS[I][V.IWYY] + (2.e0 * FS[I][V.IWXY] * FS[I][V.IRXY] - 2.e0 * FS[I][V.IWYZ] * FS[I][V.IRYZ]) * DTAU;
        FS[I][V.IWZZ] = FS[I][V.IWZZ] + (2.e0 * FS[I][V.IWXZ] * FS[I][V.IRXZ] + 2.e0 * FS[I][V.IWYZ] * FS[I][V.IRYZ]) * DTAU;
        FS[I][V.IWXY] = FS[I][V.IWXY] + ((FS[I][V.IWXX] - FS[I][V.IWYY]) * FS[I][V.IRXY] - FS[I][V.IWXZ] *
            FS[I][V.IRYZ] - FS[I][V.IWYZ] * FS[I][V.IRXZ]) * DTAU;
        FS[I][V.IWXZ] = FS[I][V.IWXZ] + ((FS[I][V.IWXX] - FS[I][V.IWZZ]) * FS[I][V.IRXZ] + FS[I][V.IWXY] *
            FS[I][V.IRYZ] - FS[I][V.IWYZ] * FS[I][V.IRXY]) * DTAU;
        FS[I][V.IWYZ] = FS[I][V.IWYZ] + ((FS[I][V.IWYY] - FS[I][V.IWZZ]) * FS[I][V.IRYZ] + FS[I][V.IWXY] *
            FS[I][V.IRXZ] - FS[I][V.IWXZ] * FS[I][V.IRXY]) * DTAU;
        //!change of Uij components at rotation
        FS[I][V.IUXX] = FS[I][V.IUXX] + (-2.e0 * FS[I][V.IUXY] * FS[I][V.IRXY] - 2.e0 * FS[I][V.IUXZ] * FS[I][V.IRXZ]) * DTAU;
        FS[I][V.IUYY] = FS[I][V.IUYY] + (2.e0 * FS[I][V.IUXY] * FS[I][V.IRXY] - 2.e0 * FS[I][V.IUYZ] * FS[I][V.IRYZ]) * DTAU;
        FS[I][V.IUZZ] = FS[I][V.IUZZ] + (2.e0 * FS[I][V.IUXZ] * FS[I][V.IRXZ] + 2.e0 * FS[I][V.IUYZ] * FS[I][V.IRYZ]) * DTAU;
        
        FS[I][V.IUXY] = FS[I][V.IUXY] + ((FS[I][V.IUXX] - FS[I][V.IUYY]) * FS[I][V.IRXY] - FS[I][V.IUXZ] *
            FS[I][V.IRYZ] - FS[I][V.IUYZ] * FS[I][V.IRXZ]) * DTAU;
       
        FS[I][V.IUXZ] = FS[I][V.IUXZ] + ((FS[I][V.IUXX] - FS[I][V.IUZZ]) * FS[I][V.IRXZ] + FS[I][V.IUXY] *
            FS[I][V.IRYZ] - FS[I][V.IUYZ] * FS[I][V.IRXY]) * DTAU;
        
        FS[I][V.IUYZ] = FS[I][V.IUYZ] + ((FS[I][V.IUYY] - FS[I][V.IUZZ]) * FS[I][V.IRYZ] + FS[I][V.IUXY] *
            FS[I][V.IRXZ] - FS[I][V.IUXZ] * FS[I][V.IRXY]) * DTAU;
        //!stress deviators
        FS[I][V.IULL] = FS[I][V.IUXX] + FS[I][V.IUYY] + FS[I][V.IUZZ];
        FS[I][V.ISXY] = 2.e0 * FS[I][V.IG] * (FS[I][V.IUXY] - FS[I][V.IWXY]);
        FS[I][V.ISXZ] = 2.e0 * FS[I][V.IG] * (FS[I][V.IUXZ] - FS[I][V.IWXZ]);
        FS[I][V.ISYZ] = 2.e0 * FS[I][V.IG] * (FS[I][V.IUYZ] - FS[I][V.IWYZ]);
        FS[I][V.ISXX] = 2.e0 * FS[I][V.IG] * (FS[I][V.IUXX] - C1D3 * FS[I][V.IULL] - FS[I][V.IWXX]);
        FS[I][V.ISYY] = 2.e0 * FS[I][V.IG] * (FS[I][V.IUYY] - C1D3 * FS[I][V.IULL] - FS[I][V.IWYY]);
        FS[I][V.ISZZ] = 2.e0 * FS[I][V.IG] * (FS[I][V.IUZZ] - C1D3 * FS[I][V.IULL] - FS[I][V.IWZZ]);
        FS[I][V.IMisSTR] = sqrt(pow((FS[I][V.ISXX] - FS[I][V.ISYY]), 2) + pow((FS[I][V.ISYY] - FS[I][V.ISZZ]), 2) + pow((FS[I][V.ISZZ] - FS[I][V.ISXX]), 2) +
            6.e0 * (pow(FS[I][V.ISXY], 2) + pow(FS[I][V.ISXZ], 2) + pow(FS[I][V.ISYZ], 2))) / sqrt(2.e0);
        MNO2 = 1.e0;
        FS[I][V.IDNS] = FS[I][V.IDNS] + FS[I][V.IDDNS] * DTAU * MNO2;
        FS[I][V.IU] = FS[I][V.IU] + FS[I][V.IDU] * DTAU * MNO2;    //!internal energy

        //!EOS BRASS
        

        if (PAR.metall == 1) SUBZ = 13.e0; //AL
        if (PAR.metall == 2) SUBZ = 29.e0; //CU
        DNS_EOS = FS[I][V.IDNS];
        T_EOS = FS[I][V.IT];
        E_EOS = FS[I][V.IU];
        //P_EOS = FS[I][V.IP];
        P_EOS = FS[I][V.IP];
        CS_EOS = FS[I][V.ICS];
        CV_EOS = FS[I][V.ICV];

        EOS_KH(DNS_EOS, E_EOS, T_EOS, P_EOS, CS_EOS, CV_EOS, SUBZ);
        

        FS[I][V.IP] = P_EOS;
        FS[I][V.IT] = T_EOS;
        FS[I][V.ICS] = CS_EOS;
        
        DWDQ_OUT = DWDQ;
        RR2_OUT = RR2;

        FS[I][V.IKK] = FS[I][V.IDNS] * pow(FS[I][V.ICS], 2);
        FS[I][V.IG] = 1.5e0 * FS[I][V.IKK] * (1.e0 - 2.e0 * FS[I][V.ICP]) / (1.e0 + FS[I][V.ICP]);
        

        for (int ALF = 0; ALF < V.NDIM; ALF++)
        {
            FS[I][V.IDX[ALF]] = FS[I][V.IDX[ALF]] + FS[I][V.IVREALV[ALF]] * DTAU;
            FS[I][V.IXX[ALF]] = FS[I][V.IXX[ALF]] + FS[I][V.IVREALV[ALF]] * DTAU;
            FS[I][V.IVREALV[ALF]] = FS[I][V.IVREALV[ALF]] + FS[I][V.IACS[ALF]] * DTAU;
        }


        /*std::cout
            << "   IG=" << FS[I][V.IG] <<std::endl
            << "   IKK=" << FS[I][V.IKK] << std::endl
            << "   IDNS=" << FS[I][V.IDNS] << std::endl
            << "   ICS=" << FS[I][V.ICS] << std::endl
            << "   IUXX=" << FS[I][V.IUXX] << std::endl
            << "   IULL=" << FS[I][V.IULL] << std::endl
            << "   IWXX=" << FS[I][V.IWXX] << std::endl
            << "   MSTR=" << FS[I][V.IMisSTR] << std::endl
            << std::endl;
        getchar();*/
        FS[I][V.ICV] = FS[I][V.ICV];
        //FS[I][V.ICV] = CV_EOS;
        

       

        //!Rescale
        for (int ALF = 0; ALF < 3; ALF++) //do ALF = 1, 3 //!num_axi
        {
            DX_RESCALE = -FS[I][V.IXX[ALF]] * V.vepsXYZ[ALF] * DTAU;
            FS[I][V.IDX[ALF]] = FS[I][V.IDX[ALF]] + DX_RESCALE;
            FS[I][V.IVV[ALF]] = FS[I][V.IVREALV[ALF]] + DX_RESCALE / DTAU;
            FS[I][V.IXX[ALF]] = FS[I][V.IXX[ALF]] + DX_RESCALE;
        }//enddo
        /*
        //!Boundary condition
        for (int ALF = 0; ALF < 3; ALF++) //do ALF = 1, 3 //!num_axi   !TRI - Axial
        {
            if (FS[I][V.IXX[ALF]] < CMINMAX[ALF][0]) //then
            {
                DX_SHIFT = CMINMAX[ALF][0] - FS[I][V.IXX[ALF]];
                FS[I][V.IDX[ALF]] = FS[I][V.IDX[ALF]] + DX_SHIFT;
                FS[I][V.IVV[ALF]] = FS[I][V.IVV[ALF]] + DX_SHIFT / DTAU;
                FS[I][V.IXX[ALF]] = CMINMAX[ALF][0];
            }//endif
            if (FS[I][V.IXX[ALF]] > CMINMAX[ALF][1]) //then
            {
                DX_SHIFT = CMINMAX[ALF][1] - FS[I][V.IXX[ALF]];
                FS[I][V.IDX[ALF]] = FS[I][V.IDX[ALF]] + DX_SHIFT;
                FS[I][V.IVV[ALF]] = FS[I][V.IVV[ALF]] + DX_SHIFT / DTAU;
                FS[I][V.IXX[ALF]] = CMINMAX[ALF][1];
            }//endif
        }//enddo
        */
    }//enddo

    TIME = TIME + DTAU;
    return TIME;
}//end subroutine SPHEP_MOVE
//!=====================================END_SPHEP_MOVE================================================ =