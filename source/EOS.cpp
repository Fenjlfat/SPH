#include "header.h"
#include "parametrs.h"
#include "EOS.h"

void EOS_KHM::EOS_KHT(double RO, double &U, double &T, double &P, double &CS, double &CV, double SBL, int KE) 
{
    parametrs par;
    int NSB = 0;
    bool found = false;
    
    // Find substance index
    for (NSB = 0; NSB < NSUB0; NSB++) 
    {
        if (static_cast<int>(SBL) == static_cast<int>(ES[NSB][0])) 
        {
            found = true;
            break;
        }
    }
    
    if (!found) 
    {
        std::cerr << "EOS_1 ERROR!!! : No DATA for SB=" << SBL << std::endl;
        exit(1);
    }

    // Extract parameters
    double ro0 = ES[NSB][1];
    double Lam = ES[NSB][2];
    double Pk = ES[NSB][3];
    double Tk = ES[NSB][4];
    double rok = ES[NSB][5];
    double r = ES[NSB][6];
    double A = ES[NSB][7];
    double K = ES[NSB][8];
    double Btk = ES[NSB][9];
    double m = ES[NSB][10];
    double n = ES[NSB][11];
    double Gam = ES[NSB][12];
    double l = ES[NSB][13];
    double kl = ES[NSB][14];
    double al = ES[NSB][15];
    double C0 = ES[NSB][16];
    double C1 = ES[NSB][17];
    double C2 = ES[NSB][18];
    double a1 = ES[NSB][19];
    double b1 = ES[NSB][20];
    double a2 = ES[NSB][21];
    double b2 = ES[NSB][22];
    double miu = ES[NSB][23];
    double ZN = ES[NSB][24];
    double Tmlt = ES[NSB][25];
    
    double tau_mlt = Tmlt / Tk;
    double Rm32 = 1.5 * RG / miu;
    double Ek = Rm32 * Tk;
    
    // Calculations
    double fi = RO / rok;
    double tau = T / Tk;
    double del = fi / r;
    double zz = l * tau / std::pow(fi, kl);
    double del_m = std::pow(del, m * par.C1D3);
    double del_n = std::pow(del, n * par.C1D3);
    double bet = Btk / std::pow(fi,  par.C2D3);
    double cx = std::exp(al*(fi-1.0));
    double cex = fi * cx;
    double Z = (0.0 + cex) / (1.0 + cex); // Zc = 0, alf = 0 as per original
    double mn1 = A / (m - n);
    double Cz1 = 1.0 / (1.0 + zz);
    double Cz12 = Cz1 * Cz1;
    double btZ = bet * tau / Z;
    double TH_btZ = std::tanh(btZ);
    double dBdf = - par.C2D3 * bet / fi;
    double dZdf = (1.0 - 0.0) * (1.0 + al * fi) * cx / (1.0 + fi * cx); // Zc = 0
    
    // Internal energy
    double eps_u = A + mn1 * (n * del_m - m * del_n);
    double eps_t = (2.0 + zz) * Cz1 * tau;
    double LCH_btZ = btZ;
    if (btZ < 1.0e2) LCH_btZ = std::log(std::cosh(btZ));
    double eps_e = Z * Z * LCH_btZ / bet;
    double eps_i = 0.0; // alf = 0
    double eps = eps_u + eps_t + eps_e + eps_i;
    U = eps * Ek;
    
    // Pressure
    double pi_u = r * K * mn1 *  par.C1D3 * m * n * del * (del_m - del_n);
    double pi_t = (Gam + zz *  par.C1D3) / (1 + 0.5 * zz) * K * fi * eps_t;
    double pi_e =  par.C2D3 * K * fi * eps_e;
    double pi = pi_u + pi_t + pi_e;
    P = Pk * pi;
    
    // Heat capacity
    double teps_t = 1 + Cz12;
    double teps_e = Z * TH_btZ;
    double teps = teps_t + teps_e; // teps_i = 0
    CV = Rm32 * teps;
    
    // Speed of sound
    double feps_u = mn1 * n * m * (del_m - del_n) / (3.0 * r * del);
    double feps_t = kl * zz * tau * Cz12 / fi;
    double feps_e = 2.0 * eps_e * dZdf / Z - eps_e * dBdf / bet + TH_btZ * tau * (dBdf * Z / bet - dZdf);
    double feps = feps_u + feps_t + feps_e; // feps_i = 0
    double ftauS = (Pk / (Tk * rok * Rm32) * pi / (fi * fi) - feps) / teps;
    
    double tpi_t =  par.C2D3 * K * fi * (zz * zz + 2.0 * zz + 3.0 * Gam) * Cz12;
    double tpi_e =  par.C2D3 * K * fi * teps_e;
    double tpi = tpi_t + tpi_e;
    
    double fpi_u = K * mn1 *  par.C1D3 * m * n * del * (( par.C1D3 * m + 1) * del_m - ( par.C1D3 * n + 1) * del_n);
    double fpi_t = K * tau * (kl * zz * (2.0 * Gam -  par.C2D3) * Cz12 + (2.0 * Gam +  par.C2D3 * zz) * Cz1);
    double fpi_e =  par.C2D3 * K * (eps_e + fi * feps_e);
    double fpi = fpi_u + fpi_t + fpi_e;
    
    double CS2 = Pk / rok * (fpi + tpi * ftauS);
    CS = std::sqrt(std::abs(CS2));
    
    // Binodal boundary densities and saturated vapor pressure
    double RO_liq = 0.0;
    double RO_vap = 0.0;
    double P_vap = -1.0;
    
    if (tau <= 1.0 && KE == 1) 
    {
        double fi_liq = 1.0 + a1 * (1.0 - tau) + b1 * std::sqrt(1.0 - tau);
        double fi_vap = std::exp(-std::pow(fi_liq - 1.0, 2) * (a2 + b2 * std::pow((fi_liq - 1.0) / (fi_liq - r), 2)));
        double pi_vap = std::exp((1.0 - 1.0 / tau) * (C0 + C1 * tau + C2 * tau * tau));
        
        if (static_cast<int>(SBL) == 13) // AL
        { 
            fi_liq *= 1.2;
        }
        
        RO_liq = rok * fi_liq;
        RO_vap = rok * fi_vap;
        P_vap = Pk * pi_vap;
        
        if ((T > Tmlt) && (fi > fi_vap) && (fi < fi_liq)) 
        {
            double xx = (fi_liq / fi - 1.0) / (fi_liq / fi_vap - 1.0);
            if (P < 0.0) P = P_vap;
            
            double U_liq, Pliq, CS_liq, CV_liq;
            double U_vap, Pvap, CS_vap, CV_vap;
            
            EOS_KHT(RO_liq, U_liq, T, Pliq, CS_liq, CV_liq, SBL, 0);
            EOS_KHT(RO_vap, U_vap, T, Pvap, CS_vap, CV_vap, SBL, 0);
            
            U = U_liq * (1.0 - xx) + U_vap * xx;
            CS2 = std::pow(CS_liq, 2) * (1.0 - xx) + std::pow(CS_vap, 2) * xx;
            CS = std::sqrt(std::abs(CS2));
            CV = CV_liq * (1.0 - xx) + CV_vap * xx;
        }
    }
}

void EOS_KHM::EOS_KH(double RO, double U, double &T, double &P, double &CS, double &CV, double SUBS) 
{
    const double eps_DT = 1.0e-6;
    
    if (RO <= 0.0) return;
    
    if (std::abs(U) > 0.0) 
    {
        double DT = T;
        double DT1 = DT;
        double UT;
        
        while (std::abs(DT) > eps_DT) 
        {
            EOS_KHT(RO, UT, T, P, CS, CV, SUBS, 1);
            
            DT1 = DT;
            DT = (U - UT) / CV;
            
            if (DT * DT1 < 0.0 && std::abs(DT) >= std::abs(DT1)) 
            {
                DT = 0.5 * DT1;
            }
            
            if (DT < 0.0 && std::abs(DT) >= T) 
            {
                DT = -0.5 * T;
            }
            
            T += DT;
        }
    } 
    else 
    {
        EOS_KHT(RO, U, T, P, CS, CV, SUBS, 1);
    }
}

void EOS_KHM::SLRelax(double &RO, double &U, double T,  double &P, double &CS, double &CV, double SBL) 
{
    double DRO;
    U = 0.0;
    EOS_KHT(RO, U, T, P, CS, CV, SBL, 1);
    
    while (std::abs(P) > 1.0e-3) 
    {
        DRO = -P / std::pow(CS, 2);
        RO += DRO;
        U = 0.0;
        EOS_KHT(RO, U, T, P, CS, CV, SBL, 1);
    }
}