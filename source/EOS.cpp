void EOS(double &E_EOS, double &P_EOS, double &T_EOS, double &CV_EOS, double &CS_EOS, double &DNS_EOS)
{
    E_EOS = 0.00;
    P_EOS = 0.00;
    T_EOS = 0.00;
    CV_EOS = 0.00;
    CS_EOS = 0.00;
    DNS_EOS = 0.00;
}
#include <iostream>
#include <cmath>
#include <vector>

const double C1D3 = 1.0 / 3.0;
const double C2D3 = 2.0 / 3.0;
const double RG = 8.314E3; // J/(kmol*K)
const int NSUB0 = 13;
const int NPAR0 = 26;

class EOS_KHM 
{
private:
    double SUBZ;
    
    // Substance parameters storage
    std::vector<std::vector<double>> ES = {
        // LI
        {3.0,  0.547e3, 22.950e6,  689.0e5, 3223.0, 0.105e3, 5.210, 3.694,  8.830, 0.575, 4.0, 2.0, 1.013, 2.579,
        0.879,  2.628, 6.403, -3.3760,  7.2890,  1.4040, 2.751, 0.904, 2.580e-2,    6.939,  3.0,  453.54},
        // NA
        {11.0,  1.014e3,  4.728e6,  275.0e5, 2573.0, 0.206e3, 4.922, 3.387, 10.455, 0.650, 5.0, 2.0, 1.104,  4.288,
        1.247, 2.683, 4.784, -0.4325, -0.0837,  1.2163, 2.639, 0.447, 6.040e-2,   22.990, 11.0,  370.90},
        // K
        {19.0,  0.903e3,  2.311e6,  152.0e5, 2223.0, 0.194e3, 4.655, 3.259,  9.050, 0.832, 5.0, 2.0, 1.066,  4.590,
        1.006, 2.153, 4.268,  0.5510, -0.8674,  0.1136, 3.420, 0.383, 3.670e-2,   39.098, 19.0,  336.50},
        // RB
        {37.0,  1.613e3,  0.972e6,  159.0e5, 2093.0, 0.346e3, 4.662, 3.183,  6.646, 0.897, 5.0, 2.0, 1.204,  5.454,
        0.747, 1.729, 4.871, -2.8720,  3.4480,  0.3755, 3.195, 0.382, 3.350e-2,   85.468, 37.0,  312.49},
        // CS
        {55.0,  2.000e3,  0.600e6,  144.0e5, 2057.0, 0.428e3, 4.673, 3.108,  5.737, 1.027, 5.0, 2.0, 1.028,  4.482,
        0.619, 1.383, 4.509, -1.8630,  2.5090, -0.8070, 4.314, 0.228, 2.040e-2,  132.905, 55.0,  301.50},
        // CU
        {29.0,  9.020e3,  5.296e6, 7460.0e5, 8390.0, 2.390e3, 3.774, 3.216,  5.276, 1.810, 6.0, 4.0, 1.917,  4.871,
        1.479, 1.940, 4.682, -1.6400,  3.9360, -0.2290, 2.934, 1.020, 3.260e-2,   63.546, 29.0, 1356.00},
        // ZN
        {30.0,  7.320e3,  1.990e6, 2630.0e5, 3190.0, 2.290e3, 3.197, 3.269,  5.298, 0.409, 7.0, 5.0, 2.370,  3.625,
        2.061, 6.385, 4.726, -1.5660,  3.6520,  0.5750, 1.650, 2.841, 1.090e-2,   65.380, 30.0,  692.50},
        // CD
        {48.0,  8.830e3,  0.995e6, 1600.0e5, 2790.0, 2.740e3, 3.220, 3.214,  5.301, 0.455, 8.0, 6.0, 2.179,  4.055,
        2.905, 6.543, 4.634, -1.5940,  2.5880, -0.3900, 2.641, 1.375, 1.210e-2,  112.410, 48.0,  593.90},
        // HG
        {80.0, 14.400e3,  0.322e6, 1510.0e5, 1763.0, 4.200e3, 3.429, 2.938,  3.049, 0.318, 9.0, 6.0, 2.944,  7.370,
        2.515, 5.185, 4.181, -0.6970,  0.9900,  0.8995, 1.636, 1.025, 0.025e-2,  200.590, 80.0,  234.10},
        // AL
        {13.0,  2.710e3,  9.210e6, 4470.0e5, 8000.0, 0.640e3, 4.230, 2.490,  5.294, 1.522, 6.0, 4.0, 2.121, 20.820,
        2.431, 2.361, 4.384, -0.6920, -0.6930, -0.6920, 2.375, 0.515, 7.300e-2,   26.981, 13.0,  933.24},
        // PB
        {82.0, 11.310e3,  0.945e6, 1840.0e5, 4980.0, 3.250e3, 3.480, 3.153,  5.294, 1.374, 7.0, 5.0, 2.693,  7.253,
        2.048, 2.946, 4.679, -1.1030,  5.6810,  0.2245, 2.163, 2.656, 0.709e-2,  207.200, 82.0,  600.44},
        // BI
        {83.0, 10.430e3,  1.058e6, 1260.0e5, 4200.0, 2.660e3, 3.920, 4.222,  5.291, 1.437, 6.0, 4.0, 2.048,  2.800,
        1.604, 1.750, 4.993, -1.7430,  3.6230, -1.0220, 3.822, 0.479, 1.680e-2,  208.980, 83.0,  490.40},
        // FE
        {26.0,  7.820e3,  7.427e6, 8250.0e5, 9600.0, 2.030e3, 3.850, 3.465,  5.275, 1.905, 6.0, 4.0, 1.805,  3.810,
        1.547, 1.759, 4.832, -3.6390,  5.0610, -0.5360, 3.320, 0.755, 2.480e-2,   55.847, 26.0, 1811.00}
    };

public:
    void EOS_KHT(double RO, double& U, double& T, double& P, double& CS, double& CV, double SBL, int KE) 
    {
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
        double del_m = std::pow(del, m*C1D3);
        double del_n = std::pow(del, n*C1D3);
        double bet = Btk / std::pow(fi, C2D3);
        double cx = std::exp(al*(fi-1.0));
        double cex = fi * cx;
        double Z = (0.0 + cex) / (1.0 + cex); // Zc = 0, alf = 0 as per original
        double mn1 = A / (m - n);
        double Cz1 = 1.0 / (1.0 + zz);
        double Cz12 = Cz1 * Cz1;
        double btZ = bet * tau / Z;
        double TH_btZ = std::tanh(btZ);
        double dBdf = -C2D3 * bet / fi;
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
        double pi_u = r * K * mn1 * C1D3 * m * n * del * (del_m - del_n);
        double pi_t = (Gam + zz * C1D3) / (1 + 0.5 * zz) * K * fi * eps_t;
        double pi_e = C2D3 * K * fi * eps_e;
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
        
        double tpi_t = C2D3 * K * fi * (zz * zz + 2.0 * zz + 3.0 * Gam) * Cz12;
        double tpi_e = C2D3 * K * fi * teps_e;
        double tpi = tpi_t + tpi_e;
        
        double fpi_u = K * mn1 * C1D3 * m * n * del * ((C1D3 * m + 1) * del_m - (C1D3 * n + 1) * del_n);
        double fpi_t = K * tau * (kl * zz * (2.0 * Gam - C2D3) * Cz12 + (2.0 * Gam + C2D3 * zz) * Cz1);
        double fpi_e = C2D3 * K * (eps_e + fi * feps_e);
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
    
    void EOS_KH(double RO, double U, double& T, double& P, double& CS, double& CV, double SUBS) 
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
    
    void SLRelax(double SBL, double& RO, double T, double& U, double& P, double& CS, double& CV) 
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
};