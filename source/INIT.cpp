#include "header.h"
#include "parametrs.h"
#include "particles.h"

//initialization vector of structure particle
void INIT(double &E_EOS, double &P_EOS, double &T_EOS, double &CV_EOS, double &CS_EOS, double &DNS_EOS, std::vector<particles> &particle)
{
    parametrs par;
    
    int idParticle = 0; //number particle
    //double *xx;
    //double *yy;
    //double *zz;
    

    for (int i = 0; i < par.npaax; i++)
    {
        for (int j = 0; j < par.npaay; j++)
        {
            for (int k = 0; k < par.npaaz; k++)
            {
                double xx = i * par.dbp;
                double yy = j * par.dbp;
                double zz = k * par.dbp;
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
				particle[idParticle].IMAS = particle[idParticle].IDNS * pow(par.dbp, 3);
				particle[idParticle].IKK = particle[idParticle].IDNS * pow(particle[idParticle].ICS, 2);
				particle[idParticle].ICP = par.nuP;
				particle[idParticle].IG = 1.5e0 * particle[idParticle].IKK * 
                    (1.e0 - 2.e0 * particle[idParticle].ICP) / (1.e0 + particle[idParticle].ICP);
				particle[idParticle].IHS = 1.e0 * par.dbp;
				//std::cout << " IHS=" << FS[id_Particle][VAR.IHS] << std::endl;
				//getchar();

                idParticle += 1;
            }
            
        }
        
    }
}