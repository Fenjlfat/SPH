#include "header.h"
#include "parametrs.h"
#include "particles.h"

//initialization vector of structure particle
void INIT(std::vector<particles> &particle)
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

                idParticle += 1;
            }
            
        }
        
    }
}