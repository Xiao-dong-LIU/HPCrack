/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "sigma.h"



void Sigma_g(grid2d & sigma, const grid2d & epsilon, const double K_g, const double G_g, 
const double tr_epsilon)
{
    for (int c=0; c<3; c++)
        for (int m=0; m<3; m++)
        {
            if(c==m)
            {
                sigma(c,m)=K_g*tr_epsilon+2.0*G_g*(epsilon(c,m)-1.0/3.0*tr_epsilon);
            }
            else
            sigma(c,m)=2.0*G_g*epsilon(c,m);
        }
}