/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _REALDOMAIN_H
#define _REALDOMAIN_H
#include <iostream>
#include <vector>
#include "mpi_struct.h"
using namespace std;

template <class T>
vector <T> realdomain (const vector <T> & inD, const MPI_Setting & M)
{
    vector <T> outD(inD.size()*2,0);
    for (int d=0; d<inD.size(); d++)
    {
        outD[2*d+1] = inD[d];
        if(M.dims[d]>1)
        {
            if (M.coordinates[d]==0) 
            {
                outD[2*d+1] = inD[d]-1;
            }
            else if (M.coordinates[d]==M.dims[d]-1) 
            {
                outD[2*d] = 1;
            }
            else 
            {
                outD[2*d] = 1;
                outD[2*d+1] = inD[d]-1;
            }
        }
    }
    return outD;
}

#endif
