/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#ifndef _MPI_ENVIRONMENT_SETUP_H
#define _MPI_ENVIRONMENT_SETUP_H
#include "mpi_struct.h"
#include <iostream>
#include <vector>

// MPI_Setting
void Set_MPI(MPI_Setting *M, std::vector<int> const& np);



#endif