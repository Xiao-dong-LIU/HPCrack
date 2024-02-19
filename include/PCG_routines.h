/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _PCG_Routines_H
#define _PCG_Routines_H
#include "grid.h"
#include "mpi_struct.h"

double vv_product(const grid<double>& r, const grid<double>& z, const MPI_Setting & M);
double grid_norm(const grid<double>& r, const MPI_Setting & M);


#endif