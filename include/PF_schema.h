/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _PF_schema_h
#define _PF_schema_h

#include "mg.h"
#include "structure_df.h"
#include "mpi_struct.h"
//------------ the phase field method 
void phm(Stack *U, const mgdouble & bulK, const mgdouble & G, const mgdouble & gc, 
MPI_Setting & M, const MG & mgp_u, const MG & mgp_d, const double lc,
const int myid, const int nbprocs, double voxel_size);

#endif
