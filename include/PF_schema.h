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
#include "parameters_config.h"
//------------ the phase field method 
void phm(Stack *U, const mgdouble & bulK, const mgdouble & G, const mgdouble & gc, 
MPI_Setting & M, const Parameters & para, const int myid, const int nbprocs);

#endif
