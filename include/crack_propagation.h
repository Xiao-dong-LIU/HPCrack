/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _CRACK_PROPAGATION_H
#define _CRACK_PROPAGATION_H

#include "mg.h"
#include "stack_and_level.h"
#include "mpi_struct.h"
#include "configuration.h"
//------------ the phase field method 
void crack_propagation(Stack *U, const mgdouble & bulK, const mgdouble & G, const mgdouble & gc, 
MPI_Setting & M, const Parameters & para, const int myid, const int nbprocs);

#endif
