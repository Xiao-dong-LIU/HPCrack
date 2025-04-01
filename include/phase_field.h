/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _PHASE_FIELD_H
#define _PHASE_FIELD_H

#include "mg.h"
#include "stack_and_level.h"
#include "mpi_struct.h"
//------------ the phase field method 
void phase_field(Stack *U, mg<double> &d, mg<double> &fd, const mg<double> & H, const mg<double> &gc, 
MPI_Setting & M, double lc, const MG & mgp_d, int myid, int nbprocs,int t, double & fd_norm);

#endif
