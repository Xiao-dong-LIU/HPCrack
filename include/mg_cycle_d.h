/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _mg_cycle_d_h
#define _mg_cycle_d_h

#include "mg.h"
#include "mpi_struct.h"
#include "structure_df.h"

void cycle_d(Stack *U, mg<double> &d, mg<double> &fd, mg<double> &dold, const mg<double> &H, 
const mg<double> &gc, mg<double> &r, int l, const double lc, const MG & mgp_d, 
MPI_Setting & M, const int myid, const int nbprocs, const int t, double & fd_norm);



#endif //