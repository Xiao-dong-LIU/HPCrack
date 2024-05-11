/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _d_total_h
#define _d_total_h
#include "mpi_struct.h"
#include "structure_df.h"
#include "mg.h"
//------------- phase field routine for computing d
void d_entire(Stack *U, mg<double> &d, mg<double> &fd, const mg<double> & H, const mg<double> &gc, 
MPI_Setting & M, double lc, const MG & mgp_d, int myid, int nbprocs,int t, double & fd_norm);

#endif