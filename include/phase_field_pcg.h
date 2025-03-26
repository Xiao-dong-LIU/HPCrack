/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _PHASE_FIELD_PCG_H
#define _PHASE_FIELD_PCG_H
#include "grid2d.h"
#include "grid.h"
#include "mpi_struct.h"


void Md_element(grid <double> & Md, const grid <double> & H, const grid <double> & gc, 
const grid <double> & DN, grid2d & N, const double lc, const double wg, 
const int i, const int j, const int k);

void Precondition_M(grid <double> & Md, const grid <double> & H, const grid <double> & gc, 
const grid <double> & DN, grid2d& N, const double lc, const double wg);

void Fld_element(grid <double> & Fi, const grid <double> & d, const grid <double> & H, 
const grid <double> & gc, const grid <double> & DN, const grid2d & N, const double lc, 
const double wg, const int i, const int j, const int k);

void Fld(grid <double> & Fi, const grid <double> & d, const grid <double> & H, 
const grid <double> & gc, const grid <double> & DN, const grid2d & N, const double lc, 
const double wg);

double initial_residual(Stack *U, grid <double> & R,  grid <double> & d, 
const grid <double> & fd,  const grid <double> & H, const grid <double> & gc, 
const int l, const double lc, const MPI_Setting & M);

void pcg_d(Stack *U, grid<double> & d, grid<double> & fd, grid<double> & R, 
const grid<double> & H, const grid<double> & gc, MPI_Setting & M, const int l, 
const double lc,  const int NB_relax, const int t, double & fd_norm);

#endif