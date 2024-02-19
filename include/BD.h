/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _BD_h
#define _BD_H
#include "grid2d.h"
#include "mpi_struct.h"
#include "gdvec.h"


void Dirichlet_BD(Stack *U, gdvecdouble & u, const MPI_Setting & M, const int l, const double Ut);
void Dirichlet_force_u(Stack* U, gdvecdouble& f, const gdvecdouble & Fi, const MPI_Setting & M, 
const int l);
void Dirichlet_residual_u(Stack* U, gdvecdouble & r, const MPI_Setting & M, const int l);

void Dirichlet_BD_d(Stack* U, griddouble & d, const MPI_Setting & M, const int l);
void Dirichlet_force_d(Stack* U, griddouble& fd, const griddouble & Fi,  const MPI_Setting & M, const int l);
void Dirichlet_residual_d(Stack* U, griddouble & r, const MPI_Setting & M, const int l);


#endif // _BD_H
