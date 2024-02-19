/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _initialisation_H
#define _initialisation_H

#include "structure_df.h"
#include "mpi_struct.h"
#include <iostream>
#include <vector>
using namespace std;

/// ------------- initialize values in datastructure
void initialize(Stack *U, const inputdomain *IDM, const int maxlevel, vector <int> const& element_nb, 
vector <double> const& X_start, const MPI_Setting &M);
/*
///------------initialisation of u,f,d,H------------///
void init_ufdH(Stack *U, multigrid3d & u, multigrid3d & v, multigrid3d & w,
			   multigrid3d & fu, multigrid3d & fv, multigrid3d & fw,
               multigrid3d & d, multigrid3d & fd, multigrid3d & H,
               int l,MPI_Setting & M);
*/
#endif // _initialisation_H
