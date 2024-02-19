/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _mgcycle_H
#define _mgcycle_H
#include "mg.h"
#include "mgvec.h"
#include "mpi_struct.h"
#include "structure_df.h"

///------------MULTIGRID DRIVING ROUTINES------------///

void cycle_u(Stack *U, mgvecdouble &u, mgvecdouble & ru, mgvecdouble & fu, mgvecdouble & uold, const mgshort & MP_choice_d, 
const mgdouble & bulK,  const mgdouble & G,  const mgdouble & d,   const mgdouble & K_element, 
const mgdouble & G_element, int l, const MG & mgp_u,  MPI_Setting & M, const int myid, 
const int nbprocs, const double small_k, const double Ut, const int t, double & f_norm);




#endif // _mgcycle_H
