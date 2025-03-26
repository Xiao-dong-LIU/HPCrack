/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _MECHANICAL_DISPLACEMENT_H
#define _MECHANICAL_DISPLACEMENT_H

#include "mgvec.h"
#include "mg.h"
#include "mpi_struct.h"

void u_entire(Stack *U, mgvecdouble &u, mgvecdouble &fu, const mgdouble &bulK, const mgdouble &G, 
const mgdouble &d, MPI_Setting & M, const MG & mgp_u, const int myid, const int nbprocs, 
const double small_k, const double Ut, const int t, double & f_norm);

#endif // MECHANICAL_DISPLACEMENT_H