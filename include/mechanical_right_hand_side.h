/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _MECHANICAL_RIGHT_HAND_SIDE_H
#define _MECHANICAL_RIGHT_HAND_SIDE_H

#include "gdvec.h"
#include "mpi_struct.h"
#include "stack_and_level.h"

// ----- nodal externe force 
void externe_force_nodal(Stack *U, gdvecdouble & fu, gdvecdouble & u, const gridshort & MP_choice_d, 
const griddouble & bulK, const griddouble & G, const griddouble & d,  const griddouble & K_element, 
const griddouble & G_element, const MPI_Setting & M, const int l, const double Ut, const double small_k);
// ----- externe force on traction direction
double externe_force_sum(griddouble & f, const MPI_Setting & M);

#endif
