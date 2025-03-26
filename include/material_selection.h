/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _MATERIAL_SELECTION_H
#define _MATERIAL_SELECTION_H

#include "grid2d.h"
#include "mgvec.h"
#include "mg.h"

void MP_choice_element(gridshort & MP_choice, const gdvecdouble & u, const griddouble & bulK, 
const griddouble & G, const griddouble & d,  const griddouble & DN, const grid2d & N, 
int l, const int i, const int j, const int k);
void MP_choice(Stack *U, mgshort & MP_choice_d, const mgvecdouble & u, const mgdouble & bulK, 
const mgdouble & G, const mgdouble & d);

// void account_nb(Stack *U,multiinput3d& MP_choice_d,MPI_Setting & M);

void MP_element(griddouble & K_element, griddouble & G_element, const gdvecdouble & u, 
const griddouble & bulK, const griddouble & G, const griddouble & d, const gridshort & MP_choice_d, 
const griddouble & DN, const grid2d & N, int i, int j, int k, double small_k);

void MP_all_element(Stack*U, const gdvecdouble & u,  const griddouble & bulK, 
const griddouble & G, const griddouble & d, const gridshort & MP_choice_d, griddouble & K_element, 
griddouble & G_element, const int l, const double small_k);
// void account_nb_zeros(Stack *U, multigrid3d& MP, MPI_Setting & M);
#endif
