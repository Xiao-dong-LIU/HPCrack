/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _MATERIAL_GAUSS_POINTS_H
#define _MATERIAL_GAUSS_POINTS_H
#include "grid2d.h"
#include "grid.h"

//        point from nodal material property
double gauss_point_MP(const grid <double> & MP, const grid2d & N, const int i, const int j, 
const int k, const int g);
//------------- compute K_d or G_d at each gauss integration point
double K_G_d_g(const griddouble &d, const griddouble &MP, const grid2d & N, const int i, const int j, const int k, 
const int g, const double small_k);
//-------------function to obtain (1-k)*g(d)+k
double g_d_node(const griddouble &d, const int i, const int j, const int k, const double small_k);
//------------- compute K_d and G_d at each node
double MP_d_node(const griddouble &d, const griddouble &MP, 
const int i, const int j, const int k, const double small_k);
//------------- compute K_d and G_d at nodes of the fineest grid
void MP_d(griddouble & MP_d, const griddouble &d, const griddouble & MP, const double small_k);
#endif // MATERIAL_GAUSS_POINTS_H
