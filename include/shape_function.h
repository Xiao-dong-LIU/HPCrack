/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _SHAPE_FUNCTION_H
#define _SHAPE_FUNCTION_H
#include "grid2d.h"
#include "grid.h"

///    Shape function  and Derivative of shape function
void Dshape_function(grid2d& N, grid <double>& DN, const double hx, const double hy, const double hz);
#endif // SHAPE_FUNCTION_H