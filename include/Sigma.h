/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _Sigma_h
#define _Sigma_h

#include "grid2d.h"

void Sigma_g(grid2d & sigma, const grid2d & epsilon, const double K_g, const double G_g, 
const double tr_epsilon);

#endif // _Sigma_h