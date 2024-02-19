/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#ifndef _d_sum_h
#define _d_sum_h

#include "grid.h"

// find the maximum variation of d over the entier domain.
double delta_d_max(const grid <double> & d, const grid <double> & d_old);

#endif