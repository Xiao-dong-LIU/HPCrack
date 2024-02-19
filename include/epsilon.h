/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _EPSILON_h
#define _EPSILON_H
#include "grid2d.h"

#include "gdvec.h"

//---------------- Compute epsilon at one gauss point
void epsilon_routine(grid2d& epsilon, const gdvecdouble & u, const griddouble & DN,  
const int i, const int j, const int k, const int g) ;
//---------------- Compute the trace of epsilon at one gauss point
double trace_epsilon(const grid2d& epsilon);
//---------------- Double dot product of the epsilon_dev at gauss integration point
double epsilon_dev(const grid2d& epsilon);

#endif // _EPSILON_H
