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
#include "config.h"
using namespace std;

/// ------------- initialize values in datastructure
void initialize(Stack *U, const Parameters & para, const MPI_Setting &M, const int maxlevel);

#endif // _initialisation_H
