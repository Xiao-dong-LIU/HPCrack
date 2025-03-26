/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _MATERIAL_COARSE_ELEMENT_H
#define _MATERIAL_COARSE_ELEMENT_H

#include "mg.h"
#include "mpi_struct.h"

// ------- compute voigt average on coarse gird
void coarse_voigt(griddouble & Hc, const griddouble & Hf, const MPI_Setting & M);
// ------- compute reuss average on coarse gird
void coarse_reuss(griddouble& Hc, const griddouble& Hf, const MPI_Setting & M);
//------ compute voigt-reuss type element material property on coarse girds
void VR_all_level(mgdouble & H, MPI_Setting & M);
//------ compute voigt type element material property on coarse girds
void Voigt_all_level(mgdouble & H, MPI_Setting & M);

#endif