/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _write_H_h
#define _write_H_h

#include "grid.h"
#include "mpi_struct.h"
#include "structure_df.h"
// --------- write solution to binary data to avoid previous time step computing
void write_vector_3D(griddouble & H, const string & filename, const int myid, const int t);
// --------- read solution from binary data stocked by the previous computation
void read_vector_3D(griddouble & H, const string & filename, const int myid, const int t);

#endif