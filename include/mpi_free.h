/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _MPI_FREE_H
#define _MPI_FREE_H
#include "mpi_struct.h"

// Frees the datatype
void mpi_data_free (MPI_Setting & M);

// Marks the communicator object for deallocation
void mpi_comm_free(MPI_Setting & M);


#endif