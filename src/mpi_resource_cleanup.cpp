/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "mpi_resource_cleanup.h"




void mpi_data_free (MPI_Setting & M)
{
    MPI_Type_free(&M.type_x);
    MPI_Type_free(&M.type_y);
    MPI_Type_free(&M.type_z);
}


void mpi_comm_free(MPI_Setting & M)
{
    MPI_Comm_free(&M.X_row_comm);
    MPI_Comm_free(&M.Y_row_comm);
    MPI_Comm_free(&M.Z_row_comm);
}