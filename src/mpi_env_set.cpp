/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "mpi_env_set.h"



///------------ MPI Routines ------------///

void Set_MPI(MPI_Setting *M, std::vector<int> const& np){
    /// number of processors in each direction
    int dims[3];
	for (int d=0; d<3; d++)
	{
		M->dims[d]=np[d];
		dims[d]=np[d];
	}
	/// control the reorder of rank number
	int reorder = 1;
	int wrap_around[3]={1,1,1};
	/// Communicator for the MPI Topology
	MPI_Comm grid_comm;
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, wrap_around,reorder,&grid_comm);
    /// Get the coordinates of each processor
    int coordinates[3];
    int my_grid_id,grid_id;
    MPI_Comm_rank(grid_comm, &my_grid_id);
    MPI_Cart_coords(grid_comm, my_grid_id, 3, coordinates);
    MPI_Cart_rank(grid_comm, coordinates, &grid_id);
    for(int i = 0 ; i <3;i++)
    {
		M->coordinates[i]=coordinates[i];
	}
    /// set row coommunicator paralleled with direction Z
    int free_coords[3] = {0,0,1};
    MPI_Cart_sub(grid_comm, free_coords, &M->Z_row_comm);
    MPI_Comm_rank(M->Z_row_comm, &M->Z_id);
    /// set row coommunicator paralleled with direction y
    free_coords[0] = 0,free_coords[1] = 1,free_coords[2] = 0;
    MPI_Cart_sub(grid_comm, free_coords, &M->Y_row_comm);
    MPI_Comm_rank(M->Y_row_comm, &M->Y_id);
    /// set row coommunicator paralleled with direction x
    free_coords[0] = 1,free_coords[1] = 0,free_coords[2] = 0;
    MPI_Cart_sub(grid_comm, free_coords, &M->X_row_comm);
    MPI_Comm_rank(M->X_row_comm, &M->X_id);
    /// Setting processor to communication
    M->idup = M->Z_id + 1;
	M->iddown = M->Z_id - 1;
	M->idright = M->Y_id + 1;
	M->idleft = M->Y_id - 1;
	M->idfront = M->X_id + 1;
	M->idbehind = M->X_id - 1;
	if(M->X_id==0) M->idbehind = MPI_PROC_NULL;
	if(M->X_id==(M->dims[0]-1)) M->idfront = MPI_PROC_NULL;
	if(M->Y_id==0) M->idleft = MPI_PROC_NULL;
	if(M->Y_id==(M->dims[1]-1)) M->idright = MPI_PROC_NULL;
	if(M->Z_id==0) M->iddown = MPI_PROC_NULL;
	if(M->Z_id==(M->dims[2]-1)) M->idup = MPI_PROC_NULL;
}
