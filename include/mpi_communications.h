/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#ifndef _MPI_COMMUNICATIONS_H
#define _MPI_COMMUNICATIONS_H
#include "mpi_struct.h"
#include "grid.h"

///------------ MPI Routines ------------///


void derived_datatypes(const vector <size_t> & table_dim, MPI_Setting *M);


// mpi communication 
// nboflayers : the nomber of layers needed to be communicated >=1
template <class T>
void SendRecv (grid<T> &u, const MPI_Setting & M)
{
		int x=u.axis()[0];
		int y=u.axis()[1];
		int z=u.axis()[2];
	//	derived_datatypes<T>(x,y,z,&M);
		MPI_Status status;
		/// Communication in X direction
		MPI_Ssend(&u(x-2,0,0),1,M.type_x,M.idfront,0,M.X_row_comm);
		MPI_Recv(&u(0,0,0),1,M.type_x,M.idbehind,0,M.X_row_comm,&status);
		MPI_Ssend(&u(1,0,0),1,M.type_x,M.idbehind,0,M.X_row_comm);
		MPI_Recv(&u(x-1,0,0),1,M.type_x,M.idfront,0,M.X_row_comm,&status);
		/// Communication in Y direction
		MPI_Ssend(&u(0,y-2,0),1,M.type_y,M.idright,0,M.Y_row_comm);
		MPI_Recv(&u(0,0,0),1,M.type_y,M.idleft,0,M.Y_row_comm,&status);
		MPI_Ssend(&u(0,1,0),1,M.type_y,M.idleft,0,M.Y_row_comm);
		MPI_Recv(&u(0,y-1,0),1,M.type_y,M.idright,0,M.Y_row_comm,&status);
		/// Communication in Z direction
		MPI_Ssend(&u(0,0,z-2),1,M.type_z,M.idup,0,M.Z_row_comm);
		MPI_Recv(&u(0,0,0),1,M.type_z,M.iddown,0,M.Z_row_comm,&status);
		MPI_Ssend(&u(0,0,1),1,M.type_z,M.iddown,0,M.Z_row_comm);
		MPI_Recv(&u(0,0,z-1),1,M.type_z,M.idup,0,M.Z_row_comm,&status);
};

#endif