/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#ifndef _MPI_STRUCT__H
#define _MPI_STRUCT__H
#include "mpi.h"

typedef struct
{
	int X_id;			   ///Processor is in direction X
	int Y_id;			   ///Processor is in direction Y
	int Z_id;			   ///Processor is in direction Z
	int dims[3];           ///Number of processors in each direction
	int coordinates[3];    ///Processor coordinates with Topology
    int idup;              /// processor id in direction Z
	int iddown;            /// processor id in direction Z
	int idright;           /// processor id in direction Y
	int idleft;            /// processor id in direction Y
	int idfront;           /// processor id in direction X
	int idbehind;          /// processor id in direction X
	MPI_Datatype type_x;   /// Derived datatypes for comminucation in direction x
	MPI_Datatype type_y;   /// Derived datatypes for comminucation in direction y
	MPI_Datatype type_z;   /// Derived datatypes for comminucation in direction z
	MPI_Comm X_row_comm;   ///Communicator for X firection
	MPI_Comm Y_row_comm;   ///Communicator for Y firection
	MPI_Comm Z_row_comm;   ///Communicator for Z firection
} MPI_Setting;

#endif