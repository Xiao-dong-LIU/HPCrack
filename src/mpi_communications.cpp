/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "mpi_communications.h"



void derived_datatypes(const vector <size_t> & table_dim, MPI_Setting *M)
{
  int x = table_dim[0];
  int y = table_dim[1];
  int z = table_dim[2];
  MPI_Type_contiguous(x*y,MPI_DOUBLE,&M->type_z);
  MPI_Type_commit(&M->type_z);
  MPI_Type_vector(z,x,(x*y),MPI_DOUBLE,&M->type_y);
  MPI_Type_commit(&M->type_y);
  vector <int> block_length(y*z,1);
  // int block_length[(y*z)];
  vector <int> displacement(y*z);
  // int displacement[(y*z)];
	for (int j =0;j<y;j++)
	  for (int k =0;k<z;k++)
	  {
		  displacement[j+k*y] = x*j + k*(x*y);
	  }
	MPI_Type_indexed(z*y,block_length.data(),displacement.data(),MPI_DOUBLE,&M->type_x);
  MPI_Type_commit(&M->type_x);
}
