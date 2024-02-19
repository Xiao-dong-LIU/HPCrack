/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _OUTPUT_RESULT_H
#define _OUTPUT_RESULT_H

#include <string>
#include "structure_df.h"
#include "mpi_struct.h"
#include "gdvec.h"




void open_file(const char *filename, int l, const int myid, int t);
void close_file(void);
void write_header(Stack *U, const MPI_Setting & M, const int l);
void write_pointdata_head(const string & output_option, const char * field_name);
void write_scalar(Stack *U, const griddouble& d, const MPI_Setting & M, const int l);
void write_vector(Stack *U, const gdvecdouble & u, const MPI_Setting & M, 
const int l, double u_proportion);
void write_pointdata_end(void);
void write_celltdata(void);
void write_celltdata_H(Stack *U, const griddouble &H, const MPI_Setting & M, const int l, 
const char *field_name);
void write_end(void);
void write_pvtk(Stack *U, const char *filename, const MPI_Setting & M, const string output_option, 
const char *field_name, const int l, const int myid, const int nbprocs, const int t);
//-------- output scalar field: d or E
void write_total_de(Stack *U, const griddouble & de, const MPI_Setting & M, 
const string & output_option, const char *filename, const char *field_name, 
const int l, const int myid, const int nbprocs, const int t);
//-------- output scalar field: H
void write_total_H(Stack *U, const griddouble &H, const MPI_Setting & M, const string & output_option, 
const char *filename, const char *field_name, const int l, const int myid, const int nbprocs, 
const int t);
//-------- output vector field: u
void write_total_u(Stack *U, const gdvecdouble &u, const MPI_Setting & M, const int l, 
const int myid, const int nbprocs, const char *filename, const char *field_name, const int t,
const double u_proportion);
void write_total_Eu(Stack *U, const gdvecdouble &u, const griddouble &de, const int l, 
const MPI_Setting & M, const int myid, const int nbprocs, const char *filename,
const char * field_name, const int t,const double u_proportion);


#endif // _OUTPUT_RESULT_H
