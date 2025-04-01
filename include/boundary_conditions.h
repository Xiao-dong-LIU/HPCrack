/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _BOUNDARY_CONDITIONS_H
#define _BOUNDARY_CONDITIONS_H
#include "mpi_struct.h"
#include "gdvec.h"
#include "stack_and_level.h"



// get the coordinate of the point in the global 
// @param origin: the origin of the domain
// @param h: the mesh size
// @param index: the index of the point in the local domain
// @param procIndex: the index of the processor in the direction 
// @param dim: the number of processors in the direction
// @param total: the number of points in the direction of the whole domain
// @param nbcolo: the number of points in the local domain
double get_coordinate(double origin, double h, int index, int procIndex, int dim, int total, int nbcolo);
// Dirichlet boundary conditions for the displacement field
// @param U: the stack structure
// @param u: the displacement field
// @param M: the MPI setting
// @param l: the level
// @param Ut: the prescribed displacement
void Dirichlet_BD(Stack *U, gdvecdouble & u, const MPI_Setting & M, const int l, const double Ut);

// Dirichlet force for the mechanical external force
// @param U: the stack structure
// @param f: the mechanical external force
// @param Fi: the internal force
// @param M: the MPI setting
// @param l: the level
void Dirichlet_force_u(Stack* U, gdvecdouble& f, const gdvecdouble & Fi, const MPI_Setting & M, 
const int l);
// Dirichlet residual for the residue 
// @param U: the stack structure
// @param r: the residue
// @param M: the MPI setting
// @param l: the level
void Dirichlet_residual_u(Stack* U, gdvecdouble & r, const MPI_Setting & M, const int l);


// Dirichlet boundary conditions for the phase field
// @param U: the stack structure
// @param d: the phase field
// @param M: the MPI setting
// @param l: the level
void Dirichlet_BD_d(Stack* U, griddouble & d, const MPI_Setting & M, const int l);void Dirichlet_force_d(Stack* U, griddouble& fd, const griddouble & Fi,  const MPI_Setting & M, const int l);
void Dirichlet_residual_d(Stack* U, griddouble & r, const MPI_Setting & M, const int l);


#endif // BOUNDARY_CONDITIONS_H
