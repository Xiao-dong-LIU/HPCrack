/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#ifndef _MULTIGRID_INTERGRID_H
#define _MULTIGRID_INTERGRID_H
#include "grid2d.h"
#include "grid.h"
#include "gdvec.h"
#include "mpi_struct.h"

//------------- Injection of solution on level l to level l-1
void solution_injection(grid <double> & uc, const grid <double> & uf, const MPI_Setting & M);

//-------------- injection for d
void coarsen_d(grid<double> &dc, grid<double> &doldc, const grid<double> &df, MPI_Setting &M);
//------------- Injection for displacement
void coarsen_u(Stack *U, gdvecdouble &uc, gdvecdouble &ucold, const gdvecdouble &uf,
			   const griddouble &bulK, MPI_Setting &M, const int lc, double Ut);
// standard restriction operator
void coarsen_f_standard_routine(grid<double> &fc, const grid<double> &rf,
								const grid<double> &Fic, const MPI_Setting &M);
// restriction for d
void coarsen_f_d(Stack *U, grid<double> &fc, const grid<double> &dc, grid<double> & rf,
const grid<double> &gcc, const grid<double> &Hc, MPI_Setting &M, const double lc, const int lf);
/// restriction for u
void coarsen_fu(Stack *U, gdvecdouble &fc, const gdvecdouble &uc, const gridshort &MP_choice_dc, 
gdvecdouble &rf, const griddouble &bulK, const griddouble &G, const griddouble &d, 
const griddouble &K_element, const griddouble &G_element, MPI_Setting &M, const int lf, 
double small_k);
//------------- standard prologation operator
void refine_standard_routine(grid<double> &uf, const grid<double> &uc,
const grid<double> &ucold, const MPI_Setting &M);
//------------- prologate d
void refine_d(grid<double> &df, const grid<double> &dc, const grid<double> doldc, MPI_Setting &M);
//-------------- prolongate u
void refine_u(Stack *U, gdvecdouble &uf, const gdvecdouble &uc, const gdvecdouble &ucold, 
MPI_Setting &M, const double Ut, const double lfine);

#endif // MULTIGRID_INTERGRID_H
