/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef strain_history_H
#define strain_history_H
#include "grid2d.h"
#include "gdvec.h"
#include "mpi_struct.h"


//---------------- Compute phi_c at Gauss point
double Phi_c(const double gc_g, const double lc) ;
//---------------- initial crack 
void initial_crack(Stack *U, griddouble & H, const MPI_Setting & M, const double H_crack, const int l);
//---------------- Compute elastic energy at Gauss point
double W(const grid2d& sigma, const grid2d& epsilon) ;
//---------------- Compute the elastic displacement limit
double U_e(Stack* U, const gdvecdouble& u, const griddouble & bulK, const griddouble& G, 
const griddouble& gc, const int l, const double lc, const MPI_Setting & M);
//---------------- Compute Positive energy part at Gauss point
double Phi_positive(const double K_g, const double G_g, const grid2d& epsilon) ;
//---------------- Compute Negative energy part at Gauss point
double Phi_negative(const double K_g, const grid2d& epsilon);
//---------------- Compute H for each element 
void strain_history(Stack* U, griddouble & H, const gdvecdouble & u, const griddouble & bulK, 
const griddouble& G, const griddouble& gc, const MPI_Setting & M, const int l, 
const double lc, const double u_proportion) ;


#endif // !strain_history_H