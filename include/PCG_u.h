/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _PCG_U_H
#define _PCG_U_H
#include "grid2d.h"
#include "mgvec.h"


void MP_all_GAUSS_Point(Stack *U, vector <double> & K_d_g, vector <double> & G_d_g, 
const gridshort & MP_choice_d, const griddouble & bulK, const griddouble & G, 
const griddouble & d, const griddouble & K_element, const griddouble & G_element, 
const grid2d & N, const int l, const int i, const int j, const int k, const double small_k);
void PreM_u_element(Stack *U, gdvecdouble & preM, const gridshort & MP_choice_d, 
const griddouble & bulK, const griddouble & G, const griddouble & d, const griddouble & K_element, 
const griddouble & G_element, const griddouble & DN, const grid2d & N, const double wg, const int l, 
const int i,const int j,const int k, const double small_k);

void PreM_u(Stack*U, gdvecdouble & preM, const gridshort & MP_choice_d, 
const griddouble & bulK, const griddouble & G, const griddouble & d, const griddouble & K_element, 
const griddouble & G_element, const griddouble & DN, const grid2d & N, const double wg, const int l, 
const double small_k);
void Internal_Force_element(Stack *U, gdvecdouble & Fi, const gdvecdouble & u,  
 const gridshort & MP_choice_d, const griddouble & bulK, const griddouble & G, 
 const griddouble & d,  const griddouble & K_element, const griddouble & G_element, 
 const griddouble & DN, const grid2d & N, const double wg, const int l, 
const int i,const int j,const int k, const double small_k);

// --------- The left hand 
void Internal_Force(Stack *U, gdvecdouble & Fi, const gdvecdouble & u,  
const gridshort & MP_choice_d, const griddouble & bulK, const griddouble & G, 
const griddouble & d,  const griddouble & K_element, const griddouble & G_element, 
const griddouble & DN, const grid2d & N, const double wg, const int l, 
const double small_k);
// --------- initial residual for displacement 
double Initial_residual_u(Stack *U, gdvecdouble & ru, const gdvecdouble & u,  
const gdvecdouble & fu, const gridshort & MP_choice_d, const griddouble & bulK, 
const griddouble & G, const griddouble & d,  const griddouble & K_element, 
const griddouble & G_element, MPI_Setting & M, const int l,  const double small_k);
// --------- pcg solveur for displacement 
void pcg_u(Stack *U, gdvecdouble & u, gdvecdouble & ru, gdvecdouble & fu,
const gridshort & MP_choice_d, const griddouble & bulK, const griddouble & G, 
const griddouble & d,  const griddouble & K_element, const griddouble & G_element,
MPI_Setting & M, const int l,  const double small_k, const double Ut, const int NB_relax, 
const int t, double & f_norm);

#endif