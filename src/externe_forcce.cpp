/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "externe_forcce.h"
#include "omp.h"
#include "mpi.h"
#include "BD.h"
#include "PCG_u.h"
#include "shape_function.h"
#include "realdomain.h"


void externe_force_nodal(Stack *U, gdvecdouble & fu, gdvecdouble & u, const gridshort & MP_choice_d, 
const griddouble & bulK, const griddouble & G, const griddouble & d,  const griddouble & K_element, 
const griddouble & G_element, const MPI_Setting & M, const int l, const double Ut, const double small_k)
{
       Level* L;
       L = U->Ll + l;
        //  boundary condition of u
       Dirichlet_BD(U,u,M,l,Ut);
       // define anray shape function and its gradient
       vector<size_t> DNsize = {3,8,8};
       grid <double> DN(DNsize);
       grid2d N(8, 8);
       // left hand 
       gdvecdouble Fi(3,d.axis());
       //   get shape function and its gradient
       Dshape_function(N, DN, L->hx, L->hy, L->hz);
       double wg = 0.125*L->hx*L->hy*L->hz;
       // get the initial residual
       Internal_Force(U,Fi,u,MP_choice_d,bulK,G,d,K_element,G_element,DN,N,wg,l,small_k);
       Dirichlet_force_u(U,fu,Fi,M,l);
}
// ----- externe force on traction direction
double externe_force_sum(griddouble & f, const MPI_Setting & M)
{
    double x,y,z;
	double ss = 0.0, s_per_MPI = 0.0; // sum of the entier domain and of each MPI
	vector <size_t> realD;
	realD = realdomain(f.axis(),M);
    #pragma omp parallel for reduction(+:s_per_MPI) 
	for (int k = realD[4]; k < realD[5]; k++)
		for (int j = realD[2]; j < realD[3]; j++)
			for (int i = realD[0]; i < realD[1]; i++)
            {
               if (M.coordinates[2]==M.dims[2]-1&&k==f.axis()[2]-1)
		       s_per_MPI+=f(i,j,k); // sum of each MPI
            }
    // do a sum over the whole domain
    MPI_Reduce(&s_per_MPI, &ss, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // send this sum to each mpi.
	MPI_Bcast(&ss, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return (ss);
}
