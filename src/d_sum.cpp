/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#include "d_sum.h"
#include "mpi.h"
#include "omp.h"
// find the maximum variation of d over the entier domain.
double delta_d_max(const grid <double> & d, const grid <double> & d_old)
{
    // initialisation of max d in each processor and the entire domain
    double d_max_proc=0.0, d_max=0.0;
    // find the max of d in each processor
    #pragma omp parallel for reduction(max:d_max_proc)
    for (int k = 0; k < d.axis()[2]; k++)
		 for (int j = 0; j <d.axis()[1]; j++)
			for (int i = 0; i <d.axis()[0]; i++)
            {
                if (d_max_proc<(d(i,j,k)-d_old(i,j,k)))
                    d_max_proc=d(i,j,k)-d_old(i,j,k);
	        // if (d_max_proc>0) std::cout<<"d_max_proc";
            }
//	std::cout<<d_max_proc;
  	MPI_Reduce(&d_max_proc, &d_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&d_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	if (d_max>0) std::cout<<d_max;
	return d_max;
}

/*
double d_sum(Stack *U, multigrid3d& d, int l, const MPI_Setting & M)
{
    Level* L;
    L=U->Ll+l;
	double ss = 0.0, s_per_MPI = 0.0; // sum of the entier domain and of each MPI
	// The start point in each mpi processer
	int sx = 1, sy = 1, sz = 1;
	// The end point in each mpi processer
	int ex = L->nbcolox - 2, ey = L->nbcoloy - 2, ez = L->nbcoloz - 2;
	if (M.dims[0] == 1) {    //one processor case
		sx = 0;
		ex = L->nbcolox-1;
	}
	else {
		if (M.coordinates[0] == 0) sx = 0;
		if (M.coordinates[0] == M.dims[0] - 1) ex = L->nbcolox-1;
	}
	if (M.dims[1] == 1) {
		sy = 0;
		ey = L->nbcoloy-1;
	}
	else {
		if (M.coordinates[1] == 0) sy = 0;
		if (M.coordinates[1] == M.dims[1] - 1) ey = L->nbcoloy-1;
	}
	if (M.dims[2] == 1) {
		sz = 0;
		ez = L->nbcoloz-1;
	}
	else {
		if (M.coordinates[2] == 0) sz = 0;
		if (M.coordinates[2] == M.dims[2] - 1) ez = L->nbcoloz-1;
	}

    // compute sum in each mpi task. 
    #pragma omp parallel for reduction(+:s_per_MPI)
	for (int i = sx; i <= ex; i++)
		for (int j = sy; j <= ey; j++)
			for (int k = sz; k <= ez; k++)
            {
                s_per_MPI+=d(l-1,i,j,k); // sum of each MPI
            }
    // do a sum over the whole domain
    MPI_Reduce(&s_per_MPI, &ss, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // send this sum to each mpi.
	MPI_Bcast(&ss, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return (ss);
}
*/



