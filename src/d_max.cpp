/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#include "d_max.h"
#include "mpi.h"
#include "omp.h"
#include <algorithm>

// find the maximum variation of d over the entier domain.
double delta_d_max(const grid <double> & d, const grid <double> & d_old)
{
    double delatad_max_proc=0.0, delatad_max=0.0;
    // find the max of the variation of d in each processor
    #pragma omp parallel for reduction(max:delatad_max_proc)
    for (int k = 0; k < d.axis()[2]; k++)
		 for (int j = 0; j <d.axis()[1]; j++)
			for (int i = 0; i <d.axis()[0]; i++)
            {
                if (delatad_max_proc<(d(i,j,k)-d_old(i,j,k)))
                    delatad_max_proc=d(i,j,k)-d_old(i,j,k);
	        // if (d_max_proc>0) std::cout<<"d_max_proc";
            }
//	std::cout<<d_max_proc;
  	MPI_Reduce(&delatad_max_proc, &delatad_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&delatad_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	if (d_max>0) std::cout<<d_max;
	return delatad_max;
}

// find max of d over the whole domain
double d_max (const grid <double> & d)
{
	double d_max=0.,d_max_proc=0.;
	    // find the max of d in each processor
    #pragma omp parallel for reduction(max:d_max_proc)
    for (int k = 0; k < d.axis()[2]; k++)
		 for (int j = 0; j <d.axis()[1]; j++)
			for (int i = 0; i <d.axis()[0]; i++)
            {
                if (d_max_proc<d(i,j,k))
                    d_max_proc=d(i,j,k);
            }
  	MPI_Reduce(&d_max_proc, &d_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&d_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return d_max;
}



