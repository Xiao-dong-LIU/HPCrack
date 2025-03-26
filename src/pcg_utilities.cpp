/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "pcg_utilities.h"
#include "omp.h"
#include "realdomain.h"



double vv_product(const grid<double>& r, const grid<double>& z, const MPI_Setting & M)
{
	double ss = 0.0, s_per_MPI = 0.0; // sum of the entier domain and of each MPI
	vector <size_t> realD;
	realD = realdomain(r.axis(),M);
    #pragma omp parallel for reduction(+:s_per_MPI)
	for (int k = realD[4]; k < realD[5]; k++)
		for (int j = realD[2]; j < realD[3]; j++)
			for (int i = realD[0]; i < realD[1]; i++)
            {
                s_per_MPI+=r(i,j,k)*z(i,j,k); // sum of each MPI
            }
    // do a sum over the whole domain
    MPI_Reduce(&s_per_MPI, &ss, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // send this sum to each mpi.
	MPI_Bcast(&ss, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return (ss);
}

double grid_norm(const grid<double>& r, const MPI_Setting & M)
{
	double ss = 0.0, s_per_MPI = 0.0; // sum of the entier domain and of each MPI
	vector <size_t> realD;
	realD = realdomain(r.axis(),M);
	/*
	for (int x=0; x<realD.size(); x++)
		cout<<realD[x]<< " ";
	cout<<endl;
	*/
    #pragma omp parallel for reduction(+:s_per_MPI)
	for (int k = realD[4]; k < realD[5]; k++)
		for (int j = realD[2]; j < realD[3]; j++)
			for (int i = realD[0]; i < realD[1]; i++)
            {
                s_per_MPI+=r(i,j,k)*r(i,j,k); // sum of each MPI
		//		cout<<"i="<<i<<",j="<<j<<",k="<<k<<",r="<<r(i,j,k)<<endl;
            }
    // do a sum over the whole domain
    MPI_Reduce(&s_per_MPI, &ss, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // send this sum to each mpi.
	MPI_Bcast(&ss, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return (ss);
}
