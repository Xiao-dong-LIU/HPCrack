/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "PCG_routines.h"
#include "mpi.h"
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

/*
double vector_vector_product_multigrid3d(Stack *U, multigrid3d& r, grid3d& z_matrix, int l, MPI_Setting & M)
{
    Level* L;
    L=U->Ll+l;
    double x,y,z;
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
    #pragma omp parallel for reduction(+:s_per_MPI)
	for (int i = sx; i <= ex; i++)
		for (int j = sy; j <= ey; j++)
			for (int k = sz; k <= ez; k++)
            {
                s_per_MPI+=r(l-1,i,j,k)*z_matrix(i,j,k); // sum of each MPI
            }
    // do a sum over the whole domain
    MPI_Reduce(&s_per_MPI, &ss, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // send this sum to each mpi.
	MPI_Bcast(&ss, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return (ss);
}

double vector_vector_product_grid3d(Stack *U, grid3d& p, grid3d& w, int l, MPI_Setting & M)
{
    Level* L;
    L=U->Ll+l;
    double x,y,z;
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
    #pragma omp parallel for reduction(+:s_per_MPI)
	for (int i = sx; i <= ex; i++)
		for (int j = sy; j <= ey; j++)
			for (int k = sz; k <= ez; k++)
            {
                s_per_MPI+=(p(i,j,k)*w(i,j,k)); // sum of each MPI
				// if(w(i,j,k)!=0)
				// std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",w="<<w(i,j,k)<<",p="<<p(i,j,k)<<std::endl;
            }
    // do a sum over the whole domain
    MPI_Reduce(&s_per_MPI, &ss, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // send this sum to each mpi.
	MPI_Bcast(&ss, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return (ss);
}

void solution_update(Stack *U, multigrid3d& x, grid3d& p, double alpha, int l)
{
    Level *L;
    L=U->Ll+l;
    #pragma omp parallel for 
    for (int i=0;i<=L->nbcolox-1;i++)
        for (int j=0;j<=L->nbcoloy-1;j++)
            for (int k=0;k<=L->nbcoloz-1;k++)
			{
				x(l-1,i,j,k)=x(l-1,i,j,k)+alpha*p(i,j,k);
				// std::cout<<"x="<<x(l-1,i,j,k)<<std::endl;
			}
}

void p_update(Stack *U, grid3d& p, grid3d& z, double beta, int l)
{
    Level *L;
    L=U->Ll+l;
    #pragma omp parallel for 
    for (int i=0;i<=L->nbcolox-1;i++)
        for (int j=0;j<=L->nbcoloy-1;j++)
            for (int k=0;k<=L->nbcoloz-1;k++)
			{
				// if(i==17&&j==16&&k==24)
				// std::cout<<"before i="<<i<<",j="<<j<<",k="<<k<<",p="<<p(i,j,k)<<",z="<<z(i,j,k)<<std::endl;
				p(i,j,k)=z(i,j,k)+beta*p(i,j,k);
				// if(i==17&&j==16&&k==24)
				// std::cout<<"after i="<<i<<",j="<<j<<",k="<<k<<",p="<<p(i,j,k)<<",z="<<z(i,j,k)<<std::endl;
			}
}

void R_residual(Stack *U, grid3d &Fi, multigrid3d &f, multigrid3d & r, int l)
{
    Level* L;
	L = U->Ll + l;
    #pragma omp parallel for
    for (int i=0; i<=L->nbcolox-1;i++)
        for (int j=0; j<=L->nbcoloy-1;j++)
            for(int k=0;k<=L->nbcoloz-1;k++)
			{
				r(l-1,i,j,k)=f(l-1,i,j,k)-Fi(i,j,k);
				// std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",r="<<r(i,j,k)<<std::endl;
			}
                
}

void R_update(Stack *U, multigrid3d& r, grid3d& w, double alpha, int l)
{
    Level *L;
    L=U->Ll+l;
    #pragma omp parallel for 
    for (int i=0;i<=L->nbcolox-1;i++)
        for (int j=0;j<=L->nbcoloy-1;j++)
            for (int k=0;k<=L->nbcoloz-1;k++)
			{
				//std::cout<<"r_old="<<r(i,j,k)<<",w="<<w(i,j,k);
				r(l-1,i,j,k)=r(l-1,i,j,k)-alpha*w(i,j,k);
				//std::cout<<" r_new="<<r(i,j,k)<<std::endl;
			}
}

void z_solver(Stack *U, grid3d& M, multigrid3d& r, grid3d& z, int l)
{
    Level *L;
    L=U->Ll+l;
    #pragma omp parallel for
    for (int i=0; i<=L->nbcolox-1;i++)
        for (int j=0; j<=L->nbcoloy-1;j++)
            for(int k=0;k<=L->nbcoloz-1;k++)
			{
				z(i,j,k)=r(l-1,i,j,k)/M(i,j,k);
				// std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",M="<<M(i,j,k)<<",z="<<z(i,j,k)<<std::endl;
			}
}
*/