/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "PCG_d.h"
#include "shape_function.h"
#include "BD.h"
#include <cmath>
#include "mpi_communications.h"
#include "gauss_point_material_property.h"
#include "omp.h"
#include "mpi_free.h"
#include "PCG_routines.h"

void Md_element(grid <double> & Md, const grid <double> & H, const grid <double> & gc, 
const grid <double> & DN, grid2d & N, const double lc, const double wg, 
const int i, const int j, const int k)
{
    ///l is the level
	int Cx[8] = { 0,1,1,0,0,1,1,0 };
	int Cy[8] = { 0,0,1,1,0,0,1,1 };
	int Cz[8] = { 0,0,0,0,1,1,1,1 };
	double gc_g;
    for (int sf = 0; sf < 8; sf++)  /// Shape Function
		{
			double Md_ite = 0.0;
			for (int g = 0; g < 8; g++)    /// Gauss Points
			{
				// get material property at gauss point
				gc_g = gauss_point_MP(gc, N, i, j, k, g);
				// compute diagonal of stiffness matrix
				Md_ite += wg * ((2. * H(i, j, k) + gc_g / lc) * N(g, sf) * N(g, sf) +
						 		 gc_g * lc * (DN(0, g, sf) * DN(0, g, sf) +
									  		 DN(1, g, sf) * DN(1, g, sf) +
									  		 DN(2, g, sf) * DN(2, g, sf)));
            }
            #pragma omp critical  // avoid addition problem of openmp
			{
				Md(Cx[sf] + i, Cy[sf] + j, Cz[sf] + k) += Md_ite;
			}
        }
}

void Precondition_M(grid <double> & Md, const grid <double> & H, const grid <double> & gc, 
const grid <double> & DN, grid2d& N, const double lc, const double wg)
{
    
    Md.Assign_zero();
    #pragma omp parallel for
	for (int k = 0; k <= Md.axis()[2] - 2; k++) 
		for (int j = 0; j <= Md.axis()[1] - 2; j++)
			for (int i = 0; i <= Md.axis()[0] - 2; i++)	  /// Nb. of Elements
                Md_element(Md,H,gc,DN,N,lc,wg,i,j,k);
}

void Fld_element(grid <double> & Fi, const grid <double> & d, const grid <double> & H, 
const grid <double> & gc, const grid <double> & DN, const grid2d & N, const double lc, 
const double wg, const int i, const int j, const int k)
{
    ///l is the level
	int Cx[8] = { 0,1,1,0,0,1,1,0 };
	int Cy[8] = { 0,0,1,1,0,0,1,1 };
	int Cz[8] = { 0,0,0,0,1,1,1,1 };
	double gc_g;
	for (int sf = 0; sf < 8; sf++)  /// Shape Function
		{
			double Fi_ite = 0.0;
			for (int g = 0; g < 8; g++)    /// Gauss Points
			{
				// get material property at gauss point
				gc_g = gauss_point_MP(gc, N, i, j, k, g);		
				// compute the left hand of equation 
				for (int tf = 0; tf < 8; tf++)   ///  Test function
				{
					Fi_ite += wg * ((2.0 * H( i, j, k) + gc_g / lc) * N(g, sf) * N(g, tf) +
						  gc_g * lc * (DN(0, g, sf) * DN(0, g, tf) +
									   DN(1, g, sf) * DN(1, g, tf) +
									   DN(2, g, sf) * DN(2, g, tf))) * d(Cx[tf] + i, Cy[tf] + j, Cz[tf] + k);
				}
			}
            #pragma omp critical  // avoid addition problem of openmp
			{
				Fi(Cx[sf] + i, Cy[sf] + j, Cz[sf] + k) += Fi_ite;
			}
		}
}

void Fld(grid <double> & Fi, const grid <double> & d, const grid <double> & H, 
const grid <double> & gc, const grid <double> & DN, const grid2d & N, const double lc, 
const double wg)
{
    Fi.Assign_zero();
    #pragma omp parallel for
	for (int k = 0; k <= Fi.axis()[2] - 2; k++)
		for (int j = 0; j <= Fi.axis()[1] - 2; j++)
			for (int i = 0; i <= Fi.axis()[0] - 2; i++)  /// Nb. of Elements
                Fld_element(Fi,d,H,gc,DN,N,lc,wg,i,j,k);
}


double initial_residual(Stack *U, grid <double> & R, grid <double> & d, 
const grid <double> & fd,  const grid <double> & H, const grid <double> & gc, 
const int l, const double lc, const MPI_Setting & M)
{
         Level* L;
        L = U->Ll + l;
    // boundary condition of d
//	Dirichlet_BD_d(U,d,M,l);
   //     Dirichlet_BD_d(U, d, l, M);
        // define anray shape function and its gradient
		vector<size_t> DNsize = {3,8,8};
        grid <double> DN(DNsize);
        grid2d N(8, 8);
        // residual matrix
        grid<double> Fi(R); // interne force
        double f_norm=0.0,r_norm=0.0;
        //   get shape function and its gradient
        Dshape_function(N, DN, L->hx, L->hy, L->hz);
		double wg = 0.125 * L->hx * L->hy * L->hz;
    	// get the initial residual
    	Fld(Fi,d,H,gc,DN,N,lc,wg);
        R = fd - Fi;
   //     Dirichlet_residual_d(U,R,l,M);
        r_norm=grid_norm(R,M);
        return r_norm;
}

void pcg_d(Stack *U, grid<double> & d, grid<double> & fd, grid<double> & R, 
const grid<double> & H, const grid<double> & gc, MPI_Setting & M, const int l, 
const double lc,  const int NB_relax, const int t, double & fd_norm)
{
    Level* L;
	L = U->Ll + l;
    // boundary condition of d
//	Dirichlet_BD_d(U,d,M,l);
	// define anray shape function and its gradient
	vector<size_t> DNsize = {3,8,8};
    grid <double> DN(DNsize);
	grid2d N(8, 8);
	// preconditionned matrix 
	grid <double> preM(d.axis());
	// residual matrix
	grid <double> Fi(preM);
	grid <double> z(preM);
	grid <double> p(preM);

	double alpha=0.,beta=0.;
	double pw=0.0, rz_old=0.0, rz_new=0.0;
	double f_norm=0.0,r_norm=0.0;
	double wg = 0.125 * L->hx * L->hy * L->hz;
	//   get shape function and its gradient
	Dshape_function(N, DN, L->hx, L->hy, L->hz);
    // compute M
    Precondition_M(preM,H,gc,DN,N,lc,wg);
    // get the initial residual
    Fld(Fi,d,H,gc,DN,N,lc,wg);
	// Dirichlet_force_d(U,Fi,fd,l,M);
	fd_norm=grid_norm(fd,M);
	R = fd - Fi;
//	Dirichlet_residual_d(U,R,M,l);
	// compute residual norm
	/*
    r_norm=multigrid3d_norm(U,R,l,M);
    if (M.coordinates[0] == 0 && M.coordinates[1] == 0 && M.coordinates[2] == 0) 
        {
            std::cout<<"Level "<<l<<" Initial Relative residual="<<r_norm/f_norm<<std::endl;
        }
	*/
	derived_datatypes(z.axis(),&M);
    for (int i=0;i<NB_relax;i++)
	{
		// get the solution by preconditionnor
    	z = R/preM;
		// update this solution
		SendRecv(z,M);
		// compute rz_new
		rz_new=vv_product(R,z,M);
		if (i==0) beta=0.0;
		else beta=rz_new/rz_old;
		// update p
		p = z + beta*p;
		// compute w
		Fld(Fi,p,H,gc,DN,N,lc,wg);
		// compute pw
		pw=vv_product(p,Fi,M);
		// compute alpha
		alpha=rz_new/pw;
		d = d + alpha * p;
		R = R - alpha * Fi;
//		Dirichlet_residual_d(U,R,M,l);
		rz_old=rz_new;
		// compute residual norm
//		if (l==U->maxlevel)
		{
    	r_norm=grid_norm(R,M);
    	if (M.coordinates[0] == 0 && M.coordinates[1] == 0 && M.coordinates[2] == 0) 
        {
            std::cout<<"Level "<<l<<" Rd ="<<r_norm<<",fd="<<fd_norm<<std::endl;
        }
		/*
		if (r_norm/f_norm<=9.e-12) 
        {
          break;
        }
		*/
		}
	}
	mpi_data_free(M);
}
