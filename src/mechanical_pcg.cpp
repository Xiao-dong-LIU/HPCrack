/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "material_gauss_points.h"
#include "strain_history.h"
#include "shape_function.h"
#include "pcg_utilities.h"
#include "mechanical_stress.h"
#include "material_selection.h"
#include "material_coarse_element.h"
#include "mechanical_multigrid.h"
#include "mechanical_pcg.h"
#include "mechanical_stress.h"
#include "mechanical_strain.h"
#include "mpi_resource_cleanup.h"
#include "mpi_communications.h"
#include "boundary_conditions.h"

//------------ get material property at gauss point from the diplacement of previous time step
void MP_all_GAUSS_Point(Stack *U, vector <double> & K_d_g, vector <double> & G_d_g, 
const gridshort & MP_choice_d, const griddouble & bulK, const griddouble & G, 
const griddouble & d, const griddouble & K_element, const griddouble & G_element, 
const grid2d & N, const int l, const int i, const int j, const int k, const double small_k)
{
	double K_g,G_g;   // bulK and G at gauss point
	for (int g=0;g<8;g++)    /// Gauss Points
		{
			//--- get material property at gausse point 
			K_g = gauss_point_MP(bulK, N, i, j, k, g);
			G_g = gauss_point_MP(G, N, i, j, k, g);
			//--- get K_d and G_d for different cases
			// for K
			if (l==U->maxlevel)
			{
				if (MP_choice_d(i,j,k)==0)
				K_d_g[g]=K_G_d_g(d,bulK,N,i,j,k,g,small_k);
				else
                       		K_d_g[g]=K_g;
			}
			else 
			{
				 K_d_g[g]=K_element(i,j,k);	
			}
			// for G
			if (l==U->maxlevel)
			{
				if (MP_choice_d(i,j,k)==2)
                       		 G_d_g[g]=G_g;	
				else 
				 G_d_g[g]=K_G_d_g(d,G,N,i,j,k,g,small_k);	
			}
			else 
				G_d_g[g]=G_element(i,j,k);	
		}
}

void PreM_u_element(Stack *U, gdvecdouble & preM, const gridshort & MP_choice_d, 
const griddouble & bulK, const griddouble & G, const griddouble & d, const griddouble & K_element, 
const griddouble & G_element, const griddouble & DN, const grid2d & N, const double wg, const int l, 
const int i,const int j,const int k, const double small_k)
{
    ///l is the level
	int Cx[8]={0,1,1,0,0,1,1,0};
	int Cy[8]={0,0,1,1,0,0,1,1};
	int Cz[8]={0,0,0,0,1,1,1,1};
	vector <double> K_d_g(8,0.);
	vector <double> G_d_g(8,0.);
	MP_all_GAUSS_Point(U,K_d_g,G_d_g,MP_choice_d,bulK,G,d,K_element,G_element,N,l,i,j,k,small_k);
	for (int sf=0;sf<8;sf++)  /// shape Function
	{
		double stifu_ite=0.0,stifv_ite=0.0,stifw_ite=0.0;
		for (int g=0;g<8;g++)    /// Gauss Points
		{
			// compute the diagnol value of stiffness matrix
			stifu_ite+=wg*((K_d_g[g]-2./3.*G_d_g[g])*DN(0,g,sf)*DN(0,g,sf)+
							2.*G_d_g[g]*DN(0,g,sf)*DN(0,g,sf)+
							G_d_g[g]*(DN(1,g,sf)*DN(1,g,sf)+DN(2,g,sf)*DN(2,g,sf)));
			stifv_ite+=wg*((K_d_g[g]-2./3.*G_d_g[g])*DN(1,g,sf)*DN(1,g,sf)+
							2.*G_d_g[g]*DN(1,g,sf)*DN(1,g,sf)+
							G_d_g[g]*(DN(0,g,sf)*DN(0,g,sf)+DN(2,g,sf)*DN(2,g,sf)));
			stifw_ite+=wg*((K_d_g[g]-2./3.*G_d_g[g])*DN(2,g,sf)*DN(2,g,sf)+
							2.*G_d_g[g]*DN(2,g,sf)*DN(2,g,sf)+
							G_d_g[g]*(DN(0,g,sf)*DN(0,g,sf)+DN(1,g,sf)*DN(1,g,sf)));
		}
		//if(i==3&&j==3&&k==3)
		//std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",Fi_w_ite="<<Fi_w_ite<<std::endl;
		#pragma omp critical
		{
		preM.getgrid(0)(Cx[sf]+i,Cy[sf]+j,Cz[sf]+k)+=stifu_ite;
		preM.getgrid(1)(Cx[sf]+i,Cy[sf]+j,Cz[sf]+k)+=stifv_ite;
		preM.getgrid(2)(Cx[sf]+i,Cy[sf]+j,Cz[sf]+k)+=stifw_ite;
		}
	}
}

void PreM_u(Stack*U, gdvecdouble & preM, const gridshort & MP_choice_d, 
const griddouble & bulK, const griddouble & G, const griddouble & d, const griddouble & K_element, 
const griddouble & G_element, const griddouble & DN, const grid2d & N, const double wg, const int l, 
const double small_k)
{
    #pragma omp parallel for
	for (int k = 0; k <= d.axis()[2] - 2; k++)   /// Nb. of Elements
		for (int j = 0; j <= d.axis()[1] - 2; j++)
			for (int i = 0; i <= d.axis()[0] - 2; i++)
                PreM_u_element(U,preM,MP_choice_d,bulK,G,d,K_element,G_element,DN,N,wg,l,i,j,k,small_k);
}

void Internal_Force_element(Stack *U, gdvecdouble & Fi, const gdvecdouble & u,  
 const gridshort & MP_choice_d, const griddouble & bulK, const griddouble & G, 
 const griddouble & d,  const griddouble & K_element, const griddouble & G_element, 
 const griddouble & DN, const grid2d & N, const double wg, const int l, 
const int i,const int j,const int k, const double small_k)
{
	///l is the level
	int Cx[8]={0,1,1,0,0,1,1,0};
	int Cy[8]={0,0,1,1,0,0,1,1};
	int Cz[8]={0,0,0,0,1,1,1,1};
	
	grid2d epsilon(3, 3);
	grid2d sigma(3, 3);
	vector<size_t> sigmasize = {8,3,3};
	griddouble sigma_all_gauss(sigmasize);	   // stock sigma on 8 Gauss points to compute left hand
	vector <double> K_d_g(8,0.);
	vector <double> G_d_g(8,0.);
	double tr_epsilon,K_d_g_s,G_d_g_s;
	for (int g=0;g<8;g++)    /// Gauss Points
		{
			MP_all_GAUSS_Point(U,K_d_g,G_d_g,MP_choice_d,bulK,G,d,K_element,G_element,N,l,i,j,k,small_k);
			epsilon_routine(epsilon, u, DN, i, j, k, g);
			tr_epsilon = trace_epsilon(epsilon);
			K_d_g_s=K_d_g[g];
			G_d_g_s=G_d_g[g];
		
			Sigma_g(sigma,epsilon,K_d_g_s,G_d_g_s,tr_epsilon);
			for (int c=0; c<3; c++)
				for (int m=0; m<3; m++)
				{
					sigma_all_gauss(g,c,m)=sigma(c,m);
			
				}	
		}
	for (int sf=0;sf<8;sf++)  /// test Function
	{
		double Fi_u_ite=0.0,Fi_v_ite=0.0,Fi_w_ite=0.0;
		for (int g=0;g<8;g++)    /// Gauss Points
		{           
		   // compute the left hand of equation
			Fi_u_ite += wg*(sigma_all_gauss(g,0,0)*DN(0,g,sf)+
				  			sigma_all_gauss(g,1,0)*DN(1,g,sf)+
							sigma_all_gauss(g,2,0)*DN(2,g,sf));
			Fi_v_ite += wg*(sigma_all_gauss(g,0,1)*DN(0,g,sf)+
				  			sigma_all_gauss(g,1,1)*DN(1,g,sf)+
							sigma_all_gauss(g,2,1)*DN(2,g,sf));
			Fi_w_ite += wg*(sigma_all_gauss(g,0,2)*DN(0,g,sf)+
				  			sigma_all_gauss(g,1,2)*DN(1,g,sf)+
							sigma_all_gauss(g,2,2)*DN(2,g,sf));
		}
		
		#pragma omp critical
		{
		Fi.getgrid(0)(Cx[sf]+i,Cy[sf]+j,Cz[sf]+k)+=Fi_u_ite;
		Fi.getgrid(1)(Cx[sf]+i,Cy[sf]+j,Cz[sf]+k)+=Fi_v_ite;
		Fi.getgrid(2)(Cx[sf]+i,Cy[sf]+j,Cz[sf]+k)+=Fi_w_ite;
		}
	}
}

void Internal_Force(Stack *U, gdvecdouble & Fi, const gdvecdouble & u,  
const gridshort & MP_choice_d, const griddouble & bulK, const griddouble & G, 
const griddouble & d,  const griddouble & K_element, const griddouble & G_element, 
const griddouble & DN, const grid2d & N, const double wg, const int l, 
const double small_k)
{
	for (int di=0; di<3; di++)
    	Fi.getgrid(di).Assign_zero();
    #pragma omp parallel for
	for (int k = 0; k <= d.axis()[2] - 2; k++)   /// Nb. of Elements
		for (int j = 0; j <= d.axis()[1] - 2; j++)
			for (int i = 0; i <= d.axis()[0] - 2; i++)
			{
				Internal_Force_element(U,Fi,u,MP_choice_d,bulK,G,d,K_element,G_element,DN,N,wg,
				l,i,j,k,small_k);
		
			}
}

double Initial_residual_u(Stack *U, gdvecdouble & ru, const gdvecdouble & u,  
const gdvecdouble & fu, const gridshort & MP_choice_d, const griddouble & bulK, 
const griddouble & G, const griddouble & d,  const griddouble & K_element, 
const griddouble & G_element, MPI_Setting & M, const int l,  const double small_k)
{
	Level* L;
	L = U->Ll + l;
	gdvecdouble Fi(3,d.axis());
	vector<size_t> DNsize = {3,8,8};
    grid <double> DN(DNsize);
	grid2d N(8, 8);
	double r_norm=0.0;
	
	//   get shape function and its gradient
	Dshape_function(N, DN, L->hx, L->hy, L->hz);
	double wg = 0.125*L->hx*L->hy*L->hz;

	Internal_Force(U,Fi,u,MP_choice_d,bulK,G,d,K_element,G_element,DN,N,wg,l,small_k);

	for (int di=0; di<3; di++)
		ru.getgrid(di) = fu.getgrid(di) - Fi.getgrid(di);
	Dirichlet_residual_u(U,ru,M,l);
	// compute residual norm
	for (int di=0; di<3; di++)
		r_norm += grid_norm(ru.getgrid(di),M);
	
	return r_norm;
}


void pcg_u(Stack *U, gdvecdouble & u, gdvecdouble & ru, gdvecdouble & fu,
const gridshort & MP_choice_d, const griddouble & bulK, const griddouble & G, 
const griddouble & d,  const griddouble & K_element, const griddouble & G_element,
MPI_Setting & M, const int l,  const double small_k, const double Ut, const int NB_relax, 
const int t, double & f_norm)
{
    
	Level* L;
	L = U->Ll + l;
	// define anray shape function and its gradient
	vector<size_t> DNsize = {3,8,8};
    grid <double> DN(DNsize);
	grid2d N(8, 8);
	// preconditionned matrix 
	gdvecdouble preM(3,d.axis());
	gdvecdouble Fi(3,d.axis());
	gdvecdouble z(3,d.axis());
	gdvecdouble p(3,d.axis());

	double alpha = 0.,beta = 0.;
    double rz_old = 0.;
	double r_norm = 0.;
	double wg = 0.125*L->hx*L->hy*L->hz;
	//   get shape function and its gradient
	Dshape_function(N, DN, L->hx, L->hy, L->hz);
	// boundary condition of u
    Dirichlet_BD(U,u,M,l,Ut);   
    // get the initial residual
    Internal_Force(U,Fi,u,MP_choice_d,bulK,G,d,K_element,G_element,DN,N,wg,l,small_k);
    if (t==0&&l==U->maxlevel)
	{
		Dirichlet_force_u(U,fu,Fi,M,l);
		f_norm = 0.;
		for (int di=0; di<3; di++)
			f_norm += grid_norm(fu.getgrid(di),M);
	}
	
//	std::cout<<"fnorm="<<f_norm<<std::endl;
	for (int di=0; di<3; di++)
		ru.getgrid(di) = fu.getgrid(di) - Fi.getgrid(di);

	Dirichlet_residual_u(U,ru,M,l);


	// compute stiff
    PreM_u(U,preM,MP_choice_d,bulK,G,d,K_element,G_element,DN,N,wg,l,small_k);
	derived_datatypes(d.axis(),&M);
    for (int i=0;i<NB_relax;i++)
	{
        
		// get the solution by preconditionnor
		for (int di=0; di<3; di++)
		{
			z.getgrid(di) = ru.getgrid(di) / preM.getgrid(di);
			// update this solution
			SendRecv(z.getgrid(di),M);
		}
		
		// compute rz_new
		double rz_new = 0.;
		for (int di=0; di<3; di++)
			rz_new += vv_product(ru.getgrid(di),z.getgrid(di),M);

		//std::cout<<"rzu="<<rz_new_u<<",rzv="<<rz_new_v<<",rzw="<<rz_new_w<<std::endl;
		if (i==0) beta=0.0;
		else beta=rz_new/rz_old;
	//	cout<<"beta="<<beta;
		// update p
		for (int di=0; di<3; di++)
			p.getgrid(di) = z.getgrid(di) + beta*p.getgrid(di);
		// compute w
		Internal_Force(U,Fi,p,MP_choice_d,bulK,G,d,K_element,G_element,DN,N,wg,l,small_k);
		// compute pw
		double pw = 0.;
        for (int di=0; di<3; di++)
			pw += vv_product(p.getgrid(di),Fi.getgrid(di),M);
	//	cout<<" pw="<<pw;
		// compute alpha
		alpha=rz_new/pw;
	//	cout<<" alpha="<<alpha<<endl;
		for (int di=0; di<3; di++)
		{
			u.getgrid(di) = u.getgrid(di) + alpha*p.getgrid(di);
			ru.getgrid(di) = ru.getgrid(di) - alpha*Fi.getgrid(di);
			
		}
		Dirichlet_residual_u(U,ru,M,l);
		rz_old=rz_new;
		// compute residual norm
//		if (l==U->maxlevel)
		{
			r_norm = 0.;
			for (int di=0; di<3; di++)
				r_norm += grid_norm(ru.getgrid(di),M);
			if (M.coordinates[0] == 0 && M.coordinates[1] == 0 && M.coordinates[2] == 0) 
			{
				std::cout<<"Level "<<l<<" Ru ="<<r_norm<<",f="<<f_norm<<std::endl;
			}
		}
		/*
		if (r_norm_s/f_norm_s<=9.e-12) 
        {
          break;
        }
		*/
	}
	mpi_data_free(M);
}
