/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "strain_history.h"
#include "epsilon.h"
#include "gauss_point_material_property.h"
#include "shape_function.h"
#include "Sigma.h"
#include "omp.h"
#include <cmath>
#include <algorithm>

using namespace std;

//---------------- Compute phi_c at Gauss point
double Phi_c(const double gc_g, const double lc) 
{
	double phi_c; 
	phi_c = 0.5*gc_g/lc;
	return phi_c;
}

//---------------- initial crack 
void initial_crack(Stack *U, griddouble & H, const MPI_Setting & M, const double H_crack, const int l)
{
Level *L,*Lmax;
L =U->Ll+l;
Lmax=U->Ll+U->maxlevel;
double x,y,z;
#pragma omp parallel for  private(x,y,z)
for (int k=0;k<=L->nbcoloz-1;k++)
{
  	/// z in each processor
  if (M.coordinates[2]==0) z=U->za+k*L->hz;
  else if (M.coordinates[2]==(M.dims[2]-1))
  z = U->za+(L->kk-L->nbcoloz+1+k)*L->hz;
  else
  z = U->za+(M.coordinates[2]*(L->nbcoloz-2) + k-1)*L->hz;
  for (int j=0;j<=L->nbcoloy-1;j++)
  {
  	  	/// y in each processor
	if (M.coordinates[1]==0) y=U->ya+j*L->hy;
	else if (M.coordinates[1]==(M.dims[1]-1))
	y = U->ya+(L->jj-L->nbcoloy+1+j)*L->hy;
	else
	y = U->ya+(M.coordinates[1]*(L->nbcoloy-2) + j-1)*L->hy;
    for (int i=0;i<=L->nbcolox-1;i++)
    {
		/// x in each processor
		if (M.coordinates[0]==0) x=U->xa+i*L->hx;
		else if (M.coordinates[0]==(M.dims[0]-1))
		x = U->xa+(L->ii-L->nbcolox+1+i)*L->hx;
		else
		x = U->xa+(M.coordinates[0]*(L->nbcolox-2) + i-1)*L->hx;
		// if (z==0.5*(U->za+U->zb))
		if (z==0.5*(U->za+U->zb)||z==0.5*(U->za+U->zb)-Lmax->hz)
		{
			if (y!=U->yb||x!=U->xb)
			H(i,j,k)=H_crack;
		}
	}
  }
}
}

//---------------- Compute elastic energy at Gauss point
double W(const grid2d& sigma, const grid2d& epsilon) 
{
	/*
	double tr_epsilon = 0.0;
	tr_epsilon = trace_epsilon(epsilon);
	double ddp_epsilon = 0.0; // Double dot product of the epsilon_dev
	ddp_epsilon = epsilon_dev(epsilon);
	double W_g; // Positive energy part
	W_g = 0.5 * K_g * tr_epsilon * tr_epsilon + G_g * ddp_epsilon;
	*/
	double W_g=0.;
	for (int c=0;c<3;c++)
		for (int m=0;m<3;m++)
		{
			W_g += 0.5*sigma(c,m)*epsilon(c,m);
		}
	return W_g;
}

//---------------- Compute the elastic displacement limit
double U_e(Stack* U, const gdvecdouble& u, const griddouble & bulK, const griddouble& G, 
const griddouble& gc, const int l, const double lc, const MPI_Setting & M)
{
	Level* L;
	L = U->Ll + l;
	double K_g, G_g, gc_g, phi_c, Ue_temp;
	double Ue=1e10,W0, Ue_min;
	double tr_epsilon;
	grid2d epsilon(3, 3);
	grid2d sigma(3, 3);
	vector<size_t> DNsize = {3,8,8};
    grid <double> DN(DNsize);
	grid2d N(8, 8);
	//   get shape function and its gradient
	Dshape_function(N, DN, L->hx, L->hy, L->hz);

	for (int i = 0; i <= u.getgrid(0).axis()[0] - 2; i++)
		for (int j = 0; j <= u.getgrid(0).axis()[1] - 2; j++)
			for (int k = 0; k <= u.getgrid(0).axis()[2] - 2; k++)   /// Nb. of Elements
			{  // element number loop, scince H is defined in element 
				for (int g = 0; g < 8; g++)
				{
					// compute material property at gauss points
					K_g = gauss_point_MP(bulK, N, i, j, k, g);
					G_g = gauss_point_MP(G, N, i, j, k, g);
					gc_g= gauss_point_MP(gc, N, i, j, k, g);
					// compute epsilon at gauss point 
					epsilon_routine(epsilon,u, DN, i, j, k, g);
					tr_epsilon=trace_epsilon(epsilon);
					Sigma_g(sigma,epsilon,K_g,G_g,tr_epsilon);
					// compute elastic energy at gauss point
					W0=W(sigma,epsilon);
					// Phi_c at gauss point
					phi_c=Phi_c(gc_g,lc);
					// std::cout<<"gc_g="<<gc_g<<",lc="<<lc<<std::endl;
					Ue_temp=sqrt(phi_c/W0);
					if (Ue_temp < Ue)
					{
						Ue=Ue_temp;
					}
				}
			}
	MPI_Reduce(&Ue, &Ue_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Ue_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return Ue_min;
}

//---------------- Compute Positive energy part at Gauss point
double Phi_positive(const double K_g, const double G_g, const grid2d& epsilon) 
{
	double tr_epsilon = 0.0;
	tr_epsilon = trace_epsilon(epsilon);
	double ddp_epsilon = 0.0; // Double dot product of the epsilon_dev
	ddp_epsilon = epsilon_dev(epsilon);
	double phi_positive; // Positive energy part
	if (tr_epsilon >= 0)
		phi_positive = 0.5 * K_g * tr_epsilon * tr_epsilon + G_g * ddp_epsilon;
	else
		phi_positive = G_g * ddp_epsilon;
	return phi_positive;
}

//---------------- Compute Negative energy part at Gauss point
double Phi_negative(const double K_g, const grid2d& epsilon) 
{
	double tr_epsilon = 0.0;
	tr_epsilon = trace_epsilon(epsilon);
	double phi_negative; // Positive energy part
	if (tr_epsilon >= 0)
		phi_negative = 0.0;
	else
		phi_negative = 0.5 * K_g * tr_epsilon * tr_epsilon;
	return phi_negative;
}

//---------------- Compute H for each element 
void strain_history(Stack* U, griddouble & H, const gdvecdouble & u, const griddouble & bulK, 
const griddouble& G, const griddouble& gc, const MPI_Setting & M, const int l, 
const double lc, const double u_proportion) 
{
	Level* L;
	L = U->Ll + l;
	double K_g, G_g, gc_g, phi_positive, phi_c, tmp, phi_c_tmp, H_tmp;
	grid2d epsilon(3, 3);
	vector<size_t> DNsize = {3,8,8};
    grid <double> DN(DNsize);
	grid2d N(8, 8);
	//   get shape function and its gradient
	Dshape_function(N, DN, L->hx, L->hy, L->hz);
	double tr_epsilon,W_g;
	grid2d sigma(3, 3);

	for (int i = 0; i <= L->nbcolox - 2; i++)
		for (int j = 0; j <= L->nbcoloy - 2; j++)
			for (int k = 0; k <= L->nbcoloz - 2; k++)   /// Nb. of Elements
			{  // element number loop, scince H is defined in element 
				tmp=0.;
				for (int g = 0; g < 8; g++)
				{
					// compute material property at gauss points
					K_g = gauss_point_MP(bulK, N, i, j, k, g);
					G_g = gauss_point_MP(G, N, i, j, k, g);
					gc_g= gauss_point_MP(gc, N, i, j, k, g);
					// compute epsilon at gauss point 
					epsilon_routine(epsilon, u, DN, i, j, k, g);
					
					// compute epsilon with the correct displqcement
					for (int c=0;c<3;c++)
						for (int m=0;m<3;m++)
							epsilon(c,m)=u_proportion*epsilon(c,m);
					// compute phi+ at gauss point
					 phi_positive = Phi_positive(K_g, G_g, epsilon)-0.5*gc_g/lc;
					 // phi_c=Phi_c(gc_g,lc);
					 tmp+=phi_positive;
				}
				H_tmp = 0.125*tmp;
//				H_tmp = max(H_tmp, 0.0);
				if (H( i, j, k) < H_tmp)
				{
					// find the max of H in each element for all the history 
					H( i, j, k) = H_tmp;
				}
			}
	// SendRecv(H,l-1,L->nbcolox,L->nbcoloy,L->nbcoloz,M);
}
