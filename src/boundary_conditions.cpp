/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "boundary_conditions.h"
#include <cmath>

double get_coordinate(double origin, double h, int index, int procIndex, int dim, int total, int nbcolo) {
    if (procIndex == 0) return origin + index * h;
    else if (procIndex == dim - 1) return origin + (total - nbcolo + 1 + index) * h;
    else return origin + (procIndex * (nbcolo - 2) + index - 1) * h;
}


void Dirichlet_BD(Stack *U, gdvecdouble & u, const MPI_Setting & M, const int l, const double Ut)
// delta_u is the displacement incremant 
// Ue is the elastic displacement 
// t is the time step
{
Level *L;
L =U->Ll+l;
double x,y,z;
#pragma omp parallel for  private(x,y,z)
for (int k=0;k<=L->nbcoloz-1;k++)
{
  	/// z in each processor
  z = get_coordinate(U->za, L->hz, k, M.coordinates[2], M.dims[2], L->kk, L->nbcoloz);
  for (int j=0;j<=L->nbcoloy-1;j++)
  {
  	  	/// y in each processor
	y = get_coordinate(U->ya, L->hy, j, M.coordinates[1], M.dims[1], L->jj, L->nbcoloy);
    for (int i=0;i<=L->nbcolox-1;i++)
    {
		/// x in each processor
		x = get_coordinate(U->xa, L->hx, i, M.coordinates[0], M.dims[0], L->ii, L->nbcolox);
		// block points
		
		if (z==U->za&&x==0.5*(U->xa+U->xb)&&y==U->ya)
		{
		//	if (direction==0)
			u.getgrid(0)(i,j,k)=0.0;
		}
		if (z==U->za&&x==U->xb&&y==0.5*(U->ya+U->yb))
		{
		//	if (direction==1)
			u.getgrid(1)(i,j,k)=0.0;
		}
		// prescribed displacement en Z
		if (z==U->za)
		{
		//	if (direction==2)
			u.getgrid(2)(i,j,k)=0.0;			
		}
		if (z==U->zb)
		{
		//	if (direction==2)
			u.getgrid(2)(i,j,k)=Ut;
			//std::cout<<"ok"<<std::endl;			
		}
		}
	}
}
}
void Dirichlet_force_u(Stack* U, gdvecdouble& f, const gdvecdouble & Fi, const MPI_Setting & M, 
const int l)
{
	Level* L;
	L = U->Ll + l;
	double x,y,z;
	#pragma omp parallel for  private(x,y,z)
	for (int k=0;k<=L->nbcoloz-1;k++)
	{
		/// z in each processor
		z = get_coordinate(U->za, L->hz, k, M.coordinates[2], M.dims[2], L->kk, L->nbcoloz);
		for (int j=0;j<=L->nbcoloy-1;j++)
		{
				/// y in each processor
			y = get_coordinate(U->ya, L->hy, j, M.coordinates[1], M.dims[1], L->jj, L->nbcoloy);
			for (int i=0;i<=L->nbcolox-1;i++)
			{
				/// x in each processor
				x = get_coordinate(U->xa, L->hx, i, M.coordinates[0], M.dims[0], L->ii, L->nbcolox);
				if (z==U->za&&x==0.5*(U->xa+U->xb)&&y==U->ya)
				{
				//	if (direction==0)
					f.getgrid(0)(i,j,k)=Fi.getgrid(0)(i,j,k);
				}
				if (z==U->za&&x==U->xb&&y==0.5*(U->ya+U->yb))
				{
				//	if (direction==1)
					f.getgrid(1)(i,j,k)=Fi.getgrid(1)(i,j,k);
				}
				if (z==U->za)
				{
				//	if (direction==2)
					f.getgrid(2)(i,j,k)=Fi.getgrid(2)(i,j,k);		
				}
				if (z==U->zb)
				{
				//	if (direction==2)
					f.getgrid(2)(i,j,k)=Fi.getgrid(2)(i,j,k);
				}	
			}
		}
	}
}
void Dirichlet_residual_u(Stack* U, gdvecdouble & r, const MPI_Setting & M, const int l)
{
	Level* L;
	L = U->Ll + l;
	double x,y,z;
	#pragma omp parallel for  private(x,y,z)
	for (int k=0;k<=L->nbcoloz-1;k++)
	{
		/// z in each processor
		z = get_coordinate(U->za, L->hz, k, M.coordinates[2], M.dims[2], L->kk, L->nbcoloz);
		for (int j=0;j<=L->nbcoloy-1;j++)
		{
				/// y in each processor
			y = get_coordinate(U->ya, L->hy, j, M.coordinates[1], M.dims[1], L->jj, L->nbcoloy);
			for (int i=0;i<=L->nbcolox-1;i++)
			{
				/// x in each processor
				x = get_coordinate(U->xa, L->hx, i, M.coordinates[0], M.dims[0], L->ii, L->nbcolox);

				if (z==U->za&&x==0.5*(U->xa+U->xb)&&y==U->ya)
				{
				//	if (direction==0)
					r.getgrid(0)(i,j,k)=0.0;
				}
				if (z==U->za&&x==U->xb&&y==0.5*(U->ya+U->yb))
				{
				//	if (direction==1)
					r.getgrid(1)(i,j,k)=0.0;
				}
				if (z==U->za)
				{
				//	if (direction==2)
					r.getgrid(2)(i,j,k)=0.0;		
				}
				if (z==U->zb)
				{
				//	if (direction==2)
					r.getgrid(2)(i,j,k)=0.0;
				}	
			}
		}
	}
}


void Dirichlet_BD_d(Stack* U, griddouble & d, const MPI_Setting & M, const int l)
{
	/*
	Level* L;
	L = U->Ll + l;
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
				
				if (x<=(U->xa+0.1*(U->xb-U->xa)))
				{
					if(z==0.5*(U->za+U->zb)||z==0.5*(U->za+U->zb)+L->hz||z==0.5*(U->za+U->zb)-L->hz)
					d(i,j,k)=1.0;
				}
				
			}
		}
	}
	*/
}

void Dirichlet_force_d(Stack* U, griddouble& fd, const griddouble & Fi,  const MPI_Setting & M, const int l)
{
	/*
	Level* L;
	L = U->Ll + l;
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
				
				if (x<=(U->xa+0.1*(U->xb-U->xa)))
				{
					if(z==0.5*(U->za+U->zb)||z==0.5*(U->za+U->zb)+L->hz||z==0.5*(U->za+U->zb)-L->hz)
					fd(i, j, k)=Fi( i, j, k);
				}	
				
			}
		}
	}
	*/
}
void Dirichlet_residual_d(Stack* U, griddouble & r, const MPI_Setting & M, const int l)
{
	/*
	Level* L;
	L = U->Ll + l;
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
				
				if (x<=(U->xa+0.1*(U->xb-U->xa)))
				{
					if(z==0.5*(U->za+U->zb)||z==0.5*(U->za+U->zb)+L->hz||z==0.5*(U->za+U->zb)-L->hz)
					r(i,j,k)=0.;
				}	
				
			}
		}
	}
	*/
}

