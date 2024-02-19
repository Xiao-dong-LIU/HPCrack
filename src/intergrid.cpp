/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "BD.h"
#include "intergrid.h"
#include "mpi_communications.h"
#include "shape_function.h"
#include "PCG_u.h"
#include "PCG_d.h"
#include "omp.h"
#include <stdlib.h>
#include <cmath>
#include "mpi_free.h"


using namespace std;

//------------- Injection of solution on level l to level l-1
void solution_injection(grid<double> &uc, const grid<double> &uf, const MPI_Setting &M)
{
	/// u represent solutions
	/// compute initial approximation on level l-1
	/// in coarsening step from level l
	/// imput value direction (0,1,2) is for vector u. set to be 0 for d.
	/// MP is material property

	int x, y, z, a;
#pragma omp parallel for private(x, y, z)
	for (int kc = 0; kc < uc.axis()[2]; kc++)
		for (int jc = 0; jc < uc.axis()[1]; jc++)
			for (int ic = 0; ic < uc.axis()[0]; ic++)
			{
				if (M.dims[0] == 1)
					x = 2 * ic;
				else
				{
					if (M.coordinates[0] == 0)
						x = 2 * ic;
					else
						x = 2 * ic - 1;
				}
				if (M.dims[1] == 1)
					y = 2 * jc;
				else
				{
					if (M.coordinates[1] == 0)
						y = 2 * jc;
					else
						y = 2 * jc - 1;
				}
				if (M.dims[2] == 1)
					z = 2 * kc;
				else
				{
					if (M.coordinates[2] == 0)
						z = 2 * kc;
					else
						z = 2 * kc - 1;
				}
				if (x < 0)
					x = 0;
				if (y < 0)
					y = 0;
				if (z < 0)
					z = 0;
				// skip nodes without material property
				/*
				a=1;
				for(int q=max(0,z-1);q<=min(z+1,L->nbcoloz-1);q++)
				for(int o=max(0,x-1);o<=min(x+1,L->nbcolox-1);o++)
					for(int p=max(0,y-1);p<=min(y+1,L->nbcoloy-1);p++)
						{
							if(fabs(MP(l-1,o,p,q))<1.0e-6)
							{
								a=0;
								break;
							}
						}
				if (a==1)
				*/
				uc(ic, jc, kc) = uf(x, y, z);
		//		cout<<"i= "<<ic<<" j= "<<jc<<" k= "<<kc<<" u= "<<uc(ic, jc, kc)<<endl;
			}
}
//-------------- injection for d
void coarsen_d(grid<double> &dc, grid<double> &doldc, const grid<double> &df, MPI_Setting &M)
{
	solution_injection(dc, df, M);
	// reset dirichelet boundary condition
	// Dirichlet_BD_d(U,dc,M,l);
	// update gohst points
	derived_datatypes(dc.axis(),&M);
	SendRecv(dc, M);
	/// store coarse grid solution in uold array
	doldc = dc;
	mpi_data_free(M);
}

//------------- Injection for displacement
void coarsen_u(Stack *U, gdvecdouble &uc, gdvecdouble &ucold, const gdvecdouble &uf,
			   const griddouble &bulK, MPI_Setting &M, const int lc, double Ut)
{
	for (int di = 0; di < 3; di++)
		solution_injection(uc.getgrid(di), uf.getgrid(di), M);
	// reset dirichelet boundary condition
	Dirichlet_BD(U, uc, M, lc, Ut);
	// update gohst points
	derived_datatypes(uc.getgrid(0).axis(),&M);
	for (int di = 0; di < 3; di++)
		SendRecv(uc.getgrid(di), M);
	/// store coarse grid solution in uold array
	for (int di = 0; di < 3; di++)
		ucold.getgrid(di) = uc.getgrid(di);
	mpi_data_free(M);
}

// standard restriction operator
void coarsen_f_standard_routine(grid<double> &fc, const grid<double> &rf,
								const grid<double> &Fic, const MPI_Setting &M)
{
	double r11, r12, r13, r14, r15, r16, r17, r18;
	double r21, r22, r23, r24, r25, r26, r27, r28, r29, r210, r211, r212;
	double r41, r42, r43, r44, r45, r46;
	double r8;
	double juck = 1.0e-6;
	int x, y, z;
	int Xlimit = rf.axis()[0] - 2;
	int Ylimit = rf.axis()[1] - 2;
	int Zlimit = rf.axis()[2] - 2;

	int sx, sy, sz;
	if (M.dims[0] == 1)
		sx = 0;
	else
	{
		if (M.coordinates[0] == 0)
			sx = 0;
		else
			sx = 1;
	}
	if (M.dims[1] == 1)
		sy = 0;
	else
	{
		if (M.coordinates[1] == 0)
			sy = 0;
		else
			sy = 1;
	}
	if (M.dims[2] == 1)
		sz = 0;
	else
	{
		if (M.coordinates[2] == 0)
			sz = 0;
		else
			sz = 1;
	}
#pragma omp parallel for private(x, y, z, r11, r12, r13, r14, r15, r16, r17, r18,               \
								 r21, r22, r23, r24, r25, r26, r27, r28, r29, r210, r211, r212, \
								 r41, r42, r43, r44, r45, r46, r8)
	for (int ic = sx; ic < fc.axis()[0]; ic++)
		for (int jc = sy; jc < fc.axis()[1]; jc++)
			for (int kc = sz; kc < fc.axis()[2]; kc++)
			{
				if (M.dims[0] == 1)
					x = 2 * ic;
				else
				{
					if (M.coordinates[0] == 0)
						x = 2 * ic;
					else
						x = 2 * ic - 1;
				}
				if (M.dims[1] == 1)
					y = 2 * jc;
				else
				{
					if (M.coordinates[1] == 0)
						y = 2 * jc;
					else
						y = 2 * jc - 1;
				}
				if (M.dims[2] == 1)
					z = 2 * kc;
				else
				{
					if (M.coordinates[2] == 0)
						z = 2 * kc;
					else
						z = 2 * kc - 1;
				}
				if (x < 1 || y < 1 || z < 1)
					r11 = 0.0;
				else
					r11 = rf(x - 1, y - 1, z - 1);

				if (x < 1 || y < 1 || z > Zlimit)
					r12 = 0.0;
				else
					r12 = rf(x - 1, y - 1, z + 1);

				if (x < 1 || y > Ylimit || z < 1)
					r13 = 0.0;
				else
					r13 = rf(x - 1, y + 1, z - 1);

				if (x < 1 || y > Ylimit || z > Zlimit)
					r14 = 0.0;
				else
					r14 = rf(x - 1, y + 1, z + 1);

				if (x > Xlimit || y < 1 || z < 1)
					r15 = 0.0;
				else
					r15 = rf(x + 1, y - 1, z - 1);

				if (x > Xlimit || y < 1 || z > Zlimit)
					r16 = 0.0;
				else
					r16 = rf(x + 1, y - 1, z + 1);

				if (x > Xlimit || y > Ylimit || z < 1)
					r17 = 0.0;
				else
					r17 = rf(x + 1, y + 1, z - 1);

				if (x > Xlimit || y > Ylimit || z > Zlimit)
					r18 = 0.0;
				else
					r18 = rf(x + 1, y + 1, z + 1);

				if (x < 1 || y < 1)
					r21 = 0.0;
				else
					r21 = rf(x - 1, y - 1, z);

				if (x < 1 || y > Ylimit)
					r22 = 0.0;
				else
					r22 = rf(x - 1, y + 1, z);

				if (x > Xlimit || y < 1)
					r23 = 0.0;
				else
					r23 = rf(x + 1, y - 1, z);

				if (x > Xlimit || y > Ylimit)
					r24 = 0.0;
				else
					r24 = rf(x + 1, y + 1, z);

				if (x < 1 || z < 1)
					r25 = 0.0;
				else
					r25 = rf(x - 1, y, z - 1);

				if (x < 1 || z > Zlimit)
					r26 = 0.0;
				else
					r26 = rf(x - 1, y, z + 1);

				if (x > Xlimit || z < 1)
					r27 = 0.0;
				else
					r27 = rf(x + 1, y, z - 1);

				if (x > Xlimit || z > Zlimit)
					r28 = 0.0;
				else
					r28 = rf(x + 1, y, z + 1);

				if (y < 1 || z < 1)
					r29 = 0.0;
				else
					r29 = rf(x, y - 1, z - 1);

				if (y < 1 || z > Zlimit)
					r210 = 0.0;
				else
					r210 = rf(x, y - 1, z + 1);

				if (y > Ylimit || z < 1)
					r211 = 0.0;
				else
					r211 = rf(x, y + 1, z - 1);

				if (y > Ylimit || z > Zlimit)
					r212 = 0.0;
				else
					r212 = rf(x, y + 1, z + 1);

				if (z < 1)
					r41 = 0.0;
				else
					r41 = rf(x, y, z - 1);

				if (z > Zlimit)
					r42 = 0.0;
				else
					r42 = rf(x, y, z + 1);

				if (y < 1)
					r43 = 0.0;
				else
					r43 = rf(x, y - 1, z);

				if (y > Ylimit)
					r44 = 0.0;
				else
					r44 = rf(x, y + 1, z);

				if (x < 1)
					r45 = 0.0;
				else
					r45 = rf(x - 1, y, z);

				if (x > Xlimit)
					r46 = 0.0;
				else
					r46 = rf(x + 1, y, z);

				r8 = rf(x, y, z);

				/// FAS coarse grid right hand side
				/// with full weighting of residuals

				fc(ic, jc, kc) = Fic(ic, jc, kc) + 0.125 * 
				(r11 + r12 + r13 + r14 + r15 + r16 + r17 + r18 + 
				2.0 * (r21 + r22 + r23 + r24 + r25 + r26 + r27 + r28 + r29 + r210 + r211 + r212) + 
				4.0 * (r41 + r42 + r43 + r44 + r45 + r46) + 
				8.0 * r8);
			//	cout<<"i= "<<ic<<" j= "<<jc<<" k= "<<kc<<" f= "<<fc(ic, jc, kc)
			//	<<" Fi="<<Fic(ic, jc, kc)<<endl;
			}
}

// restriction for d
void coarsen_f_d(Stack *U, grid<double> &fc, const grid<double> &dc, grid<double> & rf,
const grid<double> &gcc, const grid<double> &Hc, MPI_Setting &M, const double lc, const int lf)
{
	/// compute coarse grid right hand side on level l-1
	/// in coarsening step from level l

	Level *Lc;

	Lc = U->Ll + lf - 1;
	/// update residual on fine grid
	derived_datatypes(rf.axis(),&M);
	SendRecv(rf, M);
	mpi_data_free(M);
	grid<double> Fic(dc.axis());
	vector<size_t> DNsize = {3, 8, 8};
	grid<double> DN(DNsize);
	grid2d N(8, 8);
	// coarse grid shape function
	Dshape_function(N, DN, Lc->hx, Lc->hy, Lc->hz);
	/// compute interne force on coarse grid
	double wg = 0.125*Lc->hx * Lc->hy * Lc->hz;
	Fld(Fic, dc, Hc, gcc, DN, N, lc, wg);
	coarsen_f_standard_routine(fc, rf, Fic, M);
	derived_datatypes(fc.axis(),&M);
	SendRecv(fc, M);
	mpi_data_free(M);
}

/// restriction for u
void coarsen_fu(Stack *U, gdvecdouble &fc, const gdvecdouble &uc, const gridshort &MP_choice_dc, 
gdvecdouble &rf, const griddouble &bulK, const griddouble &G, const griddouble &d, 
const griddouble &K_element, const griddouble &G_element, MPI_Setting &M, const int lf, 
double small_k)
{
	/// compute coarse grid right hand side on level l-1
	/// in coarsening step from level l

	Level *Lc;
	Lc = U->Ll + lf-1;
	double wg = 0.125*Lc->hx * Lc->hy * Lc->hz;
	vector<size_t> DNsize = {3, 8, 8};
	grid<double> DN(DNsize);
	grid2d N(8, 8);
	Dshape_function(N, DN, Lc->hx, Lc->hy, Lc->hz);

	/// update residual on fine grid
	derived_datatypes(rf.getgrid(0).axis(),&M);
	for (int di = 0; di < 3; di++)
		SendRecv(rf.getgrid(di), M);
	mpi_data_free(M);
	/// Internal force on coarse grid
	gdvecdouble Fic(3, d.axis());
	Internal_Force(U, Fic, uc, MP_choice_dc, bulK, G, d, K_element, G_element, DN, N, wg, lf-1, small_k);

	// Standard MG restriction operators
	for (int di = 0; di < 3; di++)
		coarsen_f_standard_routine(fc.getgrid(di), rf.getgrid(di), Fic.getgrid(di), M);
	//	coarsen_f_crack(U,ru,fu,d,Fi_u,l,M);
	/*
		// Alcoiff type MG restriction operators
		coarsen_f_alcouff_routine(U,ru,fu,G_element,Fi_u,l,M);
		coarsen_f_alcouff_routine(U,rv,fv,G_element,Fi_v,l,M);
		coarsen_f_alcouff_routine(U,rw,fw,G_element,Fi_w,l,M);
	*/
	derived_datatypes(fc.getgrid(0).axis(),&M);
	for (int di = 0; di < 3; di++)
		SendRecv(fc.getgrid(di), M);
	mpi_data_free(M);
}

//------------- standard prologation operator
void refine_standard_routine(grid<double> &uf, const grid<double> &uc,
							 const grid<double> &ucold, const MPI_Setting &M)
{
	/// Interpolation and addition of coarse grid correction from grid l-1
	/// to grid l

	int sx, sy, sz;
	if (M.dims[0] == 1)
		sx = 0;
	else
	{
		if (M.coordinates[0] == 0)
			sx = 0;
		else
			sx = 1;
	}
	if (M.dims[1] == 1)
		sy = 0;
	else
	{
		if (M.coordinates[1] == 0)
			sy = 0;
		else
			sy = 1;
	}
	if (M.dims[2] == 1)
		sz = 0;
	else
	{
		if (M.coordinates[2] == 0)
			sz = 0;
		else
			sz = 1;
	}
	int x, y, z;
#pragma omp parallel for private(x, y, z)
	for (int ic = sx; ic < uc.axis()[0]; ic++)
		for (int jc = sy; jc < uc.axis()[1]; jc++)
			for (int kc = sz; kc < uc.axis()[2]; kc++)
			{
				if (M.dims[0] == 1)
					x = 2 * ic;
				else
				{
					if (M.coordinates[0] == 0)
						x = 2 * ic;
					else
						x = 2 * ic - 1;
				}
				if (M.dims[1] == 1)
					y = 2 * jc;
				else
				{
					if (M.coordinates[1] == 0)
						y = 2 * jc;
					else
						y = 2 * jc - 1;
				}
				if (M.dims[2] == 1)
					z = 2 * kc;
				else
				{
					if (M.coordinates[2] == 0)
						z = 2 * kc;
					else
						z = 2 * kc - 1;
				}

				uf(x, y, z) += (uc(ic, jc, kc) - ucold(ic, jc, kc));
				if (x > 0)
					uf(x - 1, y, z) += 0.5 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
											  (uc(ic - 1, jc, kc) - ucold(ic - 1, jc, kc)));
				if (y > 0)
					uf(x, y - 1, z) += 0.5 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
											  (uc(ic, jc - 1, kc) - ucold(ic, jc - 1, kc)));
				if (z > 0)
					uf(x, y, z - 1) += 0.5 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
											  (uc(ic, jc, kc - 1) - ucold(ic, jc, kc - 1)));

				if (x > 0 && y > 0)
					uf(x - 1, y - 1, z) += 0.25 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
												   (uc(ic, jc - 1, kc) - ucold(ic, jc - 1, kc)) +
												   (uc(ic - 1, jc, kc) - ucold(ic - 1, jc, kc)) +
												   (uc(ic - 1, jc - 1, kc) - ucold(ic - 1, jc - 1, kc)));
				if (x > 0 && z > 0)
					uf(x - 1, y, z - 1) += 0.25 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
												   (uc(ic - 1, jc, kc - 1) - ucold(ic - 1, jc, kc - 1)) +
												   (uc(ic, jc, kc - 1) - ucold(ic, jc, kc - 1)) +
												   (uc(ic - 1, jc, kc) - ucold(ic - 1, jc, kc)));
				if (y > 0 && z > 0)
					uf(x, y - 1, z - 1) += 0.25 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
												   (uc(ic, jc - 1, kc - 1) - ucold(ic, jc - 1, kc - 1)) +
												   (uc(ic, jc, kc - 1) - ucold(ic, jc, kc - 1)) +
												   (uc(ic, jc - 1, kc) - ucold(ic, jc - 1, kc)));
				if (x > 0 && y > 0 && z > 0)
					uf(x - 1, y - 1, z - 1) += 0.125 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
														(uc(ic, jc - 1, kc) - ucold(ic, jc - 1, kc)) +
														(uc(ic - 1, jc, kc) - ucold(ic - 1, jc, kc)) +
														(uc(ic - 1, jc - 1, kc) - ucold(ic - 1, jc - 1, kc)) +
														(uc(ic, jc, kc - 1) - ucold(ic, jc, kc - 1)) +
														(uc(ic - 1, jc, kc - 1) - ucold(ic - 1, jc, kc - 1)) +
														(uc(ic, jc - 1, kc - 1) - ucold(ic, jc - 1, kc - 1)) +
														(uc(ic - 1, jc - 1, kc - 1) - ucold(ic - 1, jc - 1, kc - 1)));
			}
}

//------------- prologate d
void refine_d(grid<double> &df, const grid<double> &dc, const grid<double> doldc, MPI_Setting &M)
{

	refine_standard_routine(df, dc, doldc, M);
	//----reset direchelet boundary condition
	// Dirichlet_BD_d(U,d,l,M);
	//----update gohst points
	derived_datatypes(df.axis(),&M);
	SendRecv(df, M);
	mpi_data_free(M);
}

//-------------- prolongate u
void refine_u(Stack *U, gdvecdouble &uf, const gdvecdouble &uc, const gdvecdouble &ucold, 
MPI_Setting &M, const double Ut, const double lfine)
{
	//---- standard proogation
	for (int d = 0; d<3; d++)
		refine_standard_routine(uf.getgrid(d), uc.getgrid(d), ucold.getgrid(d), M);
	//----reset direchelet boundary condition
	Dirichlet_BD(U, uf, M, lfine, Ut);
	derived_datatypes(uf.getgrid(0).axis(),&M);
	//----update gohst points
	for (int d = 0; d<3; d++)
		SendRecv(uf.getgrid(d), M);
	mpi_data_free(M);
}

/*

//------------- standard prologation operator for u
void refine_standard_routine_u(griddouble &uf, const griddouble &uc, const griddouble &ucold, 
const griddouble & df, const MPI_Setting &M)
{
	/// Interpolation and addition of coarse grid correction from grid l-1
	/// to grid l
	int sx, sy, sz;
	if (M.dims[0] == 1)
		sx = 0;
	else
	{
		if (M.coordinates[0] == 0)
			sx = 0;
		else
			sx = 1;
	}
	if (M.dims[1] == 1)
		sy = 0;
	else
	{
		if (M.coordinates[1] == 0)
			sy = 0;
		else
			sy = 1;
	}
	if (M.dims[2] == 1)
		sz = 0;
	else
	{
		if (M.coordinates[2] == 0)
			sz = 0;
		else
			sz = 1;
	}
	int x, y, z, a;
#pragma omp parallel for private(x, y, z, a)
	for (int kc = sz; kc < uc.getgrid(0).axis()[2]; kc++)
		for (int jc = sy; jc < uc.getgrid(0).axis()[1]; jc++)
			for (int ic = sx; ic < uc.getgrid(0).axis()[0]; ic++)
			{
				if (M.dims[0] == 1)
					x = 2 * ic;
				else
				{
					if (M.coordinates[0] == 0)
						x = 2 * ic;
					else
						x = 2 * ic - 1;
				}
				if (M.dims[1] == 1)
					y = 2 * jc;
				else
				{
					if (M.coordinates[1] == 0)
						y = 2 * jc;
					else
						y = 2 * jc - 1;
				}
				if (M.dims[2] == 1)
					z = 2 * kc;
				else
				{
					if (M.coordinates[2] == 0)
						z = 2 * kc;
					else
						z = 2 * kc - 1;
				}

				// avoid nodes where d>0.9
				if (d(x, y, z) <= 0.9)
				{
					uf(x, y, z) += (uc(ic, jc, kc) - ucold(ic, jc, kc));
					if (x > 0)
						uf(x - 1, y, z) += 0.5 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
												  (uc(ic - 1, jc, kc) - ucold(ic - 1, jc, kc)));
					if (y > 0)
						uf(x, y - 1, z) += 0.5 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
												  (uc(ic, jc - 1, kc) - ucold(ic, jc - 1, kc)));
					if (z > 0)
						uf(x, y, z - 1) += 0.5 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
												  (uc(ic, jc, kc - 1) - ucold(ic, jc, kc - 1)));

					if (x > 0 && y > 0)
						uf(x - 1, y - 1, z) += 0.25 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
													   (uc(ic, jc - 1, kc) - ucold(ic, jc - 1, kc)) +
													   (uc(ic - 1, jc, kc) - ucold(ic - 1, jc, kc)) +
													   (uc(ic - 1, jc - 1, kc) - ucold(ic - 1, jc - 1, kc)));
					if (x > 0 && z > 0)
						uf(x - 1, y, z - 1) += 0.25 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
													   (uc(ic - 1, jc, kc - 1) - ucold(ic - 1, jc, kc - 1)) +
													   (uc(ic, jc, kc - 1) - ucold(ic, jc, kc - 1)) +
													   (uc(ic - 1, jc, kc) - ucold(ic - 1, jc, kc)));
					if (y > 0 && z > 0)
						uf(x, y - 1, z - 1) += 0.25 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
													   (uc(ic, jc - 1, kc - 1) - ucold(ic, jc - 1, kc - 1)) +
													   (uc(ic, jc, kc - 1) - ucold(ic, jc, kc - 1)) +
													   (uc(ic, jc - 1, kc) - ucold(ic, jc - 1, kc)));
					if (x > 0 && y > 0 && z > 0)
						uf(x - 1, y - 1, z - 1) += 0.125 * ((uc(ic, jc, kc) - ucold(ic, jc, kc)) +
															(uc(ic, jc - 1, kc) - ucold(ic, jc - 1, kc)) +
															(uc(ic - 1, jc, kc) - ucold(ic - 1, jc, kc)) +
															(uc(ic - 1, jc - 1, kc) - ucold(ic - 1, jc - 1, kc)) +
															(uc(ic, jc, kc - 1) - ucold(ic, jc, kc - 1)) +
															(uc(ic - 1, jc, kc - 1) - ucold(ic - 1, jc, kc - 1)) +
															(uc(ic, jc - 1, kc - 1) - ucold(ic, jc - 1, kc - 1)) +
															(uc(ic - 1, jc - 1, kc - 1) - ucold(ic - 1, jc - 1, kc - 1)));
				}
			}
}
//--------------- Alcouff type prolongation operator
void refine_alcouff_routine(Stack *U, multigrid3d &u, multigrid3d &uold,
multigrid3d &MP, int l,MPI_Setting & M)
{
/// Interpolation and addition of coarse grid correction from grid l-1
/// to grid l

Level  *Lc;
Lc =U->Ll+l-1;
double juck=1.0e-6; // to avoid 0 in denominator
int sx,sy,sz;
if(M.dims[0]==1) sx =0;
else{
	if (M.coordinates[0]==0) sx =0;
	else sx =1;}
if(M.dims[1]==1) sy =0;
else{
	if (M.coordinates[1]==0) sy =0;
	else sy =1;}
if(M.dims[2]==1) sz =0;
else{
	if (M.coordinates[2]==0) sz =0;
	else sz =1;}
int x,y,z;
#pragma omp parallel for private(x,y,z)
for (int ic=sx;ic<=Lc->nbcolox-1;ic++)
  for (int jc=sy;jc<=Lc->nbcoloy-1;jc++)
	  for (int kc=sz;kc<=Lc->nbcoloz-1;kc++)
	  {
		if(M.dims[0]==1) x = 2*ic;
		else{
			if (M.coordinates[0]==0) x = 2*ic;
			else x = 2*ic-1;}
		if(M.dims[1]==1) y = 2*jc;
		else{
			if (M.coordinates[1]==0) y = 2*jc;
			else y = 2*jc-1;}
		if(M.dims[2]==1) z = 2*kc;
		else{
			if (M.coordinates[2]==0) z = 2*kc;
			else z = 2*kc-1;}

		u(l-1, x  , y  , z  )+=(u(l-2, ic, jc, kc)-uold(l-2, ic, jc, kc));
		if(x>0)
		u(l-1, x-1, y  , z  )+=(MP(l-2, ic  , jc, kc)*(u(l-2, ic  , jc, kc)-uold(l-2, ic  , jc, kc))+
								MP(l-2, ic-1, jc, kc)*(u(l-2, ic-1, jc, kc)-uold(l-2, ic-1, jc, kc)))/
								(MP(l-2, ic  , jc, kc)+MP(l-2, ic-1, jc, kc));
		if(y>0)
		u(l-1, x  , y-1, z  )+=(MP(l-2, ic, jc  , kc)*(u(l-2, ic, jc  , kc)-uold(l-2, ic, jc  , kc))+
								MP(l-2, ic, jc-1, kc)*(u(l-2, ic, jc-1, kc)-uold(l-2, ic, jc-1, kc)))/
								(MP(l-2, ic  , jc, kc)+MP(l-2, ic, jc-1, kc));
		if(z>0)
		u(l-1, x  , y  , z-1)+=(MP(l-2, ic, jc, kc)*(u(l-2, ic, jc  , kc)-uold(l-2, ic, jc  , kc))+
								MP(l-2, ic, jc, kc-1)*(u(l-2, ic, jc, kc-1)-uold(l-2, ic, jc, kc-1)))/
								(MP(l-2, ic  , jc, kc)+MP(l-2, ic, jc, kc-1));

		if(x>0&&y>0)
		u(l-1, x-1, y-1, z  )+=(MP(l-2, ic  , jc  , kc)*(u(l-2, ic  , jc  , kc)-uold(l-2, ic  , jc  , kc))+
								MP(l-2, ic  , jc-1, kc)*(u(l-2, ic  , jc-1, kc)-uold(l-2, ic  , jc-1, kc))+
								MP(l-2, ic-1, jc  , kc)*(u(l-2, ic-1, jc  , kc)-uold(l-2, ic-1, jc  , kc))+
								MP(l-2, ic-1, jc-1, kc)*(u(l-2, ic-1, jc-1, kc)-uold(l-2, ic-1, jc-1, kc)))/
								(MP(l-2, ic  , jc  , kc)+MP(l-2, ic  , jc-1, kc)+
								 MP(l-2, ic-1, jc  , kc)+MP(l-2, ic-1, jc-1, kc));
		if(x>0&&z>0)
		u(l-1, x-1, y  , z-1)+=(MP(l-2, ic  , jc, kc  )*(u(l-2, ic  , jc, kc  )-uold(l-2, ic  , jc, kc  ))+
								MP(l-2, ic-1, jc, kc-1)*(u(l-2, ic-1, jc, kc-1)-uold(l-2, ic-1, jc, kc-1))+
								MP(l-2, ic  , jc, kc-1)*(u(l-2, ic  , jc, kc-1)-uold(l-2, ic  , jc, kc-1))+
								MP(l-2, ic-1, jc, kc  )*(u(l-2, ic-1, jc, kc  )-uold(l-2, ic-1, jc, kc  )))/
								(MP(l-2, ic  , jc  , kc)+MP(l-2, ic  , jc, kc-1)+
								 MP(l-2, ic-1, jc  , kc)+MP(l-2, ic-1, jc, kc-1));
		if(y>0&&z>0)
		u(l-1, x  , y-1, z-1)+=(MP(l-2, ic, jc  , kc  )*(u(l-2, ic, jc  , kc  )-uold(l-2, ic, jc  , kc  ))+
								MP(l-2, ic, jc-1, kc-1)*(u(l-2, ic, jc-1, kc-1)-uold(l-2, ic, jc-1, kc-1))+
								MP(l-2, ic, jc  , kc-1)*(u(l-2, ic, jc  , kc-1)-uold(l-2, ic, jc  , kc-1))+
								MP(l-2, ic, jc-1, kc  )*(u(l-2, ic, jc-1, kc  )-uold(l-2, ic, jc-1, kc  )))/
								(MP(l-2, ic, jc  , kc  )+MP(l-2, ic, jc-1, kc-1)+
								 MP(l-2, ic, jc  , kc-1)+MP(l-2, ic, jc-1, kc  ));
		if(x>0&&y>0&&z>0)
		u(l-1, x-1, y-1, z-1)+=(MP(l-2, ic  , jc  , kc  )*(u(l-2, ic  , jc  , kc  )-uold(l-2, ic  , jc  , kc  ))+
								MP(l-2, ic  , jc-1, kc  )*(u(l-2, ic  , jc-1, kc  )-uold(l-2, ic  , jc-1, kc  ))+
								MP(l-2, ic-1, jc  , kc  )*(u(l-2, ic-1, jc  , kc  )-uold(l-2, ic-1, jc  , kc  ))+
								MP(l-2, ic-1, jc-1, kc  )*(u(l-2, ic-1, jc-1, kc  )-uold(l-2, ic-1, jc-1, kc  ))+
								MP(l-2, ic  , jc  , kc-1)*(u(l-2, ic  , jc  , kc-1)-uold(l-2, ic  , jc  , kc-1))+
								MP(l-2, ic-1, jc  , kc-1)*(u(l-2, ic-1, jc  , kc-1)-uold(l-2, ic-1, jc  , kc-1))+
								MP(l-2, ic  , jc-1, kc-1)*(u(l-2, ic  , jc-1, kc-1)-uold(l-2, ic  , jc-1, kc-1))+
								MP(l-2, ic-1, jc-1, kc-1)*(u(l-2, ic-1, jc-1, kc-1)-uold(l-2, ic-1, jc-1, kc-1)))/
								(MP(l-2, ic  , jc  , kc  )+MP(l-2, ic  , jc-1, kc-1)+
								 MP(l-2, ic  , jc  , kc-1)+MP(l-2, ic  , jc-1, kc  )+
								 MP(l-2, ic-1, jc  , kc  )+MP(l-2, ic-1, jc  , kc-1)+
								 MP(l-2, ic-1, jc-1, kc  )+MP(l-2, ic-1, jc-1, kc-1));
	  }
}
*/
/*
// fmg interpolation to obtain initial solution on grid l+1
void fmg_interpolate_routine(Stack *U, multigrid3d &u, multigrid3d &uconv,
							 int l, MPI_Setting &M)
{
	/// interpolation of coarse grid k-1 solution to fine grid l
	/// to serve as first approximation. bi-cubic interpolation

	Level *Lc;
	Lc = U->Ll + l - 1;

	/// store grid l-1 solution for later use in convergence check
	for (int ic = 0; ic <= Lc->nbcolox - 1; ic++)
		for (int jc = 0; jc <= Lc->nbcoloy - 1; jc++)
			for (int kc = 0; kc <= Lc->nbcoloz - 1; kc++)
				uconv(l - 2, ic, jc, kc) = u(l - 2, ic, jc, kc);

	/// first inject to points coinciding with coarse grid points
	int x, y, z;
	int sx, sy, sz;
	if (M.dims[0] == 1)
		sx = 0;
	else
	{
		if (M.coordinates[0] == 0)
			sx = 0;
		else
			sx = 1;
	}
	if (M.dims[1] == 1)
		sy = 0;
	else
	{
		if (M.coordinates[1] == 0)
			sy = 0;
		else
			sy = 1;
	}
	if (M.dims[2] == 1)
		sz = 0;
	else
	{
		if (M.coordinates[2] == 0)
			sz = 0;
		else
			sz = 1;
	}
	for (int ic = sx; ic <= Lc->nbcolox - 1; ic++)
		for (int jc = sy; jc <= Lc->nbcoloy - 1; jc++)
			for (int kc = sz; kc <= Lc->nbcoloz - 1; kc++)
			{
				if (M.dims[0] == 1)
					x = 2 * ic;
				else
				{
					if (M.coordinates[0] == 0)
						x = 2 * ic;
					else
						x = 2 * ic - 1;
				}
				if (M.dims[1] == 1)
					y = 2 * jc;
				else
				{
					if (M.coordinates[1] == 0)
						y = 2 * jc;
					else
						y = 2 * jc - 1;
				}
				if (M.dims[2] == 1)
					z = 2 * kc;
				else
				{
					if (M.coordinates[2] == 0)
						z = 2 * kc;
					else
						z = 2 * kc - 1;
				}
				u(l - 1, x, y, z) = u(l - 2, ic, jc, kc);
				if (x > 1)
				{
					u(l - 1, x - 1, y, z) = (u(l - 2, ic, jc, kc) + u(l - 2, ic - 1, jc, kc)) * 0.5;
				}
				if (y > 1)
				{
					u(l - 1, x, y - 1, z) = (u(l - 2, ic, jc, kc) + u(l - 2, ic, jc - 1, kc)) * 0.5;
				}
				if (z > 1)
				{
					u(l - 1, x, y, z - 1) = (u(l - 2, ic, jc, kc) + u(l - 2, ic, jc, kc - 1)) * 0.5;
				}
				if (x > 1 && y > 1)
				{
					u(l - 1, x - 1, y - 1, z) = (u(l - 2, ic, jc, kc) + u(l - 2, ic, jc - 1, kc) +
												 u(l - 2, ic - 1, jc, kc) + u(l - 2, ic - 1, jc - 1, kc)) *
												0.25;
				}
				if (x > 1 && z > 1)
				{
					u(l - 1, x - 1, y, z - 1) = (u(l - 2, ic, jc, kc) + u(l - 2, ic - 1, jc, kc - 1) +
												 u(l - 2, ic, jc, kc - 1) + u(l - 2, ic - 1, jc, kc)) *
												0.25;
				}
				if (y > 1 && z > 1)
				{
					u(l - 1, x, y - 1, z - 1) = (u(l - 2, ic, jc, kc) + u(l - 2, ic, jc - 1, kc - 1) +
												 u(l - 2, ic, jc, kc - 1) + u(l - 2, ic, jc - 1, kc)) *
												0.25;
				}
				if (x > 1 && y > 1 && z > 1)
				{
					u(l - 1, x - 1, y - 1, z - 1) = (u(l - 2, ic, jc, kc) + u(l - 2, ic, jc - 1, kc) +
													 u(l - 2, ic - 1, jc, kc) + u(l - 2, ic - 1, jc - 1, kc) +
													 u(l - 2, ic, jc, kc - 1) + u(l - 2, ic - 1, jc, kc - 1) +
													 u(l - 2, ic, jc - 1, kc - 1) + u(l - 2, ic - 1, jc - 1, kc - 1)) *
													0.125;
				}
			}
}
*/
/*
//------------------ fmg interpolation for d
void fmg_interpolation_d(Stack *U, multigrid3d &d, multigrid3d &dconv,
						int l,MPI_Setting & M)
{
	Level  *L;
	L =U->Ll+l;
	fmg_interpolate_routine(U,d,dconv,l,M);
	// Dirichlet_BD_d(U,d,l,M);
	SendRecv(d,l-1,L->nbcolox,L->nbcoloy,L->nbcoloz,M);
}
*/
/*
//------------------ fmg interpolation for u
void fmg_interpolation_u(Stack *U, multigrid3d &u, multigrid3d &uconv, multigrid3d &v,
						 multigrid3d &vconv, multigrid3d &w, multigrid3d &wconv,
						 int l, MPI_Setting &M, double Ut)
{
	Level *L;
	L = U->Ll + l;
	fmg_interpolate_routine(U, u, uconv, l, M);
	Dirichlet_BD(U, u, l, M, 0, Ut);
	SendRecv(u, l - 1, L->nbcolox, L->nbcoloy, L->nbcoloz, M);
	fmg_interpolate_routine(U, v, vconv, l, M);
	Dirichlet_BD(U, v, l, M, 1, Ut);
	SendRecv(v, l - 1, L->nbcolox, L->nbcoloy, L->nbcoloz, M);
	fmg_interpolate_routine(U, w, wconv, l, M);
	Dirichlet_BD(U, w, l, M, 2, Ut);
	SendRecv(w, l - 1, L->nbcolox, L->nbcoloy, L->nbcoloz, M);
}
*/
/*-------------convergence check (L2 norm) using converged solutions
			   on level k and on next coarser grid k-1
double conver(Stack *U, multigrid3d &u, multigrid3d &uconv,
			  int l, MPI_Setting &M)
{

	Level *Lc, *L;
	double err = 0.0, serr = 0.0;

	L = U->Ll + l;
	Lc = U->Ll + l - 1;
	int sx, sy, sz, ex, ey, ez;
	if (M.dims[0] == 1)
	{
		sx = 0;
		ex = Lc->nbcolox - 1;
	}
	else
	{
		if (M.coordinates[0] == 0)
			sx = 0;
		else
			sx = 1;
		if (M.coordinates[0] == M.dims[0] - 1)
			ex = Lc->nbcolox - 1;
		else
			ex = Lc->nbcolox - 2;
	}
	if (M.dims[1] == 1)
	{
		sy = 0;
		ey = Lc->nbcoloy - 1;
	}
	else
	{
		if (M.coordinates[1] == 0)
			sy = 0;
		else
			sy = 1;
		if (M.coordinates[1] == M.dims[1] - 1)
			ey = Lc->nbcoloy - 1;
		else
			ey = Lc->nbcoloy - 2;
	}
	if (M.dims[2] == 1)
	{
		sz = 0;
		ez = Lc->nbcoloz - 1;
	}
	else
	{
		if (M.coordinates[2] == 0)
			sz = 0;
		else
			sz = 1;
		if (M.coordinates[2] == M.dims[2] - 1)
			ez = Lc->nbcoloz - 1;
		else
			ez = Lc->nbcoloz - 2;
	}

	int x, y, z;
	if (l == U->maxlevel)
	{
#pragma omp parallel for private(x, y, z) reduction(+ \
													: err)
		for (int ic = sx; ic <= ex; ic++)
			for (int jc = sy; jc <= ey; jc++)
				for (int kc = sz; kc <= ez; kc++)
				{
					if (M.dims[0] == 1)
						x = 2 * ic;
					else
					{
						if (M.coordinates[0] == 0)
							x = 2 * ic;
						else
							x = 2 * ic - 1;
					}
					if (M.dims[1] == 1)
						y = 2 * jc;
					else
					{
						if (M.coordinates[1] == 0)
							y = 2 * jc;
						else
							y = 2 * jc - 1;
					}
					if (M.dims[2] == 1)
						z = 2 * kc;
					else
					{
						if (M.coordinates[2] == 0)
							z = 2 * kc;
						else
							z = 2 * kc - 1;
					}
					err += ((uconv(l - 2, ic, jc, kc) - u(l - 1, x, y, z)) * (uconv(l - 2, ic, jc, kc) - u(l - 1, x, y, z)));
				}
	}
	else
	{
#pragma omp parallel for private(x, y, z) reduction(+ \
													: err)
		for (int ic = sx; ic <= ex; ic++)
			for (int jc = sy; jc <= ey; jc++)
				for (int kc = sz; kc <= ez; kc++)
				{
					if (M.dims[0] == 1)
						x = 2 * ic;
					else
					{
						if (M.coordinates[0] == 0)
							x = 2 * ic;
						else
							x = 2 * ic - 1;
					}
					if (M.dims[1] == 1)
						y = 2 * jc;
					else
					{
						if (M.coordinates[1] == 0)
							y = 2 * jc;
						else
							y = 2 * jc - 1;
					}
					if (M.dims[2] == 1)
						z = 2 * kc;
					else
					{
						if (M.coordinates[2] == 0)
							z = 2 * kc;
						else
							z = 2 * kc - 1;
					}
					err += ((uconv(l - 2, ic, jc, kc) - uconv(l - 1, x, y, z)) * (uconv(l - 2, ic, jc, kc) - uconv(l - 1, x, y, z)));
				}
	}
	MPI_Reduce(&err, &serr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&serr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return (serr / ((Lc->ii) * (Lc->jj) * (Lc->kk)));
}
*/