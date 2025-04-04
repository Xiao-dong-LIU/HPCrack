/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "phase_field_multigrid.h"
#include "multigrid_intergrid.h"
#include "phase_field_pcg.h"



///------------MULTIGRID DRIVING ROUTINES------------///

void cycle_d(Stack *U, mg<double> &d, mg<double> &fd, mg<double> &dold, const mg<double> &H, 
const mg<double> &gc, mg<double> &r, int l, const double lc, const MG & mgp_d, 
MPI_Setting & M, const int myid, const int nbprocs, const int t, double & fd_norm)
{
/// performs coarse grid correction cycle starting on level l
/// nu1 pre-relaxations, nu2 postrelaxations, nu0 relaxations
/// on the coarsest grid, cycleindex gamma=1 for Vcycle,
/// gamma=2 for Wcycle


	if (l==0) /// NB of relaxations on the coarsest grid
		pcg_d(U,d.getLevel(l),fd.getLevel(l),r.getLevel(l),H.getLevel(l),gc.getLevel(l),M,l,lc,mgp_d.nu0,t,fd_norm);
	  //	precondition_conjugate_gradient_d(U,d,fd,r,H,gc,l,lc,M,mgp_d.nu0,t,fd_norm);
	else
	  {
		  pcg_d(U,d.getLevel(l),fd.getLevel(l),r.getLevel(l),H.getLevel(l),gc.getLevel(l),M,l,lc,mgp_d.nu1,t,fd_norm);
		  coarsen_d(d.getLevel(l-1),dold.getLevel(l-1),d.getLevel(l),M);
		  coarsen_f_d(U,fd.getLevel(l-1),d.getLevel(l-1),r.getLevel(l),gc.getLevel(l-1),H.getLevel(l-1),M,lc,l);
		  for (int j=1;j<=mgp_d.gamma;j++) 
		  {
			  cycle_d(U,d,fd,dold,H,gc,r,l-1,lc,mgp_d,M,myid,nbprocs,t,fd_norm);
		  }
		  refine_d(d.getLevel(l),d.getLevel(l-1),dold.getLevel(l-1),M);
		  for (int i=1;i<=mgp_d.nu2;i++)
		  pcg_d(U,d.getLevel(l),fd.getLevel(l),r.getLevel(l),H.getLevel(l),gc.getLevel(l),M,l,lc,mgp_d.nu2,t,fd_norm);
	  }
}