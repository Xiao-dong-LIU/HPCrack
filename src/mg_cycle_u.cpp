/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "intergrid.h"
#include "mg_cycle_u.h"
#include "PCG_u.h"
#include "output_result.h"



///------------MULTIGRID DRIVING ROUTINES------------///

void cycle_u(Stack *U, mgvecdouble &u, mgvecdouble & ru, mgvecdouble & fu, mgvecdouble & uold, 
const mgshort & MP_choice_d, const mgdouble & bulK,  const mgdouble & G,  const mgdouble & d,  
const mgdouble & K_element, const mgdouble & G_element, int l, const MG & mgp_u,  MPI_Setting & M, 
const int myid, const int nbprocs, const double small_k, const double Ut, const int t, double & f_norm)
{
/// performs coarse grid correction cycle starting on level l
/// nu1 pre-relaxations, nu2 postrelaxations, nu0 relaxations
/// on the coarsest grid, cycleindex gamma=1 for Vcycle,
/// gamma=2 for Wcycle


	if (l==0)
	pcg_u(U,u.getLevel(l),ru.getLevel(l),fu.getLevel(l),MP_choice_d.getLevel(l),bulK.getLevel(l),
	G.getLevel(l),d.getLevel(l),K_element.getLevel(l),G_element.getLevel(l),M,l,small_k,Ut,
	mgp_u.nu0,t,f_norm);
	else
	  {
		pcg_u(U,u.getLevel(l),ru.getLevel(l),fu.getLevel(l),MP_choice_d.getLevel(l),bulK.getLevel(l),
	G.getLevel(l),d.getLevel(l),K_element.getLevel(l),G_element.getLevel(l),M,l,small_k,Ut,
	mgp_u.nu1,t,f_norm);
		coarsen_u(U,u.getLevel(l-1),uold.getLevel(l-1),u.getLevel(l),bulK.getLevel(l),M,l-1,Ut);
		
		coarsen_fu(U,fu.getLevel(l-1),u.getLevel(l-1),MP_choice_d.getLevel(l-1),
		ru.getLevel(l),bulK.getLevel(l-1),G.getLevel(l-1),d.getLevel(l-1),K_element.getLevel(l-1),
		G_element.getLevel(l-1),M,l,small_k);

		for (int j=1;j<=mgp_u.gamma;j++) 
			cycle_u(U,u,ru,fu,uold,MP_choice_d,bulK,G,d,K_element,
			G_element,l-1,mgp_u,M,myid,nbprocs,small_k,Ut,t,f_norm);
		  /*
		  char filename_u[1024];
       	  sprintf(filename_u, "U_T_%d", l);
          write_total_u(U,u,v,w,M,l,myid,nbprocs,filename_u);
			*/
		refine_u(U,u.getLevel(l),u.getLevel(l-1),uold.getLevel(l-1),M,Ut,l);
/*
		 ///------------ output r ------------///
		 char filename_r[1024];
		 char field_name_r[100]="r";    // output field name in the VTK file
		 sprintf(filename_r, "%s", field_name_r);
		 write_total_u(U,ru,rv,rw,M,l,myid,nbprocs,filename_r,field_name_r,t,1.0);
*/
		pcg_u(U,u.getLevel(l),ru.getLevel(l),fu.getLevel(l),MP_choice_d.getLevel(l),bulK.getLevel(l),
	G.getLevel(l),d.getLevel(l),K_element.getLevel(l),G_element.getLevel(l),M,l,small_k,Ut,
	mgp_u.nu2,t,f_norm);
	  }
}
/*
void fmg_u(Stack *U, multigrid3d &u, multigrid3d &v, multigrid3d &w, multiinput3d& MP_choice_d,
multigrid3d &fu, multigrid3d &fv, multigrid3d &fw, multigrid3d &uold, multigrid3d &vold, 
multigrid3d &wold,multigrid3d &uconv, multigrid3d &vconv, multigrid3d &wconv, multigrid3d & bulK, 
multigrid3d & G, multigrid3d & d, multigrid3d & K_element, multigrid3d & G_element, multigrid3d &ru, 
multigrid3d &rv, multigrid3d &rw, int l, const MG & mgp_u, MPI_Setting & M, int myid, int nbprocs, 
double small_k, double Ut, int t, double & f_norm)
{
/// performs FMG with k levels and ncy cycles per level
	if (U->maxlevel==1)
	  for (int j=1;j< mgp_u.ncy;j++)
	  	precondition_conjugate_gradient_u(U,u,v,w,MP_choice_d,ru,rv,rw,fu,fv,fw,bulK,G,d,
		  K_element,G_element,l,M,small_k,Ut,mgp_u.nu0,t,f_norm);
	
	if (l==1) 
		precondition_conjugate_gradient_u(U,u,v,w,MP_choice_d,ru,rv,rw,fu,fv,fw,bulK,G,d,
		K_element,G_element,l,M,small_k,Ut,mgp_u.nu0,t,f_norm);
	else
	  {
		  fmg_u(U,u,v,w,MP_choice_d,fu,fv,fw,uold,vold,wold,uconv,vconv,wconv,bulK,G,d,
		  K_element,G_element,ru,rv,rw,l-1,mgp_u,M,myid,nbprocs,small_k,Ut,t,f_norm);
		  fmg_interpolation_u(U,u,uconv,v,vconv,w,wconv,l,M,Ut);
		  for (int j=1;j<=mgp_u.ncy;j++)
			{
			cycle_u(U,u,v,w,MP_choice_d,fu,fv,fw,uold,vold,wold,bulK,G,d,K_element,
			G_element,ru,rv,rw,l,mgp_u,M,myid,nbprocs,small_k,Ut,t,f_norm);
			if(M.coordinates[0]==0&&M.coordinates[1]==0&&M.coordinates[2]==0)printf("\n");
			}
	  }
}
*/