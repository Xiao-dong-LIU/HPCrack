/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "intergrid.h"
#include "mg_cycle_u.h"
#include "pcg_u.h"
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
	
		refine_u(U,u.getLevel(l),u.getLevel(l-1),uold.getLevel(l-1),M,Ut,l);

		pcg_u(U,u.getLevel(l),ru.getLevel(l),fu.getLevel(l),MP_choice_d.getLevel(l),bulK.getLevel(l),
	G.getLevel(l),d.getLevel(l),K_element.getLevel(l),G_element.getLevel(l),M,l,small_k,Ut,
	mgp_u.nu2,t,f_norm);
	  }
}