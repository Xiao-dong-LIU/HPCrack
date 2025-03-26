/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "material_selection.h"
#include "mechanical_displacement.h"
#include "mechanical_multigrid.h"
#include "mechanical_pcg.h"
#include "multigrid_intergrid.h"
#include "vtk_output.h"
#include "image_input.h"
#include "pcg_utilities.h"
#include "mechanical_pcg.h"
#include "boundary_conditions.h"
#include "mechanical_right_hand_side.h"
#include "material_coarse_element.h"

// ---------- the entire progrm to calculate u

void u_entire(Stack *U, mgvecdouble &u, mgvecdouble &fu, const mgdouble &bulK, const mgdouble &G, 
const mgdouble &d, MPI_Setting & M, const MG & mgp_u, const int myid, const int nbprocs, 
const double small_k, const double Ut,const int t,double & f_norm)
{
    if(myid==0)
      std::cout<<" For u "<<std::endl;
    // define vetors for computing
    mgvecdouble uold(U,U->maxlevel-1,3);
    mgvecdouble uconv(U,U->maxlevel-1,3);
    
    mgvecdouble ru(U,U->maxlevel,3);

    mgshort MP_choice_d(U,U->maxlevel);
    mgdouble K_element(U,U->maxlevel);
	  mgdouble G_element(U,U->maxlevel);
    // creat a field to control compression and traction case, 0 for tr(e)>=0, 1 for tr(e)<0, 
    // 2 for phi+<phi-
    MP_choice(U,MP_choice_d,u,bulK,G,d);
    // compute material property in element 
    MP_all_element(U,u.getLevel(U->maxlevel),bulK.getLevel(U->maxlevel),G.getLevel(U->maxlevel),
    d.getLevel(U->maxlevel),MP_choice_d.getLevel(U->maxlevel),K_element.getLevel(U->maxlevel),
    G_element.getLevel(U->maxlevel),U->maxlevel,small_k);
    // compute material property on all levels 
    VR_all_level(K_element,M);
    VR_all_level(G_element,M);
/*
    ///------------ output K_element ------------///
    char filename_K[1024];
    char field_name_K[100];    // output field name in the VTK file
    sprintf(field_name_K, "K_before_element");
    sprintf(filename_K, "%s", field_name_K);
    std::string output_option_K="C"; // output celldata
    write_total_H(U,K_element,M,output_option_K,filename_K,field_name_K,U->maxlevel,myid,nbprocs,t); 
    */
 //  account_nb_zeros(U,K_element,M);
	/*
    if (t==0)
    {
      // fmg cycles for u
      fmg_u(U,u,v,w,MP_choice_d,fu,fv,fw,uold,vold,wold,uconv,vconv,wconv,bulK,G,d,
      K_element,G_element, ru,rv,rw,U->maxlevel,mgp_u,M,myid,nbprocs,small_k,Ut,t,f_norm);
      Converged_Internal_Force(U,u,v,w,fu,fv,fw,MP_choice_d,bulK,G,d,K_element,G_element,Ut,U->maxlevel,small_k,M);

    }
    else 
    {
      */
    // V cycles for u
    double r_norm=0.0; 
    double relatif_residual=1.0; 

    r_norm=Initial_residual_u(U,ru.getLevel(U->maxlevel),u.getLevel(U->maxlevel),
    fu.getLevel(U->maxlevel),MP_choice_d.getLevel(U->maxlevel), bulK.getLevel(U->maxlevel),
    G.getLevel(U->maxlevel),d.getLevel(U->maxlevel),K_element.getLevel(U->maxlevel),
    G_element.getLevel(U->maxlevel),M,U->maxlevel,small_k);

    if(myid==0)
        std::cout<<" Initial Residual = "<<r_norm<<std::endl;
     //   exit(-1);
	  if (t>=1&&f_norm!=0)
      relatif_residual=r_norm/f_norm;
     for (int j=1;j<=mgp_u.ncy;j++)
         {
            if(relatif_residual>1.e-6)
            {
              cycle_u(U,u,ru,fu,uold,MP_choice_d,bulK,G,d,K_element,G_element,
              U->maxlevel,mgp_u,M,myid,nbprocs,small_k,Ut,t,f_norm);
              if (myid==0) std::cout<<std::endl;
              Dirichlet_residual_u(U,ru.getLevel(U->maxlevel),M,U->maxlevel);
		          r_norm = 0.;
              for (int di=0; di<3; di++)
				        r_norm += grid_norm(ru.getLevel(U->maxlevel).getgrid(di),M);

              relatif_residual=r_norm/f_norm;
              if(myid==0)
              std::cout<<"Ru="<<r_norm<<",f="<<f_norm<<",Cycle="<<j<<std::endl;
            }
            if (relatif_residual<=1.e-6) 
            {
              if(myid==0)
                std::cout<<" Converged "<<std::endl;
              break;
            }
         }
        externe_force_nodal(U, fu.getLevel(U->maxlevel), u.getLevel(U->maxlevel),
        MP_choice_d.getLevel(U->maxlevel), bulK.getLevel(U->maxlevel), G.getLevel(U->maxlevel),
        d.getLevel(U->maxlevel), K_element.getLevel(U->maxlevel), G_element.getLevel(U->maxlevel),
        M,U->maxlevel,Ut,small_k);


}

