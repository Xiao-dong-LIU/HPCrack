/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "PF_schema.h"
#include "u_total.h"
#include "d_total.h"
#include "strain_history.h"
#include "initialisation.h"

#include "output_result.h"
#include "imput_image.h"
#include "Fd.h"
#include "BD.h"
#include "element_coarse.h"
#include "write_H.h"
#include <iomanip>
#include "externe_forcce.h"
#include "d_sum.h"
#include <fstream>
#include "mgvec.h"

//------------ the phase field method 

void phm(Stack *U, const mgdouble & bulK, const mgdouble & G, const mgdouble & gc, MPI_Setting & M, const Parameters & para,
const int myid, const int nbprocs)
{
     
    double small_k =1.e-5;
    double fd_norm=0., f_norm=0.;
    double Ue=0.0,Ut=0.0;
    double delta_U;
    double u_proportion,U_unit;
    U_unit=1.;
    double FeW;
    double d_max=1.;
    mg<double> d(U,U->maxlevel);
    mg<double> H(U,U->maxlevel);
    mg<double> fd(U,U->maxlevel);
    mgvecdouble u(U,U->maxlevel,3);
    mgvecdouble fu(U,U->maxlevel,3);
    grid<double> d_old(d.getLevel(U->maxlevel));
/*
    read_vector_3D(u.getLevel(U->maxlevel).getgrid(0),"u",myid,0); 
    read_vector_3D(u.getLevel(U->maxlevel).getgrid(1),"v",myid,0); 
    read_vector_3D(u.getLevel(U->maxlevel).getgrid(2),"w",myid,0); 
*/
    std::ofstream F_U ("F_U.dat");
    if (F_U.is_open())
    {
        if (myid==0)
        std::cout<<"File opened!"<<std::endl;
    }
    else
    {
        if (myid==0)
        std::cout<<"Can not open file!"<<std::endl;
        exit(-1);
    }

	double t1, t2;
    for (int t=0;t<5000;t++)
    {
        if (myid==0)
        cout<<"***** At Time step "<<t<<endl;
        // ----- compute u
            u_entire(U,u,fu,bulK,G,d,M,para.mg_u,myid,nbprocs,small_k,U_unit,t,f_norm);   
            if (t==0)
            {
                Ue=U_e(U,u.getLevel(U->maxlevel),bulK.getLevel(U->maxlevel),G.getLevel(U->maxlevel),
                gc.getLevel(U->maxlevel),U->maxlevel,para.lc,M);
                Ut=Ue;
                delta_U=1.5e-3*Ue;
                if(myid==0)
                std::cout << "Ue="<< Ue << std::endl;
                u_proportion=Ut/U_unit;
            }
            // ----- compute H
            strain_history(U,H.getLevel(U->maxlevel),u.getLevel(U->maxlevel),bulK.getLevel(U->maxlevel),
            G.getLevel(U->maxlevel),gc.getLevel(U->maxlevel),M,U->maxlevel,para.lc,u_proportion);
            Voigt_all_level(H,M);
            // ----- compute right hand of d
        for (int l=U->maxlevel;l>=0;l--)
        {
             fd_h(U,fd.getLevel(l),H.getLevel(l),l);
        }
        // ----- compute d
        d_entire(U,d,fd,H,gc,M,para.lc,para.mg_d,myid,nbprocs,t,fd_norm);
        d_max = delta_d_max(d.getLevel(U->maxlevel),d_old);
        d_old = d.getLevel(U->maxlevel);
            
        FeW=externe_force_sum(fu.getLevel(U->maxlevel).getgrid(2),M);       
	    FeW = para.voxel_size*FeW;	
        if (t>=500)
        { 
                ///------------ wrtite history ------------///
                write_vector_3D(u.getLevel(U->maxlevel).getgrid(0),"u",myid,t); 
                write_vector_3D(u.getLevel(U->maxlevel).getgrid(1),"v",myid,t); 
                write_vector_3D(u.getLevel(U->maxlevel).getgrid(2),"w",myid,t); 
                write_vector_3D(H.getLevel(U->maxlevel),"H",myid,t); 
                write_vector_3D(d.getLevel(U->maxlevel),"d",myid,t); 
                ///------------ output d ------------///
                write_total_de(U,d.getLevel(U->maxlevel),M,"S","d","d",U->maxlevel,myid,nbprocs,t);
        }
        Ut=Ut+delta_U;
        u_proportion=Ut/U_unit;
        if (myid==0)
        {
	      F_U<<"T= "<<t<<",U= "<<std::setprecision(10)<<Ut<<" ,Fw= "<<std::setprecision(10)
          <<u_proportion*FeW<<" ,deltaU= "<<delta_U<<" ,d_max= "<<d_max<<std::endl;
	    }
        if (u_proportion*FeW<0.5)
	    break;
    }
    F_U.close();
}
