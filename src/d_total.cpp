/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "d_total.h"
#include "mg_cycle_d.h"
#include "pcg_routines.h"
#include "pcg_d.h"

void d_entire(Stack *U, mg<double> &d, mg<double> &fd, const mg<double> & H, const mg<double> &gc, 
MPI_Setting & M, double lc, const MG & mgp_d, int myid, int nbprocs,int t, double & fd_norm)
{
    if(myid==0)
      std::cout<<" For d "<<std::endl;
    mg<double> dold(U,U->maxlevel-1);
    mg<double> R(U,U->maxlevel);


    double r_norm=0.0;
    r_norm=initial_residual(U,R.getLevel(U->maxlevel),d.getLevel(U->maxlevel),
    fd.getLevel(U->maxlevel),H.getLevel(U->maxlevel),gc.getLevel(U->maxlevel),U->maxlevel,lc,M);
    double sresidual=1.0;
      for (int j=1;j<=mgp_d.ncy;j++)
      {      
        if (r_norm>1.e-12)
        {
          if(sresidual>1.e-12)
          {
            cycle_d(U,d,fd,dold,H,gc,R,U->maxlevel,lc,mgp_d,M,myid,nbprocs,t,fd_norm);
            if (myid==0) std::cout<<std::endl;
          } 
       //   Dirichlet_residual_d(U,R,U->maxlevel,M);      
          r_norm=grid_norm(R.getLevel(U->maxlevel),M);
          sresidual=r_norm/fd_norm;
          if(myid==0)
              std::cout<<"Rd="<<r_norm<<",fd="<<fd_norm<<",Cycle="<<j<<std::endl;
          // display global residual by the mmaster MPI
          if (sresidual<=1.e-12) 
          {
            if(myid==0)
              std::cout<<" Final Relative Residual = "<<sresidual<<std::endl;
            break;
          }
        }
        else 
        {
          if(myid==0)
              std::cout<<" Converged "<<std::endl;
            break;
        }
      }
}
