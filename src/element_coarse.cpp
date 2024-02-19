/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "element_coarse.h"
#include "omp.h"
#include <cmath>
#include "mpi_free.h"
#include "mpi_communications.h"

using namespace std;


// ------- compute voigt average on coarse gird
void coarse_voigt(griddouble & Hc, const griddouble & Hf, const MPI_Setting & M)
{
    int sx=1,sy=1,sz=1,ex=Hc.axis()[0]-2,ey=Hc.axis()[1]-2,ez=Hc.axis()[2]-2;
	if (M.coordinates[0]==0) sx =0;
    if (M.coordinates[1]==0) sy =0;
    if (M.coordinates[2]==0) sz =0;
	int x,y,z;
    #pragma omp parallel for private(x,y,z)
	for (int ic=sx;ic<=ex;ic++)
	  for (int jc=sy;jc<=ey;jc++)
		for(int kc=sz;kc<=ez;kc++)
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
            Hc(ic,jc,kc)=0.125*(Hf(x,y  ,z  )+Hf(x+1,y  ,z  )+
                                Hf(x,y+1,z  )+Hf(x+1,y+1,z  )+
                                Hf(x,y  ,z+1)+Hf(x+1,y  ,z+1)+
                                Hf(x,y+1,z+1)+Hf(x+1,y+1,z+1));
        }
}


// ------- compute reuss average on coarse gird
void coarse_reuss(griddouble& Hc, const griddouble& Hf, const MPI_Setting & M)
{
    int sx=1,sy=1,sz=1,ex=Hc.axis()[0]-2,ey=Hc.axis()[1]-2,ez=Hc.axis()[2]-2;
	if (M.coordinates[0]==0) sx =0;
    if (M.coordinates[1]==0) sy =0;
    if (M.coordinates[2]==0) sz =0;
	int x,y,z;
    #pragma omp parallel for private(x,y,z)
	for (int ic=sx;ic<=ex;ic++)
	  for (int jc=sy;jc<=ey;jc++)
		for(int kc=sz;kc<=ez;kc++)
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
            Hc(ic,jc,kc)=8./(1./Hf(x,y,z)+1./Hf(x+1,y,z)+
                             1./Hf(x,y+1,z)+1./Hf(x+1,y+1,z)+
                             1./Hf(x,y,z+1)+1./Hf(x+1,y,z+1)+
                             1./Hf(x,y+1,z+1)+1./Hf(x+1,y+1,z+1));
        }
}


//------ compute voigt-reuss type element material property  on coarse girds
void VR_all_level(mgdouble & H, MPI_Setting & M)
{
    for (int l = H.nbl()-1; l>0; l--)
	{  
        grid <double> Hv(H.getLevel(l-1).axis());
        grid <double> Hr(Hv);
    //    grid <double> Fvr(Fv);
        coarse_voigt(Hv,H.getLevel(l),M);
        coarse_reuss(Hr,H.getLevel(l),M);
   //     Fvr = Fv*Fr;
        H.getLevel(l-1) = 0.5*(0.5*(Hv+Hr)+(Hv*Hr).gridsqrt());
		derived_datatypes(H.getLevel(l-1).axis(),&M);
        SendRecv(H.getLevel(l-1),M);
   //     MPI_Barrier(MPI_COMM_WORLD);
        mpi_data_free(M);
    }
}


//------ compute voigt type element material property  on coarse girds
void Voigt_all_level(mgdouble & H, MPI_Setting & M)
{
    for (int l = H.nbl()-1; l>0; l--)
	{  
        grid <double> Hv(H.getLevel(l-1).axis());
        grid <double> Hr(Hv);
        coarse_voigt(H.getLevel(l-1),H.getLevel(l),M);
		derived_datatypes(H.getLevel(l-1).axis(),&M);
        SendRecv(H.getLevel(l-1),M);
     //   MPI_Barrier(MPI_COMM_WORLD);
        mpi_data_free(M);
    }
}
