/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "gauss_point_material_property.h"
#include "shape_function.h"
#include "material_choice.h"
#include "strain_history.h"
#include "output_result.h"
#include "epsilon.h"
#include "omp.h"
#include "mpi.h"

void MP_choice_element(gridshort & MP_choice, const gdvecdouble & u, const griddouble & bulK, 
const griddouble & G, const griddouble & d,  const griddouble & DN, const grid2d & N, 
int l, const int i, const int j, const int k)
{
	double K_g,G_g;   // bulK and G at gauss point
	double phi_positive=0.,phi_negative=0.,tr_epsilon=0.;
	grid2d epsilon(3, 3);
	for (int g=0;g<8;g++)    /// Gauss Points
		{
			//--- get material property at gausse point 
			K_g = gauss_point_MP(bulK, N, i, j, k, g);
			G_g = gauss_point_MP(G, N, i, j, k, g);
			//--- compute epsilon at this gausse point
			epsilon_routine(epsilon, u, DN, i, j, k, g);
			tr_epsilon += trace_epsilon(epsilon);
			//---phi+ at gauss point
			phi_positive += Phi_positive(K_g, G_g, epsilon);
			//---phi- at Gauss point
			phi_negative += Phi_negative(K_g, epsilon);
		}
   	if(tr_epsilon>=0.0)
		MP_choice(i,j,k)=0;
	if (tr_epsilon<0.0) 
       	MP_choice(i,j,k)=1;
	if (phi_positive<phi_negative)
		MP_choice(i,j,k)=2;
}

void MP_choice(Stack *U, mgshort & MP_choice_d, const mgvecdouble & u, const mgdouble & bulK, 
const mgdouble & G, const mgdouble & d)
{
    Level* L;
    // define anray shape function and its gradient
	vector<size_t> DNsize = {3,8,8};
    grid <double> DN(DNsize);
	grid2d N(8, 8);
    for (int l=U->maxlevel;l>=0;l--)
    {
        L=U->Ll+l;
        Dshape_function(N, DN, L->hx, L->hy, L->hz);
        #pragma omp parallel for                
        for (int k = 0; k <= d.getLevel(l).axis()[2] - 2; k++)   /// Nb. of Elements
			for (int i = 0; i <= d.getLevel(l).axis()[1] - 2; i++)
				for (int j = 0; j <= d.getLevel(l).axis()[0] - 2; j++)
                    MP_choice_element(MP_choice_d.getLevel(l),u.getLevel(l),bulK.getLevel(l),G.getLevel(l),d.getLevel(l),DN,N,l,i,j,k);
    }
}
/*
void account_nb(Stack *U, multiinput3d& MP_choice_d, MPI_Setting & M)
{
    Level* L;
    for (int l=U->maxlevel;l>=1;l--)
    {
        L=U->Ll+l;
		int a=0,b=0,c=0;
		int sa=0,sb=0,sc=0;
    //    #pragma omp parallel for
        for (int i = 0; i <= L->nbcolox - 2; i++)
            for (int j = 0; j <= L->nbcoloy - 2; j++)
                for (int k = 0; k <= L->nbcoloz - 2; k++)   /// Nb. of Elements
				{
					if (MP_choice_d(l-1,i,j,k)==0)
						a=a+1;
					if (MP_choice_d(l-1,i,j,k)==1)
						b=b+1;
					if (MP_choice_d(l-1,i,j,k)==2)
						c=c+1;
				}
	MPI_Reduce(&a, &sa, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sa, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Reduce(&b, &sb, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sb, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Reduce(&c, &sc, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sc, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (M.coordinates[0]==0&&M.coordinates[1]==0&&M.coordinates[2]==0)
		std::cout<<"l="<<l<<",a="<<sa<<",b="<<sb<<",c="<<sc<<std::endl;
    }
}
*/
void MP_element(griddouble & K_element, griddouble & G_element, const gdvecdouble & u, 
const griddouble & bulK, const griddouble & G, const griddouble & d, const gridshort & MP_choice_d, 
const griddouble & DN, const grid2d & N, int i, int j, int k, double small_k)
{
	double K_g,G_g;   // bulK and G at gauss point
	double K_e=0.,G_e=0.;   // bulK and G in element
	for (int g=0;g<8;g++)    /// Gauss Points
		{
			//--- get material property at gausse point 
			K_g = gauss_point_MP(bulK, N, i, j, k, g);
			G_g = gauss_point_MP(G, N, i, j, k, g);
			//--- get K_d and G_d for different cases
			// for K
			if (MP_choice_d(i,j,k)==0) 
			{
					K_e+=K_G_d_g(d,bulK,N,i,j,k,g,small_k);
			}
			else 
			K_e+=K_g;
			// for G
			if (MP_choice_d(i,j,k)==2)
			{
				G_e+=G_g;
			}
			else 
			{
				G_e+=K_G_d_g(d,G,N,i,j,k,g,small_k);
			}
		}
	K_element(i,j,k)=0.125*K_e;
	G_element(i,j,k)=0.125*G_e;
}

void MP_all_element(Stack*U, const gdvecdouble & u,  const griddouble & bulK, 
const griddouble & G, const griddouble & d, const gridshort & MP_choice_d, griddouble & K_element, 
griddouble & G_element, const int l, const double small_k)
{
	Level* L;
	L=U->Ll+l;
    // define anray shape function and its gradient
	vector<size_t> DNsize = {3,8,8};
    grid <double> DN(DNsize);
	grid2d N(8, 8);
    Dshape_function(N, DN, L->hx, L->hy, L->hz);
    //#pragma omp parallel for
    for (int i = 0; i <= L->nbcolox - 2; i++)
       for (int j = 0; j <= L->nbcoloy - 2; j++)
           for (int k = 0; k <= L->nbcoloz - 2; k++)   /// Nb. of Elements
		        MP_element(K_element,G_element,u,bulK,G,d,MP_choice_d,DN,N,i,j,k,small_k);
}
/*
void account_nb_zeros(Stack *U, multigrid3d& MP, MPI_Setting M)
{
    Level* L;
    for (int l=U->maxlevel;l>=1;l--)
    {
        L=U->Ll+l;
		int a=0;
		int sa=0;
    //    #pragma omp parallel for
        for (int i = 0; i <= L->nbcolox - 2; i++)
            for (int j = 0; j <= L->nbcoloy - 2; j++)
                for (int k = 0; k <= L->nbcoloz - 2; k++)   /// Nb. of Elements
				{
					if (MP(l-1,i,j,k)<1.e-12)
						a=a+1;
				}
	MPI_Reduce(&a, &sa, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sa, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (M.coordinates[0]==0&&M.coordinates[1]==0&&M.coordinates[2]==0)
		std::cout<<"l="<<l<<",a="<<sa<<std::endl;
    }
}
*/