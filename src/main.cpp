/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "PF_schema.h"
#include "omp.h"
#include "mpi.h"
#include "mpi_env_set.h"
#include "initialisation.h"
#include "imput_image.h"
#include "PF_schema.h"
#include "mg_cycle_parameter.h"
#include "output_result.h"
#include "mg.h"
#include "mpi_free.h"
#include "parameters_config.h"


using namespace std;

int main(int argc, char * argv[]){


    ///------------ MPI Initialization ------------///
    int nbprocs, myid;
    int level_mpi_provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &level_mpi_provided);
	MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	Parameters para;
	para.readFromFile("config/config.txt");

	///------------ Verification of MPI topology
	int np_prod=1;
	for (int d=0;d<para.np.size();d++)
		{
			np_prod*=para.np[d];
		}
	if(np_prod!=nbprocs)
	{
		if (myid==0)
			std::cerr<<"Wrong NB of MPI tasks"<<std::endl;
		exit(-1);
	}
	///------------ Creation of the Cartesian topology
	MPI_Setting  M;
	Set_MPI(&M,para.np);
	///------------ Compute maxlevel
	// max MG levels, started by 0 
	int maxlevel = 0;
	// NB of element on the finest grid
    int nb_ele_f = para.IDM.voxelnx - para.IDM.coefxa -1;
	// NB of element on coarse grid
	int nb_ele_c = para.element_nb[0];
	if (2*nb_ele_c > nb_ele_f)
		maxlevel = 0;
	else
	{
		do 
		{
			nb_ele_c = 2*nb_ele_c;
			if (nb_ele_c <= nb_ele_f)
			maxlevel += 1; 
		}while (nb_ele_c <= nb_ele_f);
	}
	///------------ Input setting ------------///
	input_setting(para.IDM,  para.element_nb, maxlevel);

	///------------ Initialization  ------------///
    Stack  U;
    initialize(&U, &para.IDM, maxlevel, para.element_nb, para.X_start ,M);
    ///------------ Image and material setting ------------///
	mg<double> bulK(&U,maxlevel);   /// bulK modulus
	mg<double> G(&U,maxlevel);		/// shear modulus G
	mg<double> gc(&U,maxlevel);		/// gc
	///------------ Assigne material property 
	input_total(&U,bulK,G,gc,M,para.IDM,para.IMGname,para.voxel_size);

	///------------ Output material property
	for (int j=maxlevel;j<=maxlevel;j++)
	{
		write_total_de(&U,bulK.getLevel(j),M,"S","bulK","bulK",j,myid,nbprocs,0);
		write_total_de(&U,gc.getLevel(j),M,"S","gc","gc",j,myid,nbprocs,0);
		write_total_de(&U,G.getLevel(j),M,"S","G","G",j,myid,nbprocs,0);
	}
	/// ----------- lc of phase field 
	Level *L;
	L=U.Ll+maxlevel;
	double lc=L->hz*3.0;

 	double t1 = MPI_Wtime();      /// get start time
	phm(&U,bulK,G,gc,M,para.mg_u,para.mg_d,lc,myid,nbprocs,para.voxel_size);
	/// get final time
	double t2 = MPI_Wtime();
	double t = t2-t1;
	if(myid==0) cout << "I am processor " <<myid<< ", Time = "<<t<<endl;
	mpi_comm_free(M);
	finalize(&U);
	MPI_Finalize();
	return 0;
}
