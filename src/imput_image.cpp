/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "imput_image.h" 
#include "mpi_free.h"
#include "omp.h"
#include "realdomain.h"
#include "mpi_communications.h"
#include "csv_read.h"


void input_setting(inputdomain & IDM, vector <int> const& element_nb, const int maxlevel)
/// Function to define the imput domian of image, x0,y0,z0 are the start voxel position
/// nx0, ny0, nz0 are the number of elements on the coarest grid 
/// _voxelnx, _voxelny, _voxelnz NB of voxels in image
{
	IDM.coefxb = IDM.coefxa+element_nb[0]*pow(2,maxlevel);			///end of domain in X, start by 0
    IDM.coefyb = IDM.coefya+element_nb[1]*pow(2,maxlevel);			///end of domain in Y, start by 0
    IDM.coefzb = IDM.coefza+element_nb[2]*pow(2,maxlevel);			///end of domain in Z, start by 0
}


///------------ raed gray level from image ------------///
void input(grid <short>  & coef3d, const inputdomain & IDM, const MPI_Setting & M, 
 const string & filename)
{
	// c++ style
	 std::ifstream is(filename, std::ios::binary | std::ios::in);
	if(!is) 
    {
      if(M.coordinates[0]==0&&M.coordinates[1]==0&&M.coordinates[2]==0)
	 	std::cerr << "Cannot open file!" << std::endl;
		 exit(-1);
    }
    else
    {
		// get length of file:
		is.seekg (0, is.end);
		long long length = is.tellg();
		is.seekg (0, is.beg);
		unsigned char* buffer = new unsigned char[length];
		if(M.coordinates[0]==0&&M.coordinates[1]==0&&M.coordinates[2]==0)
		std::cout << "Reading " << length << " characters...\n";
		std::string read;
		is.read(reinterpret_cast< char *> (buffer),length);
		if (is)
		{
		if(M.coordinates[0]==0&&M.coordinates[1]==0&&M.coordinates[2]==0)
			std::cout << "all characters read successfully.\n";
		}
		else
		{
			if(M.coordinates[0]==0&&M.coordinates[1]==0&&M.coordinates[2]==0)
			std::cout << "error: only " << is.gcount() << " could be read\n";
		}
		// ...buffer contains the entire file...
		
		for(size_t i=0; i< coef3d.axis()[0]; i++)
			for(size_t j=0; j< coef3d.axis()[1]; j++)
				for(size_t k=0; k< coef3d.axis()[2]; k++)
				{
					long long  x = IDM.coefxa+i+IDM.voxelnx*(IDM.coefya+j)+IDM.voxelnx*IDM.voxelny*(IDM.coefza+k);
					coef3d(i,j,k) = buffer[x];
				}
		delete[] buffer;
  }
}

///-------------- obtain gray level for each MPI processor 
void GL_MPI(Stack *U, grid <short> & m_coef3d, const MPI_Setting & M, const inputdomain & IDM, 
			const string & inputfilename)
{
    Level *L;
	L = U->Ll+U->maxlevel;	 
//	cout << IDM.coefxb-IDM.coefxa+1 <<" "<<IDM.coefyb-IDM.coefya+1<< " "<<IDM.coefzb-IDM.coefza+1<<endl;
	vector<size_t> imagesize = {(size_t)(IDM.coefxb-IDM.coefxa+1), 
				    			(size_t)(IDM.coefyb-IDM.coefya+1), 
				     			(size_t)(IDM.coefzb-IDM.coefza+1)};
	grid <short> coef3d(imagesize);
//	cout << coef3d.axis()[0] <<" "<<coef3d.axis()[1]<< " "<<coef3d.axis()[2]<<endl;
	input(coef3d,IDM,M,inputfilename);
	//	exit (-1);

	long x,y,z;
	
	#pragma omp parallel for private(x,y,z)
	for(size_t k =0;k<m_coef3d.axis()[2];k++)
		for(size_t j =0;j<m_coef3d.axis()[1];j++)
			for(size_t i =0;i<m_coef3d.axis()[0];i++)
	    {
			if(M.coordinates[0]==0) {x = i;}
			else if (M.coordinates[0] == M.dims[0]-1)
				{x = L->ii-L->nbcolox+1+i;}
			else
				{x = M.coordinates[0]*(L->nbcolox-2)+i-1;}
			if(M.coordinates[1]==0) {y = j;}
			else if (M.coordinates[1] == M.dims[1]-1)
				{y = L->jj-L->nbcoloy+1+j;}
			else
				{y = M.coordinates[1]*(L->nbcoloy-2)+j-1;}
			if(M.coordinates[2]==0) {z = k;}
			else if (M.coordinates[2] == M.dims[2]-1)
				{z = L->kk-L->nbcoloz+1+k;}
			else
				{z = M.coordinates[2]*(L->nbcoloz-2)+k-1;}
			m_coef3d(i,j,k)=coef3d(x,y,z);
		}	
	
}
//--------------obtain MP on the finest grid
void image_MP(Stack *U, grid<double> & bulK, grid<double> & G, grid<double> & gc, const MPI_Setting & M, 
const inputdomain & IDM, const string & inputfilename, const double voxel_size)
{
	
	grid <short> m_coef(bulK.axis());
	GL_MPI(U,m_coef,M,IDM,inputfilename);
	vector<int> graylevel;
	vector<vector<double>> MP; 
    material_property(graylevel,MP);
	#pragma omp parallel for
	for(size_t k =0;k<bulK.axis()[2];k++)
	  	for(size_t j =0;j<bulK.axis()[1];j++)
			for(size_t i =0;i<bulK.axis()[0];i++)
			{	
				int position = findPosition(graylevel,m_coef(i,j,k));
				double E  = MP[0][position] + (m_coef(i,j,k)-graylevel[position])/(graylevel[position+1]-graylevel[position])*(MP[0][position+1]-MP[0][position]);
				double nu = MP[1][position] + (m_coef(i,j,k)-graylevel[position])/(graylevel[position+1]-graylevel[position])*(MP[1][position+1]-MP[1][position]);
				gc(i,j,k) = MP[2][position] + (m_coef(i,j,k)-graylevel[position])/(graylevel[position+1]-graylevel[position])*(MP[2][position+1]-MP[2][position]);
				E = voxel_size*E;
				bulK(i,j,k)=1.0/(3.0*(1.0-2.0*nu))*E;
				G(i,j,k)=1.0/(2.0*(1.0+nu))*E;
			}
}

void voight (grid <double> & MPc, const grid <double> & MPf, const MPI_Setting & M)
{
	vector <size_t> realD;
	realD = realdomain(MPc.axis(),M);
	int x,y,z;
	int Xlimit = MPf.axis()[0] - 1;
	int Ylimit = MPf.axis()[1] - 1;
	int Zlimit = MPf.axis()[2] - 1;
	#pragma omp parallel for private(x,y,z)
	for (int ic=realD[0]; ic<realD[1]; ic++)
	  for (int jc=realD[2]; jc<realD[3]; jc++)
		for(int kc=realD[4]; kc<realD[5]; kc++)
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
			MPc(ic,jc,kc)=1.0/27.0*( MPf(max(x-1,0), 	   max(y-1,0), 		max(z-1,0))+
								     MPf(max(x-1,0), 	   max(y-1,0), 		min(z+1,Zlimit))+
									  MPf(max(x-1,0),	   min(y+1,Ylimit), max(z-1,0))+
									  MPf(max(x-1,0),	   min(y+1,Ylimit), min(z+1,Zlimit))+
									  MPf(min(x+1,Xlimit), max(y-1,0), 		max(z-1,0))+
									  MPf(min(x+1,Xlimit), max(y-1,0), 		min(z+1,Zlimit))+
									  MPf(min(x+1,Xlimit), min(y+1,Ylimit), max(z-1,0))+
									  MPf(min(x+1,Xlimit), min(y+1,Ylimit), min(z+1,Zlimit))+
							          MPf(max(x-1,0), 	   max(y-1,0), 		z  )+
							          MPf(max(x-1,0), 	   min(y+1,Ylimit), z  )+
									  MPf(min(x+1,Xlimit), max(y-1,0), 		z  )+
									  MPf(min(x+1,Xlimit), min(y+1,Ylimit), z  )+
									  MPf(max(x-1,0), 	   y,   			max(z-1,0))+
									  MPf(max(x-1,0), 	   y,   			min(z+1,Zlimit))+
									  MPf(min(x+1,Xlimit), y,   			max(z-1,0))+
									  MPf(min(x+1,Xlimit), y,   			min(z+1,Zlimit))+
									  MPf(x,   			   max(y-1,0), 		max(z-1,0))+
									  MPf(x,   			   max(y-1,0),		min(z+1,Zlimit))+
									  MPf(x,   			   min(y+1,Ylimit), max(z-1,0))+
									  MPf(x,   			   min(y+1,Ylimit), min(z+1,Zlimit))+
							          MPf(x,   			   y,   			max(z-1,0))+
							          MPf(x,   			   y,   			min(z+1,Zlimit))+
									  MPf(x,   			   max(y-1,0), 		z  )+
									  MPf(x,  			   min(y+1,Ylimit), z  )+
									  MPf(max(x-1,0), 	   y,   			z  )+
									  MPf(min(x+1,Xlimit), y,   			z  )+
							          MPf(x, y, z));
		}
}


void reuss (grid <double> & MPc, const grid <double> & MPf, const MPI_Setting & M)
{
	vector <size_t> realD;
	realD = realdomain(MPc.axis(),M);
	int x,y,z;
	int Xlimit = MPf.axis()[0] - 1;
	int Ylimit = MPf.axis()[1] - 1;
	int Zlimit = MPf.axis()[2] - 1;
	#pragma omp parallel for private(x,y,z)
	for (int ic=realD[0]; ic<realD[1]; ic++)
	  for (int jc=realD[2]; jc<realD[3]; jc++)
		for(int kc=realD[4]; kc<realD[5]; kc++)
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
			MPc(ic,jc,kc)=27.0/(1./ MPf(max(x-1,0), 	   max(y-1,0), 		max(z-1,0))+
								1./ MPf(max(x-1,0), 	   max(y-1,0), 		min(z+1,Zlimit))+
								1./ MPf(max(x-1,0),	   min(y+1,Ylimit), max(z-1,0))+
								1./ MPf(max(x-1,0),	   min(y+1,Ylimit), min(z+1,Zlimit))+
								1./ MPf(min(x+1,Xlimit), max(y-1,0), 		max(z-1,0))+
								1./ MPf(min(x+1,Xlimit), max(y-1,0), 		min(z+1,Zlimit))+
								1./ MPf(min(x+1,Xlimit), min(y+1,Ylimit), max(z-1,0))+
								1./ MPf(min(x+1,Xlimit), min(y+1,Ylimit), min(z+1,Zlimit))+
							    1./ MPf(max(x-1,0), 	   max(y-1,0), 		z  )+
							    1./ MPf(max(x-1,0), 	   min(y+1,Ylimit), z  )+
								1./ MPf(min(x+1,Xlimit), max(y-1,0), 		z  )+
								1./ MPf(min(x+1,Xlimit), min(y+1,Ylimit), z  )+
								1./ MPf(max(x-1,0), 	   y,   			max(z-1,0))+
								1./ MPf(max(x-1,0), 	   y,   			min(z+1,Zlimit))+
								1./ MPf(min(x+1,Xlimit), y,   			max(z-1,0))+
								1./ MPf(min(x+1,Xlimit), y,   			min(z+1,Zlimit))+
								1./ MPf(x,   			   max(y-1,0), 		max(z-1,0))+
								1./ MPf(x,   			   max(y-1,0),		min(z+1,Zlimit))+
								1./ MPf(x,   			   min(y+1,Ylimit), max(z-1,0))+
								1./ MPf(x,   			   min(y+1,Ylimit), min(z+1,Zlimit))+
							    1./ MPf(x,   			   y,   			max(z-1,0))+
							    1./ MPf(x,   			   y,   			min(z+1,Zlimit))+
								1./ MPf(x,   			   max(y-1,0), 		z  )+
								1./ MPf(x,  			   min(y+1,Ylimit), z  )+
								1./ MPf(max(x-1,0), 	   y,   			z  )+
								1./ MPf(min(x+1,Xlimit), y,   			z  )+
							    1./ MPf(x, y, z));
		}
}

//------ compute gral level on coarse girds
void coarse_VR(mg<double> & MP, MPI_Setting & M)
{
	///--------- obtain Voigt for all coarse level 
	for (int l = MP.nbl()-1; l>0; l--)
	{  
        grid <double> MPv(MP.getLevel(l-1).axis());
        grid <double> MPr(MPv);
    //    grid <double> Fvr(Fv);
        voight(MPv,MP.getLevel(l),M);
        reuss(MPr,MP.getLevel(l),M);
   //     Fvr = Fv*Fr;
        MP.getLevel(l-1) = 0.5*(0.5*(MPv+MPr)+(MPv*MPr).gridsqrt());
		derived_datatypes(MP.getLevel(l-1).axis(),&M);
        SendRecv(MP.getLevel(l-1),M);
	//	MPI_Barrier(MPI_COMM_WORLD);
		mpi_data_free(M);
    }
}

// ------------- obtain MP on all grids
void input_total(Stack *U, mg<double> & bulK, mg<double> & G, mg<double> & gc, MPI_Setting & M, Parameters & para)
{
	image_MP(U,bulK.getLevel(U->maxlevel),G.getLevel(U->maxlevel),gc.getLevel(U->maxlevel),M,para.IDM,para.IMGname,para.voxel_size);
	coarse_VR(bulK,M);
	coarse_VR(G,M);
	coarse_VR(gc,M);
}
