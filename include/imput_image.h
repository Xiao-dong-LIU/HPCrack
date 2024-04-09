/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _IMPUT_Image_H
#define _IMPUT_Image_H
#include "grid2d.h"
#include "iostream"
#include "grid.h"
#include "mg.h"
#include "mpi_struct.h"
using namespace std;


void input_setting(inputdomain & IDM, vector <int> const& element_nb, const int maxlevel);
///------------ raed gray level from image 
void input(grid <short>  & coef3d, const inputdomain & IDM, const MPI_Setting & M, 
 const string & filename);
///-------------- obtain gray level for each MPI processor 
void GL_MPI(Stack *U, grid <short> & m_coef3d, const MPI_Setting & M, const inputdomain & IDM, 
			const string & inputfilename);
//--------------obtain MP on the finest grid
void image_MP(Stack *U, grid<double> & bulK, grid<double> & G, grid<double> & gc, const MPI_Setting & M, 
				const inputdomain & IDM, const string & inputfilename, const double voxel_size);
// ------------- calculate Vogit ane Reuss MP on coarse grids
void voight (grid <double> & MPc, const grid <double> & MPf, const MPI_Setting & M);
void reuss (grid <double> & MPc, const grid <double> & MPf, const MPI_Setting & M);
// ------------- calculate homogenized MP on coarse grids
void coarse_VR(mg<double> & MP, MPI_Setting & M);
// ------------- obtain MP on all grids
void input_total(Stack *U, mg<double> & bulK, mg<double> & G, mg<double> & gc, 
MPI_Setting & M, const inputdomain & IDM, const string & inputfilename,const double voxel_size);
#endif
