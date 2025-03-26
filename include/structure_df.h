/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _STRUCTURE_DF_H
#define _STRUCTURE_DF_H

typedef struct
{
	double hx;             /// mesh size x
	double hy;             /// mesh size y
	double hz;			   /// mesh size z
	long    ii;             /// number of nodes x
	long    jj;             /// number of nodes y
	long    kk;             /// number of nodes z
	int    nbcolox;         /// number of nodes x in each processor
	int    nbcoloy;         /// number of nodes y in each processor
	int    nbcoloz;         /// number of nodes z in each processor
} Level;

typedef struct
{
	int    maxlevel;        /// NB of levels
	int    nx0;             /// number of nodes coarsest grid
	int    ny0;             /// number of nodes coarsest grid
	int    nz0;             /// number of nodes coarsest grid
	double xa,xb;           /// begin,end of computational domain
	double ya,yb;           /// begin,end of computational domain
	double za,zb;           /// begin,end of computational domain
	double wu;              /// work unit counter
	Level  *Ll ;            /// array of grid levels
} Stack;

//imput for image
typedef struct
{
 int coefxa;            ///begin of domain in X, start by 0
 int coefxb;			///end of domain in X, start by 0
 int coefya;			///begin of domain in Y, start by 0
 int coefyb;			///end of domain in Y, start by 0
 int coefza;			///begin of domain in Z, start by 0
 int coefzb;			///end of domain in Z, start by 0
 int voxelnx;            ///NB of voxels in image on X
 int voxelny;				///NB of voxels in image on Y
 int voxelnz;				///NB of voxels in image on Z
} inputdomain;

//Material property
typedef struct
{
	 double E_C; 				/// Young modulus of inclusion
	 double E_F; 				/// Young modulus of matrix
	 double nuC; 				/// Poisson ratio  of inclusion
	 double nuF; 				/// Poisson ratio of matrix
	 double gc_C;				/// gc of inclusion
	 double gc_F;				/// gc of matrix
	 double gc_interface;
} Material_property;

typedef struct
{
	int nu0;             /// mesh size x
	int nu1;             /// mesh size y
	int nu2;			   /// mesh size z
	int gamma;             /// number of nodes x
	int ncy;             /// number of nodes y
} MG;

#endif
