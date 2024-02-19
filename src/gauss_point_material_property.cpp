/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "gauss_point_material_property.h"
#include "epsilon.h"
#include "imput_image.h"
#include "strain_history.h"
#include "shape_function.h"

//---------------Computer material property at gauss integratiion 
//        point from nodal material property
double gauss_point_MP(const grid <double> & MP, const grid2d & N, const int i, const int j, 
const int k, const int g)
//  MP: material property, i,j,k: node number, l: level number, g: gauss point nmber
{
	int Cx[8] = { 0,1,1,0,0,1,1,0 };
	int Cy[8] = { 0,0,1,1,0,0,1,1 };
	int Cz[8] = { 0,0,0,0,1,1,1,1 };
	double coefg=0.0, tmp;
	int a = 1;
	for (int n = 0; n < 8; n++)
	{
		//   skip element which has a '0' property at its node
		if (MP(Cx[n] + i, Cy[n] + j, Cz[n] + k) == 0.0)
			a = 0;
	}
	if (a == 1)
	{
		tmp = 0.0;
		for (int n = 0; n < 8; n++)  /// Interpolation with shape function
			tmp += MP(Cx[n] + i, Cy[n] + j, Cz[n] + k) * N(g, n);
		coefg = tmp;
	}
	return coefg;
}

//------------- compute K_d or G_d at each gauss integration point
double K_G_d_g(const griddouble &d, const griddouble &MP, const grid2d & N, const int i, const int j, const int k, 
const int g, const double small_k)
// set bulK to G, it can compute G_d
{
	double d_g,coefg,MP_g, g_d;
	MP_g=gauss_point_MP(MP,N,i,j,k,g);
	d_g=gauss_point_MP(d,N,i,j,k,g);
	if (d_g<0.)
		d_g=0.;
	else if (d_g>1.)
		d_g=1.;
	g_d=(1.-small_k)*(1-d_g)*(1-d_g)+small_k;
	/*
	if (i==3&&j==3&&k>=13&&k<=19)
	std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",d="<<d(l-1,i,j,k)<<",d_g="<<d_g<<",g_d="<<g_d<<std::endl;
	*/
	coefg=g_d*MP_g;
	return coefg;
}

//-------------function to obtain (1-k)*g(d)+k
double g_d_node(const griddouble &d, const int i, const int j, const int k, const double small_k)
{
	double g_d,g_d_k;
	double c;
	c=d(i,j,k);
	if (c<0.)
		c=0.;
	else if (c>1.)
		c=1.;
	g_d=(1.0-c)*(1.0-c);
	g_d_k=(1.-small_k)*g_d+small_k;
	return(g_d_k);
}

//------------- compute K_d and G_d at each node
double MP_d_node(const griddouble &d, const griddouble &MP, 
const int i, const int j, const int k, const double small_k)
{
	double K_d_n,g_d_k;
	g_d_k=g_d_node(d,i,j,k,small_k);
	K_d_n=g_d_k*MP(i,j,k);
	return(K_d_n);
}

//------------- compute K_d and G_d at nodes of the fineest grid
void MP_d(griddouble & MP_d, const griddouble &d, const griddouble & MP, const double small_k)
{
	#pragma omp parallel for
	for (int k=0;k<d.axis()[2];k++)
		for (int j=0;j<d.axis()[1];j++)
			for (int i=0;i<d.axis()[0];i++)
				MP_d(i,j,k)=MP_d_node(d,MP,i,j,k,small_k);
}
