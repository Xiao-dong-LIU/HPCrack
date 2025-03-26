/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "grid2d.h"
#include <set>


//For grid2d
grid2d::grid2d(int _nx, int _ny):nx(_nx), ny(_ny), a(0){
    ///std::fill(val.begin(), val.end(), 0.);
    a = new double[nx*ny];
    std::fill(a, a+nx*ny, 0.); 
}
grid2d::grid2d(const grid2d &in):nx(in.nx), ny(in.ny), a(0){
	a = new double[nx*ny];
	std::copy(in.a, in.a+nx*ny, a); 
}
  double & grid2d::operator()(int i, int j)
  {
	return a[i+j*nx];
  }
  const double & grid2d::operator()(int i, int j) const
  {
	return a[i+j*nx];
  }
grid2d::~grid2d(){
	delete []a;
	a=0;
}

//For free memory of U
void finalize(Stack *U)
{
/// free memory at end of program
free(U->Ll);
}
