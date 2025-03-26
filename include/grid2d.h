/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _GRID2D_H
#define _GRID2D_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "structure_df.h"


// 2d vector class
class grid2d{
	
private:
  int nx, ny;
  ///std::vector<double > val;
  double *a;
public:
	/// Construct a dense nx*ny matrix
  grid2d(const int nx, const int ny);
	/// Copy constructor
  grid2d(const grid2d & in);
   /// Operator ()
  double & operator()(int i, int j);
  const double & operator()(int i, int j) const;
  ~grid2d();
};

// free memory at end of program
void finalize(Stack *U);

#endif
