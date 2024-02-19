/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <vector>

using namespace std;


///------------ Problen dimension
    static const int D=3; 
///------------  NB of processor in each direction
    std::vector<int> np={1,1,1};    /// product of composants must equal to nbprocs  
///------------  beginning voxel position
    std::vector<double> X_start = {-64.,-64.,-64.};
///------------ Image file name 
    string IMGname = "data/inclusion_more_129_129_129.raw";  
///------------ NB of voxels in each direction in image
    std::vector<int>  voxel_nb={129,129,129};
///------------ NB of elements on the coarest grid 
    std::vector <int> element_nb={4,4,4};
///------------ ROI setting (image) ------------///
    std::vector<int>  ROI_start={0,0,0};     /// start voxel position 
///------------ voxel size ------------///
    static const double voxel_size=2.53e-6;          /// in m

///------------ Multigrid cycle parameter
	// to compute u (mechanical part)
 	static const int ncy_u =12;       /// NB of mg Cycles
  	static const int nu0_u =50;      /// NB of relaxation of the coarest grid
	static const int nu1_u=4;        /// NB of relaxxation going up
  	static const int nu2_u=2;        /// NB of relaxxation going down
  	static const int gamma_u=1;      /// Cycle type, 1 for V Cycle, 2 for W Cycle
    // to compute d (phase field part)
	static const int ncy_d =12;       /// NB of mg Cycles
  	static const int nu0_d =10;      /// NB of relaxation of the coarest grid
	static const int nu1_d=2;        /// NB of relaxxation going up
  	static const int nu2_d=1;        /// NB of relaxxation going down
  	static const int gamma_d=1;      /// Cycle type, 1 for V Cycle, 2 for W Cycle
#endif
