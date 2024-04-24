/*=========================================================================

Copyright (c) 2024 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#ifndef CSV_raed_h
#define CSV_raed_h

#include <iostream>
#include <vector>
 
using namespace std;
///funvtion to read CSV table  
void material_property( vector<int> & graylevel, vector<vector<double>> & material_property);

int findPosition(const std::vector<int>& graylevel, int x);

#endif  