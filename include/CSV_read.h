/*=========================================================================
 * 
 *  file CSV_read.h
 *  Created on Mon Feb 12 2024 by
 *  Xiaodong LIU (xiaodong.liu@cnrs.fr)
 *  Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183
 * 
 * 
 *  Copyright (c) 2024 CNRS, Ecole Centrale de Nantes, Nantes Université
 * 
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