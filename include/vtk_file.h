/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef vtk_file_h
#define vtk_file_h
#include <iostream>
#include <string>

using namespace std;

// extension_para : 0 for vti, 1 for pvti
// l :  level of output 
// myid : proc ID
// t : time steps
string get_file_name (const string & field_name, const int extension_para, 
const int l, const int myid, const int t);

// Open output file
FILE * open_file(const string & filename);
// Close output file
void close_file(FILE * fp);

#endif