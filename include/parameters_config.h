/*=========================================================================

Copyright (c) 2024 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "structure_df.h"

using namespace std;

struct Parameters {
    // Parameters
    int D;
    vector<int> np;
    vector<double> X_start;
    string IMGname;
    vector<int> element_nb;
    double voxel_size;
    double lc; 
    MG mg_u;
    MG mg_d;
    inputdomain IDM;
    int outK,outG,outgc,outdesplacement,outH;
    int outTdata,Tdatafrequence;

    // Function to read parameters from a file
    void readFromFile(const string& filename);

private:
    // Helper function to parse vector values
    void parseVector(const string& str, vector<int>& result);
    void parseVector(const string& str, vector<double>& result);
};


#endif
