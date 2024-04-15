#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream> // For file I/O
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
    int outK,outG,outgc;

    // Function to read parameters from a file
    void readFromFile(const string& filename);

private:
    // Helper function to parse vector values
    void parseVector(const string& str, vector<int>& result);
    void parseVector(const string& str, vector<double>& result);
};


#endif
