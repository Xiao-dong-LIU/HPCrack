/*=========================================================================

Copyright (c) 2024 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "parameters_config.h"
#include <sstream> // Include for stringstream


void Parameters::readFromFile(const string& filename) {
    ifstream file(filename);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            // Ignore lines starting with '#'
            if (line.empty() || line[0] == '#') {
                continue;
            }
            size_t delimiter_pos = line.find('=');
            if (delimiter_pos != string::npos) {
                string key = line.substr(0, delimiter_pos);
                string value = line.substr(delimiter_pos + 1);
                // Parse and set parameter values
                if (key == "D") {
                    D = stoi(value);
                } else if (key == "np") {
                    parseVector(value, np);
                } else if (key == "X_start") {
                    parseVector(value, X_start);
                } else if (key == "IMGname") {
                    IMGname = value;
                } else if (key == "voxel_nb_X") {
                    IDM.voxelnx = stod(value);
                } else if (key == "voxel_nb_Y") {
                    IDM.voxelny = stod(value);
                } else if (key == "voxel_nb_Z") {
                    IDM.voxelnz = stod(value);
                } else if (key == "ROI_start_X") {
                    IDM.coefxa = stod(value);
                } else if (key == "ROI_start_Y") {
                    IDM.coefya = stod(value);
                } else if (key == "ROI_start_Z") {
                    IDM.coefza = stod(value);
                } else if (key == "voxel_size") {
                    voxel_size = stod(value);
                } else if (key == "element_nb") {
                     parseVector(value, element_nb);
                } else if (key == "ncy_u") {
                    mg_u.ncy = stoi(value);
                } else if (key == "nu0_u") {
                    mg_u.nu0 = stoi(value);
                } else if (key == "nu1_u") {
                    mg_u.nu1 = stoi(value);
                } else if (key == "nu2_u") {
                    mg_u.nu2 = stoi(value);
                } else if (key == "gamma_u") {
                    mg_u.gamma = stoi(value);
                } else if (key == "ncy_d") {
                    mg_d.ncy = stoi(value);
                } else if (key == "nu0_d") {
                    mg_d.nu0 = stoi(value);
                } else if (key == "nu1_d") {
                    mg_d.nu1 = stoi(value);
                } else if (key == "nu2_d") {
                    mg_d.nu2 = stoi(value);
                } else if (key == "gamma_d") {
                    mg_d.gamma = stoi(value);
                } else if (key == "lc") {
                    lc = stoi(value);
                } else if (key == "outK") {
                    outK = stoi(value);
                } else if (key == "outG") {
                    outG = stoi(value);
                } else if (key == "outgc") {
                    outgc = stoi(value);
                } else if (key == "outdesplacement") {
                    outdesplacement = stoi(value);
                } else if (key == "outH") {
                    outH = stoi(value);
                } else if (key == "outTdata") {
                    outTdata = stoi(value);
                } else if (key == "Tdatafrequence") {
                    Tdatafrequence = stoi(value);
                }
                
                // Add parsing for other parameters here
            }
        }
        file.close();
    } else {
        cerr << "Error: Unable to open file: " << filename << endl;
    }
}

void Parameters::parseVector(const string& str, vector<int>& result) {
    stringstream ss(str);
    string token;
    while (getline(ss, token, ',')) {
        result.push_back(stoi(token));
    }
}

void Parameters::parseVector(const string& str, vector<double>& result) {
    stringstream ss(str);
    string token;
    while (getline(ss, token, ',')) {
        result.push_back(stod(token));
    }
}

