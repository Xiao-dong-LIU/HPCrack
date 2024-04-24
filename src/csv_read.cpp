/*=========================================================================

Copyright (c) 2024 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#include <fstream>
#include <string>
#include <sstream>
#include "csv_read.h"
 
void material_property( vector<int> & graylevel, vector<vector<double>> & material_property)
{
    vector<vector<string>> content;
    vector<string> row;
    string line, word;
    //-----------------  read the CSV table
    fstream file ("data/Material_Property.csv", ios::in);
    if(file.is_open())
    {
        while(getline(file, line))
        {
            row.clear();
 
            stringstream str(line);
 
            while(getline(str, word, ','))
                row.push_back(word);
            content.push_back(row);
        }
    }
    else
    {
        std::cerr<<"Could not open the file\n";
         exit(EXIT_FAILURE);
    }
    //----------------- cast gray level of different phase to int 
    for (int i=1; i<content[0].size(); i++)
    {
        graylevel.push_back(stoi(content[0][i]));
    }
    vector<double> mpone;
    //----------------- cast material property of different phase to double 
    for(int i=1; i<content.size(); i++)
    {
        mpone.clear();
        for(int j=1; j<content[i].size(); j++)
        {
           mpone.push_back(std::stod(content[i][j]));
        }
        material_property.push_back(mpone);
    }
/*
    //----------------- Print for un check 
    for (int i=0; i<graylevel.size(); i++)
        cout << " "<<graylevel[i]<<" ";
    cout<<endl; 
    for (int i=0; i<material_property.size(); i++)
    {
        for (int j=0; j<material_property[i].size(); j++)
        {
            cout <<material_property[i][j]<<" ";
        }
        cout<<endl;
    }
*/
}


int findPosition(const std::vector<int>& graylevel, int x) 
{
    for (int i = 0; i < graylevel.size() - 1; ++i) {
        if (x >= graylevel[i] && x <= graylevel[i + 1]) {
            return i; // Return the position where x fits
        }
    }
    return -1; // If x is not within the range of any consecutive elements
}

 