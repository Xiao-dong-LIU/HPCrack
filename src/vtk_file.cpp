/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#include "vtk_file.h"

using namespace std;


string get_file_name (const string & field_name, const int extension_para, 
const int l, const int myid, const int t)
// extension_para : 0 for vti, 1 for pvti
// l :  level of output 
// myid : proc ID
// t : time steps
{
    string filename;
    string espace = "_";
    string folder = "Output/";
    filename = folder + field_name + espace + "l" + espace;
    filename += to_string(l);
    filename += (espace+"t"+espace+to_string(t));
    if (extension_para==0)
    {
        filename += (espace + "Proc" + espace);
        filename += to_string(myid);
        filename += ".vti";
    }
    else if (extension_para==1)
    {
        filename += ".pvti";
    }
    return filename;
}


FILE * open_file(const string & filename)
{
    FILE * fp;
    fp = fopen(filename.c_str(), "wb+");
    return fp;
}

void close_file(FILE * fp)
{
    fclose(fp);
    fp = NULL;
}

