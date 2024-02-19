/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "output_clean.h"
#include <dirent.h>
#include <iostream>
#include <stdio.h>

void clean_output_folder(void)
{
    // These are data types defined in the "dirent" header
    DIR *theFolder = opendir("Output");
    struct dirent *next_file;
	std::string filepath_string;
    char *filepath;

    while ( (next_file = readdir(theFolder)) != NULL )
    {
        // build the path for each file in the folder
		filepath_string=std::string("Output/")+std::string(next_file->d_name);
		// std::cout<<filepath_string<<std::endl;
		filepath=&filepath_string[0];
        remove(filepath);
    }
    closedir(theFolder);
}
