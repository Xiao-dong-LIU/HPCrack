/*=========================================================================

Copyright (c) 2024 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include <filesystem> // for C++17 and later
#include "folder_verification.h"

namespace fs = std::filesystem;

void folder_verification(const string & folderpath) {

    // Check if the folder already exists
    if (!fs::exists(folderpath)) {
        // If it doesn't exist, create it
        if (fs::create_directory(folderpath)) {
            std::cout <<folderpath<< " created successfully!" << std::endl;
        } else {
            std::cerr << "Failed to create " <<folderpath<<"!"<< std::endl;
        }
    } else {
        std::cout <<folderpath<< " already exists!" << std::endl;
    }
}
