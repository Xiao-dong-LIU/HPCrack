/*=========================================================================

Copyright (c) 2024 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include <filesystem> // for C++17 and later
#include "folder_verification.h"

namespace fs = std::filesystem;

void folder_verification() {
    std::string folderPath = "Output";

    // Check if the folder already exists
    if (!fs::exists(folderPath)) {
        // If it doesn't exist, create it
        if (fs::create_directory(folderPath)) {
            std::cout << "Output created successfully!" << std::endl;
        } else {
            std::cerr << "Failed to create Output!" << std::endl;
        }
    } else {
        std::cout << "Output already exists!" << std::endl;
    }
}
