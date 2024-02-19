
<a name="readme-top"></a>

# HPCrack


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#description">Description</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#conditions-of-use ">Conditions of use </a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>


## Description
HPCrack is a massively parallel code for the simulation of fracture problems directly from 3D tomographic images. It is a hybrid MPI/OpenMP library developped on C/C++.


## Getting Started
### Prerequisites
MPI and C++ std17 library compiler are required for the compilation.
### Installation
You can clone the project to your local directory 
```
git clone https://src.koda.cnrs.fr/xiaodong/hpcrack.git
```
## Compilation 
You can just do a make with the provided Makefile
```
make
```

## Contributing
If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement". Don't forget to give the project a star! Thanks again!
1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request


## Conditions of use 
Copyright CNRS, École Centrale de Nantes, NANTES Université,

HPCrack is provided to you free of charge. It is released under the CeCILL-B license (see LICENCE, and https://cecill.info/licences/Licence_CeCILL-B_V1-fr.html). Except for src/vtkBase64Utilities.cpp and include/vtkBase64Utilities.h, which are distributed under the BSD license.

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited
liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or 
data to be ensured and, more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL-B license and that you accept its terms.

You can acknowledge (using references [1]) the contribution of this library in any scientific publication dependent upon the use of the library.

[1] X Liu, J Réthoré, AA Lubrecht, An efficient matrix-free preconditioned conjugate gradient based multigrid method for phase field modeling of fracture in heterogeneous materials from 3D images, Computer Methods in Applied Mechanics and Engineering, 388, 114266

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Contact

## Acknowledgments
Lilou Gauthier is thanked for the generation of the image in the example. 
<p align="right">(<a href="#readme-top">back to top</a>)</p>

